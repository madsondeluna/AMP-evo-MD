#!/usr/bin/env python3
"""
Aggregate a finished or running campaign into a single, self-describing analysis dataset.

WHY THIS EXISTS
---------------
Every candidate already writes an observables.json into its own directory, and the genetic
algorithm appends one row per evaluation to the progress file. Those are the raw records, but
they are scattered across one directory per sequence and one line per evaluation. Answering any
cohort-level question ("which residues separate bacterial from host damage across 100 AMPs?")
means reading all of them and joining the sequence metadata from the progress file onto the
observables. This module does that join once and writes the result to a small set of files that
every figure then reads. Nothing downstream re-reads a raw observables.json or re-parses the
progress file, and nothing is ever recomputed from a trajectory to redraw a figure.

WHAT IS WRITTEN (under <campaign>/analysis by default)
------------------------------------------------------
  cohort.csv       one row per candidate: every scalar observable for every membrane role,
                   plus sequence descriptors (charge, hydrophobic moment, group fractions) and
                   the GA metadata (generation, fitness, raw score, orientation, n_run), plus
                   the derived damage / host-damage / selectivity.
  sequences.csv    one row per candidate: the sequence and its per-residue counts over the
                   active alphabet, for the amino-acid analyses.
  metadata.json    campaign-level provenance and a data dictionary: GROMACS/martinize2 build,
                   loaded modules, the membrane-type to role map, which candidates have a kept
                   trajectory, the selectivity definition, and the column descriptions.

The three files together are the dataset. They are plain text, diffable, and portable off the
cluster, which is the whole point: the measurements and the metadata travel together.

USAGE
-----
    python tools/amp_dataset.py <campaign_dir> [--progress output.txt] [--outdir analysis]

<campaign_dir> is the GROMACS run output directory (core.settings.output_dir), the one that
contains one subdirectory per evaluated sequence and the progress file.
"""

import argparse
import csv
import datetime
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mdtools
import ampstyle

try:
    import ga_report
except ImportError:
    ga_report = None

OBSERVABLES_FILENAME = "observables.json"
# Active alphabet: the residues the generator is allowed to use. Kept here as the analysis-side
# copy so the dataset is readable without importing the run settings; validated against the data.
ACTIVE_ALPHABET = list("AEFKLMQRSWY")
PIPELINE_VERSION = "1.0"


# --- small helpers -----------------------------------------------------------------------------
def _scalar(entry, key):
    """Pull a scalar mean out of an observable field that may be {mean,err,units}, a float, or None."""
    value = entry.get(key)
    if isinstance(value, dict):
        return value.get("mean")
    return value


def _err(entry, key):
    value = entry.get(key)
    if isinstance(value, dict):
        return value.get("err")
    return None


def _load_observables(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def _progress_metadata(progress_path):
    """
    Map each sequence to its GA metadata: first generation seen, last fitness, raw, orientation,
    max n_run. Returns ({sequence: {...}}, settings). Falls back to empty if the file or the
    parser is unavailable, so the dataset can still be built from observables alone.
    """
    per_seq = {}
    settings = {}
    if ga_report is None or not os.path.isfile(progress_path):
        return per_seq, settings, {}
    settings, iterations, _unparsed = ga_report.parse_progress(progress_path)
    for gen in sorted(iterations):
        for rec in iterations[gen]:
            seq = rec["sequence"]
            if seq not in per_seq:
                # First time we see the sequence fixes the generation it was introduced.
                per_seq[seq] = {"generation": gen, "fitness": rec["fitness"],
                                "raw": rec["raw"], "orientation": rec["orientation"],
                                "n_run": rec["n_run"]}
            else:
                # Later appearances carry the running-average fitness and the higher n_run.
                per_seq[seq]["fitness"] = rec["fitness"]
                per_seq[seq]["raw"] = rec["raw"]
                per_seq[seq]["orientation"] = rec["orientation"]
                per_seq[seq]["n_run"] = max(per_seq[seq]["n_run"], rec["n_run"])
    return per_seq, settings, iterations


def _ga_timeline(iterations):
    """
    Per-iteration monitoring of the search, from the progress file alone.

    Each row summarises one iteration: optimiser performance (best, parents, mean fitness),
    population turnover (evaluated vs unique), diversity (mean pairwise Hamming, mean positional
    entropy) and composition (mean residue-group fractions). This is the saved form of what
    ga_report draws, so the monitoring figures never re-parse the progress file.
    """
    rows = []
    for gen in sorted(iterations):
        recs = iterations[gen]
        if not recs:
            continue
        seqs = [r["sequence"] for r in recs]
        fit = np.array([r["fitness"] for r in recs], float)
        finite = fit[np.isfinite(fit)]
        order = np.argsort(-np.where(np.isfinite(fit), fit, -np.inf))
        n_parents = max(1, len(recs) // 4)  # top quartile stands in for the selected parents
        row = {
            "generation": gen,
            "n_evaluated": len(recs),
            "n_unique": len(set(seqs)),
            "fitness_best": float(np.nanmax(fit)) if finite.size else float("nan"),
            "fitness_mean": float(np.nanmean(fit)) if finite.size else float("nan"),
            "fitness_parents": float(np.nanmean(fit[order[:n_parents]])) if finite.size else float("nan"),
            "pairwise_hamming": mdtools.mean_pairwise_hamming(seqs),
            "entropy_mean": float(np.mean(mdtools.positional_entropy(seqs)))
                            if seqs else float("nan"),
        }
        totals = {g: 0 for g in mdtools.RESIDUE_GROUPS}
        nres = 0
        for s in seqs:
            for g, c in mdtools.group_counts(s).items():
                totals[g] += c
            nres += len(s)
        for g in mdtools.RESIDUE_GROUPS:
            row["frac_%s" % g] = totals[g] / nres if nres else float("nan")
        rows.append(row)
    return rows


def _positional_frequency(iterations):
    """Residue-group frequency at each position, for the final iteration's population."""
    if not iterations:
        return []
    last = max(iterations)
    seqs = [r["sequence"] for r in iterations[last]]
    freqs, length = mdtools.positional_group_frequency(seqs)
    rows = []
    for g, arr in freqs.items():
        for pos in range(length):
            rows.append({"group": g, "position": pos + 1, "frequency": float(arr[pos])})
    return rows


# --- the collector -----------------------------------------------------------------------------
def collect(campaign_dir, progress_path=None, outdir=None):
    """Build the dataset and write it; return the metadata dict describing what was written."""
    campaign_dir = os.path.abspath(campaign_dir)
    if progress_path is None:
        progress_path = os.path.join(campaign_dir, "output.txt")
    if outdir is None:
        outdir = os.path.join(campaign_dir, "analysis")
    os.makedirs(outdir, exist_ok=True)

    per_seq, settings, iterations = _progress_metadata(progress_path)

    # Discover candidate directories: any directory directly under the campaign holding an
    # observables.json. Walking one level deep is enough and avoids descending into failed_* and
    # analysis/ subtrees.
    candidates = []
    for name in sorted(os.listdir(campaign_dir)):
        cdir = os.path.join(campaign_dir, name)
        obs_path = os.path.join(cdir, OBSERVABLES_FILENAME)
        if os.path.isdir(cdir) and os.path.isfile(obs_path):
            candidates.append((name, cdir, obs_path))

    role_of = {}        # membrane_type -> role, recorded for the data dictionary
    roles_seen = set()
    provenance = None
    rows = []
    seq_rows = []
    traj_available = []

    for cid, cdir, obs_path in candidates:
        obs = _load_observables(obs_path)
        if obs is None:
            continue
        if provenance is None and isinstance(obs.get("_provenance"), dict):
            provenance = obs["_provenance"]

        sequence = cid  # candidate directory is the uppercased sequence
        meta = per_seq.get(sequence, {})

        # Sequence descriptors, all from the sequence alone (no simulation needed).
        charge = mdtools.net_charge(sequence)
        groups = mdtools.group_counts(sequence)
        length = len(sequence)
        row = {
            "id": cid,
            "sequence": sequence,
            "length": length,
            "generation": meta.get("generation"),
            "fitness": meta.get("fitness"),
            "raw": meta.get("raw"),
            "orientation_ga": meta.get("orientation"),
            "n_run": meta.get("n_run"),
            "charge": charge,
            "hydrophobic_moment": ampstyle.hydrophobic_moment(sequence),
            "mean_hydrophobicity": ampstyle.mean_hydrophobicity(sequence),
        }
        for g, c in groups.items():
            row["frac_%s" % g] = c / length if length else float("nan")

        # Per-membrane-role observables. Each membrane_type key maps to a role; columns are keyed
        # by role so target/host stay stable even if a membrane is renamed.
        damage = host_damage = None
        for mtype, entry in obs.items():
            if mtype.startswith("_") or not isinstance(entry, dict) or "role" not in entry:
                continue
            role = entry["role"]
            role_of[mtype] = role
            roles_seen.add(role)
            row["dg_%s" % role] = _scalar(entry, "delta_gamma")
            row["dg_err_%s" % role] = _err(entry, "delta_gamma")
            row["dgwork_%s" % role] = _scalar(entry, "delta_gamma_work")
            row["surften_%s" % role] = _scalar(entry, "surften")
            row["dmem_%s" % role] = _scalar(entry, "d_mem")
            row["orient_%s" % role] = entry.get("orientation")
            row["orient_frames_%s" % role] = entry.get("orientation_frames_at_interface")
            row["elip_%s" % role] = _scalar(entry, "e_lipids")
            row["esol_%s" % role] = _scalar(entry, "e_solvent")
            row["area_%s" % role] = entry.get("area_nm2")
            row["bnm_%s" % role] = entry.get("b_nm")
            if role == "target":
                damage = row["dg_%s" % role]
            elif role == "host":
                host_damage = row["dg_%s" % role]

        # Derived selectivity, documented in metadata. Reward damage to the target, penalise
        # damage to the host; host damage below zero (membrane stabilised) is not rewarded.
        row["damage"] = damage
        row["host_damage"] = host_damage
        if damage is not None and host_damage is not None:
            row["selectivity"] = damage - 1.5 * max(host_damage, 0.0)
        else:
            row["selectivity"] = None
        rows.append(row)

        # Per-residue counts for the amino-acid analyses.
        srow = {"id": cid, "sequence": sequence,
                "generation": meta.get("generation"), "fitness": meta.get("fitness"),
                "selectivity": row["selectivity"], "damage": damage, "host_damage": host_damage}
        for aa in ACTIVE_ALPHABET:
            srow["n_%s" % aa] = sequence.count(aa)
        seq_rows.append(srow)

        # Record which candidates kept a trajectory (xtc) or the persisted distance series, so
        # the per-candidate report knows what it can draw without a re-run.
        has_xtc = any(f.endswith(".xtc") for f in os.listdir(cdir))
        has_dist = any(f.endswith("_distance.xvg") for f in os.listdir(cdir))
        if has_xtc or has_dist:
            traj_available.append({"id": cid, "xtc": has_xtc, "distance_xvg": has_dist})

    if not rows:
        sys.exit("No candidates with %s found under %s. Point this at the GROMACS output "
                 "directory (core.settings.output_dir)." % (OBSERVABLES_FILENAME, campaign_dir))

    # --- write cohort.csv ---------------------------------------------------------------------
    # Union of keys across rows, with a stable preferred order first so the file is readable.
    preferred = ["id", "sequence", "length", "generation", "fitness", "raw", "orientation_ga",
                 "n_run", "charge", "hydrophobic_moment", "mean_hydrophobicity",
                 "damage", "host_damage", "selectivity"]
    all_keys = list(preferred)
    for row in rows:
        for k in row:
            if k not in all_keys:
                all_keys.append(k)
    cohort_path = os.path.join(outdir, "cohort.csv")
    _write_csv(cohort_path, rows, all_keys)

    # --- write sequences.csv ------------------------------------------------------------------
    seq_keys = ["id", "sequence", "generation", "fitness", "selectivity", "damage", "host_damage"]
    seq_keys += ["n_%s" % aa for aa in ACTIVE_ALPHABET]
    seq_path = os.path.join(outdir, "sequences.csv")
    _write_csv(seq_path, seq_rows, seq_keys)

    # --- write metadata.json ------------------------------------------------------------------
    metadata = {
        "pipeline_version": PIPELINE_VERSION,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "campaign_dir": campaign_dir,
        "progress_file": progress_path if os.path.isfile(progress_path) else None,
        "n_candidates": len(rows),
        "n_with_generation": sum(1 for r in rows if r.get("generation") is not None),
        "roles": sorted(roles_seen),
        "membrane_role_map": role_of,
        "active_alphabet": ACTIVE_ALPHABET,
        "provenance": provenance or {},
        "ga_settings": settings,
        "trajectories_available": traj_available,
        "n_with_trajectory": sum(1 for t in traj_available if t["xtc"]),
        "selectivity_definition": "selectivity = dg_target - 1.5 * max(dg_host, 0); "
                                  "reward damage to the bacterial target, penalise host damage, "
                                  "do not reward host stabilisation",
        "files": {
            "cohort.csv": "one row per candidate: scalar observables per role + descriptors + GA metadata",
            "sequences.csv": "one row per candidate: sequence and per-residue counts over the active alphabet",
        },
        "columns": {
            "dg_<role>": "delta_gamma (kJ mol-1 nm-2), gamma_free - surften, on the membrane of that role",
            "dgwork_<role>": "delta_gamma * area (kJ mol-1)",
            "dmem_<role>": "peptide-membrane COM distance along z (nm)",
            "orient_<role>": "mean orientation factor over frames (0..1)",
            "elip_<role>": "peptide-lipid interaction energy Coul-SR+LJ-SR (kJ mol-1)",
            "esol_<role>": "peptide-solvent interaction energy Coul-SR+LJ-SR (kJ mol-1)",
            "damage": "dg on the target (bacterial) membrane",
            "host_damage": "dg on the host membrane",
            "selectivity": "see selectivity_definition",
            "hydrophobic_moment": "Eisenberg hydrophobic moment per residue (ideal helix, 100 deg)",
        },
    }
    # --- write the GA timeline and positional frequency (monitoring, from the progress file) ---
    timeline = _ga_timeline(iterations)
    if timeline:
        tl_keys = ["generation", "n_evaluated", "n_unique", "fitness_best", "fitness_parents",
                   "fitness_mean", "pairwise_hamming", "entropy_mean"]
        tl_keys += ["frac_%s" % g for g in mdtools.RESIDUE_GROUPS]
        _write_csv(os.path.join(outdir, "ga_timeline.csv"), timeline, tl_keys)
        metadata["files"]["ga_timeline.csv"] = "one row per GA iteration: fitness, turnover, diversity, composition"
    positional = _positional_frequency(iterations)
    if positional:
        _write_csv(os.path.join(outdir, "positional.csv"), positional, ["group", "position", "frequency"])
        metadata["files"]["positional.csv"] = "residue-group frequency by position, final generation"
    metadata["n_generations"] = len(timeline)

    meta_path = os.path.join(outdir, "metadata.json")
    with open(meta_path, "w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)

    print("Wrote dataset for %d candidates (%d with GA generation, %d with trajectory) to %s"
          % (len(rows), metadata["n_with_generation"], metadata["n_with_trajectory"], outdir))
    print("  cohort.csv     %d rows x %d columns" % (len(rows), len(all_keys)))
    print("  sequences.csv  %d rows" % len(seq_rows))
    print("  metadata.json  roles=%s" % ", ".join(sorted(roles_seen)))
    return metadata


def _write_csv(path, rows, keys):
    """Write dict rows to CSV with a fixed column order; missing keys are left blank."""
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in keys})


# --- dataset reader (used by every report) -----------------------------------------------------
def load_cohort(path):
    """
    Read a cohort.csv (or sequences.csv) into {column: numpy array}.

    Numeric columns become float arrays (blanks become nan); non-numeric columns stay as object
    arrays of strings. This is the single reader every figure uses, so no report re-parses the
    raw campaign.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        cols = reader.fieldnames or []
        data = {c: [] for c in cols}
        for row in reader:
            for c in cols:
                data[c].append(row.get(c, ""))
    out = {}
    for c, values in data.items():
        numeric = []
        ok = True
        for v in values:
            if v == "":
                numeric.append(np.nan)
                continue
            try:
                numeric.append(float(v))
            except ValueError:
                ok = False
                break
        out[c] = np.array(numeric, float) if ok else np.array(values, dtype=object)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description="Aggregate a campaign into an analysis dataset.")
    parser.add_argument("campaign_dir", help="GROMACS run output directory (core.settings.output_dir)")
    parser.add_argument("--progress", default=None, help="progress file (default <campaign>/output.txt)")
    parser.add_argument("--outdir", default=None, help="output directory (default <campaign>/analysis)")
    args = parser.parse_args(argv)
    collect(args.campaign_dir, progress_path=args.progress, outdir=args.outdir)


if __name__ == "__main__":
    main()
