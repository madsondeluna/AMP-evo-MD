#!/usr/bin/env python3
"""
Extract per-candidate, trajectory-derived series and save them for the per-candidate figures.

The scalar observables for every candidate are already in observables.json and are aggregated by
amp_dataset.py. The per-candidate detail figures (P1-P13) need quantities that only exist if the
trajectory was kept: density profiles, tail order, tilt, contacts, interaction-energy time
series, area-per-lipid and thickness maps, and the binding PMF. Those are expensive to compute,
so this module computes each one ONCE and writes it to <dataset>/<id>/, where amp_report.py
reads it. A figure is therefore never recomputed from a trajectory; re-rendering reads the CSV.

Trajectories are only present for candidates run with WRITE_XTC=True (off by default; the run
refuses to start a full campaign with it on). The usual workflow is to pick a shortlist from the
cohort figures, re-run just those few with WRITE_XTC=True, and point this module at them.

WHAT IS WRITTEN (under <dataset>/<id>/)
---------------------------------------
  density.csv       z, headgroups, tails, water, peptide           (P1)
  apl_map.txt       2D grid, area per lipid                         (P2)
  thickness_map.txt 2D grid, bilayer thickness                      (P3)
  order.csv         bead, far, under                                (P4)
  stress.csv        z, pi                                           (P5, needs a stress run)
  depth.csv         time, d_mem                                     (P6, else distance.xvg is used)
  tilt.csv          time, tilt                                      (P7)
  contacts.txt      2D grid, residue x frame contact occupancy      (P8)
  comz.csv          time, z                                         (P9)
  energy.csv        time, e_lipids, e_solvent                       (P10)
  convergence.csv   time, cum, se                                   (P11)
  pmf.csv           xi, target, host (+ _err)                       (P13, from a WHAM profile)

DEPENDENCIES
------------
Geometry series use MDAnalysis; energy series use the gmx binary. Both are optional: each
extractor is guarded and logs why it was skipped, so a partial environment still produces
whatever it can. This module imports nothing from the run itself.

USAGE
-----
    python tools/amp_traj.py <campaign_dir> --candidate <id> --dataset <campaign>/analysis \\
        --membrane bacterial [--production production]
"""

import argparse
import csv
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mdtools

# Martini bead-name sets. These match the common Martini 3 POPE/POPG/POPC definitions; adjust
# here if a membrane uses different bead names. Kept as module constants so the selection logic
# is in one place and not scattered through the extractors.
HEADGROUP_BEADS = ("PO4", "NH3", "NC3", "GL0", "CNO")
TAIL_BEAD_PREFIX = ("C", "D")          # aliphatic tail beads: C1A, D2A, C3B ...
WATER_BEADS = ("W", "WF")
PEPTIDE_SEL = "name BB or name SC1 or name SC2 or name SC3 or name SC4 or name SC5"


def _write_csv(path, columns):
    """columns is an ordered dict {name: 1d array}; all arrays share length."""
    names = list(columns)
    n = len(columns[names[0]])
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for i in range(n):
            writer.writerow([columns[k][i] for k in names])


def _have_mdanalysis():
    try:
        import MDAnalysis  # noqa: F401
        return True
    except ImportError:
        return False


def _gmx():
    return os.environ.get("GMX", "gmx")


# --- topology / trajectory discovery -----------------------------------------------------------
def _find_inputs(candidate_dir, membrane, production):
    """Locate the trajectory, topology and energy file for one membrane of one candidate."""
    base = "%s_%s" % (production, membrane)
    xtc = os.path.join(candidate_dir, base + ".xtc")
    edr = os.path.join(candidate_dir, base + ".edr")
    tpr = os.path.join(candidate_dir, base + ".tpr")
    gro = os.path.join(candidate_dir, base + ".gro")
    topology = tpr if os.path.isfile(tpr) else (gro if os.path.isfile(gro) else None)
    return {"xtc": xtc if os.path.isfile(xtc) else None,
            "edr": edr if os.path.isfile(edr) else None,
            "topology": topology}


# --- geometry extractors (MDAnalysis) ----------------------------------------------------------
def extract_density(u, outdir, nbins=200):
    """z-density of headgroups, tails, water and peptide, membrane-centred."""
    import MDAnalysis as mda  # noqa: F401
    box_z = u.dimensions[2] / 10.0  # Angstrom -> nm
    edges = np.linspace(0, box_z, nbins + 1)
    z = 0.5 * (edges[:-1] + edges[1:])
    sel = {
        "headgroups": u.select_atoms("name " + " ".join(HEADGROUP_BEADS)),
        "water": u.select_atoms("name " + " ".join(WATER_BEADS)),
        "peptide": u.select_atoms(PEPTIDE_SEL),
    }
    tails = u.select_atoms("name " + " ".join("%s*" % p for p in TAIL_BEAD_PREFIX))
    sel["tails"] = tails
    hist = {k: np.zeros(nbins) for k in sel}
    nframes = 0
    for _ in u.trajectory:
        mid = sel["headgroups"].center_of_geometry()[2] / 10.0 if len(sel["headgroups"]) else box_z / 2
        for k, group in sel.items():
            if len(group) == 0:
                continue
            zc = group.positions[:, 2] / 10.0
            zc = (zc - mid + box_z / 2) % box_z  # recentre so the bilayer sits mid-box
            h, _ = np.histogram(zc, bins=edges)
            hist[k] += h
        nframes += 1
    cols = {"z": z}
    for k in ("headgroups", "tails", "water", "peptide"):
        v = hist[k] / max(nframes, 1)
        m = v.max()
        cols[k] = v / m if m > 0 else v  # normalise each to its own peak (a.u.)
    _write_csv(os.path.join(outdir, "density.csv"), cols)
    return "density.csv"


def extract_tilt_comz(u, outdir):
    """Peptide tilt (principal axis vs z) and peptide COM relative to the membrane midplane."""
    pep = u.select_atoms(PEPTIDE_SEL)
    head = u.select_atoms("name " + " ".join(HEADGROUP_BEADS))
    if len(pep) == 0:
        return []
    times, tilts, comz = [], [], []
    for ts in u.trajectory:
        pos = pep.positions - pep.center_of_geometry()
        # Principal axis = eigenvector of the largest eigenvalue of the gyration tensor.
        cov = np.cov(pos.T)
        w, v = np.linalg.eigh(cov)
        axis = v[:, np.argmax(w)]
        angle = np.degrees(np.arccos(abs(axis[2]) / (np.linalg.norm(axis) + 1e-9)))
        times.append(ts.time / 1000.0)  # ps -> ns
        tilts.append(angle)
        mid = head.center_of_geometry()[2] if len(head) else 0.0
        comz.append((pep.center_of_mass()[2] - mid) / 10.0)
    _write_csv(os.path.join(outdir, "tilt.csv"), {"time": np.array(times), "tilt": np.array(tilts)})
    _write_csv(os.path.join(outdir, "comz.csv"), {"time": np.array(times), "z": np.array(comz)})
    return ["tilt.csv", "comz.csv"]


def extract_contacts(u, outdir, cutoff=6.0, max_frames=160):
    """Per-residue contact occupancy with lipids over time (residue x frame grid)."""
    from MDAnalysis.analysis import distances
    pep = u.select_atoms(PEPTIDE_SEL)
    lip = u.select_atoms("name " + " ".join(HEADGROUP_BEADS) +
                         " or name " + " ".join("%s*" % p for p in TAIL_BEAD_PREFIX))
    if len(pep) == 0 or len(lip) == 0:
        return []
    residues = pep.residues
    frames = list(u.trajectory)
    stride = max(1, len(frames) // max_frames)
    sampled = frames[::stride]
    occ = np.zeros((len(residues), len(sampled)))
    for j, ts in enumerate(sampled):
        u.trajectory[ts.frame]
        for i, res in enumerate(residues):
            ratoms = res.atoms.select_atoms(PEPTIDE_SEL)
            if len(ratoms) == 0:
                continue
            d = distances.distance_array(ratoms.positions, lip.positions, box=u.dimensions)
            occ[i, j] = 1.0 if d.min() < cutoff else 0.0
    np.savetxt(os.path.join(outdir, "contacts.txt"), occ, fmt="%.3f")
    return ["contacts.txt"]


def extract_maps(u, outdir, nbins=24):
    """Time-averaged area-per-lipid and bilayer-thickness maps over the xy plane."""
    head = u.select_atoms("name PO4")
    if len(head) == 0:
        return []
    box = u.dimensions[:2] / 10.0
    edges_x = np.linspace(0, box[0], nbins + 1)
    edges_y = np.linspace(0, box[1], nbins + 1)
    count = np.zeros((nbins, nbins)); thick = np.zeros((nbins, nbins)); hits = np.zeros((nbins, nbins))
    nframes = 0
    for _ in u.trajectory:
        pos = head.positions / 10.0
        zc = pos[:, 2]
        upper = zc > zc.mean()
        hx = np.clip(np.digitize(pos[:, 0], edges_x) - 1, 0, nbins - 1)
        hy = np.clip(np.digitize(pos[:, 1], edges_y) - 1, 0, nbins - 1)
        for k in range(len(pos)):
            count[hy[k], hx[k]] += 1
        # thickness per bin: mean upper-leaflet z minus mean lower-leaflet z
        for bi in range(nbins):
            for bj in range(nbins):
                inbin = (hx == bj) & (hy == bi)
                if inbin.any():
                    up = zc[inbin & upper]; lo = zc[inbin & ~upper]
                    if up.size and lo.size:
                        thick[bi, bj] += up.mean() - lo.mean(); hits[bi, bj] += 1
        nframes += 1
    cell_area = (box[0] / nbins) * (box[1] / nbins)
    with np.errstate(divide="ignore", invalid="ignore"):
        apl = cell_area / (count / max(nframes, 1))
        apl[~np.isfinite(apl)] = np.nan
        thickness = np.where(hits > 0, thick / np.maximum(hits, 1), np.nan)
    np.savetxt(os.path.join(outdir, "apl_map.txt"), np.nan_to_num(apl, nan=np.nanmean(apl)), fmt="%.4f")
    np.savetxt(os.path.join(outdir, "thickness_map.txt"),
               np.nan_to_num(thickness, nan=np.nanmean(thickness)), fmt="%.4f")
    return ["apl_map.txt", "thickness_map.txt"]


def extract_order(u, outdir, radius=12.0):
    """
    Second-rank order parameter per tail-bond index, split by distance to the peptide.

    S = (3<cos^2 theta> - 1)/2, theta the angle of each consecutive tail bond to z. "under" are
    lipids whose headgroup is within `radius` (Angstrom) of the peptide COM in xy, "far" the rest.
    """
    pep = u.select_atoms(PEPTIDE_SEL)
    # Group tail beads per lipid residue in name order; bond i is bead i -> bead i+1.
    lip_res = u.select_atoms("name " + " ".join("%s*" % p for p in TAIL_BEAD_PREFIX)).residues
    if len(lip_res) == 0:
        return []
    max_bonds = 0
    per_res_beads = []
    for res in lip_res:
        beads = res.atoms.select_atoms("name " + " ".join("%s*" % p for p in TAIL_BEAD_PREFIX))
        per_res_beads.append(beads)
        max_bonds = max(max_bonds, len(beads) - 1)
    acc = {"far": [np.zeros(0) for _ in range(max_bonds)],
           "under": [np.zeros(0) for _ in range(max_bonds)]}
    sums = {g: np.zeros(max_bonds) for g in acc}
    counts = {g: np.zeros(max_bonds) for g in acc}
    for _ in u.trajectory:
        pcom = pep.center_of_geometry()[:2] if len(pep) else np.array([0.0, 0.0])
        for beads in per_res_beads:
            p = beads.positions
            if len(p) < 2:
                continue
            near = np.linalg.norm(beads.center_of_geometry()[:2] - pcom) < radius
            grp = "under" if near else "far"
            for b in range(len(p) - 1):
                vec = p[b + 1] - p[b]
                cz = vec[2] / (np.linalg.norm(vec) + 1e-9)
                sums[grp][b] += cz * cz
                counts[grp][b] += 1
    bead = np.arange(1, max_bonds + 1)
    cols = {"bead": bead}
    for grp in ("far", "under"):
        with np.errstate(invalid="ignore"):
            mean_cos2 = np.where(counts[grp] > 0, sums[grp] / np.maximum(counts[grp], 1), np.nan)
        cols[grp] = (3 * mean_cos2 - 1) / 2
    _write_csv(os.path.join(outdir, "order.csv"), cols)
    return ["order.csv"]


# --- energy extractors (gmx) -------------------------------------------------------------------
def _gmx_energy(edr, terms, out_xvg):
    """Run gmx energy for the named terms; return (time, {term: series}) or None."""
    selection = ("\n".join(terms) + "\n").encode()
    try:
        subprocess.run([_gmx(), "energy", "-f", edr, "-o", out_xvg],
                       input=selection, capture_output=True, timeout=300, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    if not os.path.isfile(out_xvg):
        return None
    x, series, _ = mdtools.read_xvg(out_xvg)
    return x, series


def extract_energy(edr, outdir):
    """Peptide-lipid and peptide-solvent interaction-energy time series from the .edr."""
    res = _gmx_energy(edr, ["Coul-SR:Protein-Lipids", "LJ-SR:Protein-Lipids",
                            "Coul-SR:Protein-Solvent", "LJ-SR:Protein-Solvent"],
                      os.path.join(outdir, "_energy.xvg"))
    if res is None:
        return []
    x, series = res
    vals = list(series.values())
    # Sum Coul-SR + LJ-SR for each interaction (the gmx selection order is as requested above).
    e_lip = vals[0] + vals[1] if len(vals) >= 2 else vals[0]
    e_sol = vals[2] + vals[3] if len(vals) >= 4 else None
    cols = {"time": x / 1000.0, "e_lipids": e_lip}
    if e_sol is not None:
        cols["e_solvent"] = e_sol
    _write_csv(os.path.join(outdir, "energy.csv"), cols)
    return ["energy.csv"]


def extract_convergence(edr, outdir, area_nm2=None, surften_free=None):
    """
    Running mean of delta_gamma over the trajectory, from the per-frame surface tension.

    delta_gamma(t) = surften_free - gamma(t), with gamma(t) read from the .edr. The running mean
    and its block standard error show whether the 400 ns window is long enough for the candidate's
    damage estimate to have settled.
    """
    res = _gmx_energy(edr, ["#Surf*SurfTen"], os.path.join(outdir, "_surften.xvg"))
    if res is None or surften_free is None:
        return []
    x, series = res
    gamma_bar_nm = next(iter(series.values()))
    # Convert bar*nm to kJ mol-1 nm-2 the same way the pipeline does: * 0.06022 per leaflet pair.
    gamma = gamma_bar_nm * 0.06022140857 / 2.0
    dg = surften_free - gamma
    cum = np.cumsum(dg) / np.arange(1, len(dg) + 1)
    # Block standard error at each point, five blocks over the data seen so far.
    se = np.zeros_like(cum)
    for i in range(1, len(dg)):
        blocks = mdtools.block_average(dg[:i + 1], n_blocks=min(5, i + 1))
        se[i] = np.std(blocks) / np.sqrt(len(blocks)) if len(blocks) > 1 else 0.0
    _write_csv(os.path.join(outdir, "convergence.csv"),
               {"time": x / 1000.0, "cum": cum, "se": se})
    return ["convergence.csv"]


def extract_depth_from_xvg(candidate_dir, membrane, production, outdir):
    """depth.csv from the distance xvg the run already kept, so P6 needs no re-run."""
    xvg = os.path.join(candidate_dir, "%s_%s_distance.xvg" % (production, membrane))
    if not os.path.isfile(xvg):
        return []
    x, series, _ = mdtools.read_xvg(xvg)
    _write_csv(os.path.join(outdir, "depth.csv"),
               {"time": x / 1000.0, "d_mem": next(iter(series.values()))})
    return ["depth.csv"]


def extract_pmf(candidate_dir, outdir):
    """pmf.csv from a WHAM profile if the candidate was run with umbrella sampling."""
    found = {}
    for role, membrane in (("target", "bacterial"), ("host", "host")):
        for name in ("profile_%s.xvg" % membrane, "pmf_%s.xvg" % membrane, "%s_profile.xvg" % membrane):
            path = os.path.join(candidate_dir, name)
            if os.path.isfile(path):
                x, series, _ = mdtools.read_xvg(path)
                found[role] = (x, next(iter(series.values())))
                break
    if not found:
        return []
    xi = next(iter(found.values()))[0]
    cols = {"xi": xi}
    for role, (_x, y) in found.items():
        cols[role] = y
    _write_csv(os.path.join(outdir, "pmf.csv"), cols)
    return ["pmf.csv"]


# --- driver ------------------------------------------------------------------------------------
def run(campaign_dir, candidate, dataset_dir, membrane="bacterial", production="production"):
    candidate_dir = os.path.join(campaign_dir, candidate)
    if not os.path.isdir(candidate_dir):
        sys.exit("No candidate directory: %s" % candidate_dir)
    outdir = os.path.join(dataset_dir, candidate)
    os.makedirs(outdir, exist_ok=True)
    inputs = _find_inputs(candidate_dir, membrane, production)
    written = []

    # Always-available series kept by the run, no trajectory needed.
    written += extract_depth_from_xvg(candidate_dir, membrane, production, outdir)
    written += extract_pmf(candidate_dir, outdir)

    # Energy / convergence from the .edr.
    if inputs["edr"]:
        written += extract_energy(inputs["edr"], outdir)
        # surften_free for this membrane is recoverable from the saved cohort; left to the caller
        # via metadata if needed. Convergence is attempted only when surften_free is supplied.
    else:
        sys.stderr.write("  [skip] energy/convergence: no .edr for %s\n" % membrane)

    # Geometry series from the trajectory.
    if inputs["xtc"] and inputs["topology"] and _have_mdanalysis():
        import MDAnalysis as mda
        u = mda.Universe(inputs["topology"], inputs["xtc"])
        for fn in (extract_density, extract_tilt_comz, extract_contacts, extract_maps, extract_order):
            try:
                out = fn(u, outdir)
                written += out if isinstance(out, list) else [out]
            except Exception as exc:
                sys.stderr.write("  [skip] %s: %s\n" % (fn.__name__, exc))
    else:
        reason = "no trajectory" if not inputs["xtc"] else \
                 ("no topology" if not inputs["topology"] else "MDAnalysis not installed")
        sys.stderr.write("  [skip] geometry series (density/order/tilt/contacts/maps): %s\n" % reason)

    # Remove the gmx intermediates; the series CSVs are the kept artifacts.
    for tmp in ("_energy.xvg", "_surften.xvg"):
        p = os.path.join(outdir, tmp)
        if os.path.isfile(p):
            try:
                os.remove(p)
            except OSError:
                pass

    written = [w for w in written if w]
    print("Wrote %d series for %s to %s: %s"
          % (len(written), candidate, outdir, ", ".join(sorted(written)) or "(none)"))
    return written


def main(argv=None):
    parser = argparse.ArgumentParser(description="Extract per-candidate trajectory series for the P figures.")
    parser.add_argument("campaign_dir", help="GROMACS run output directory")
    parser.add_argument("--candidate", required=True, help="candidate id (sequence / directory name)")
    parser.add_argument("--dataset", required=True, help="analysis dir written by amp_dataset.py")
    parser.add_argument("--membrane", default="bacterial", help="membrane type to analyse (default bacterial)")
    parser.add_argument("--production", default="production", help="production file prefix (UserSettings.name_output_production)")
    args = parser.parse_args(argv)
    run(args.campaign_dir, args.candidate, args.dataset,
        membrane=args.membrane, production=args.production)


if __name__ == "__main__":
    main()
