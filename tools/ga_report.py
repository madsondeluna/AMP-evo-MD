#!/usr/bin/env python3
"""
Track a running or finished campaign: optimiser performance, population diversity, composition.

Reads the progress file the genetic algorithm appends to as it runs, so it can be pointed at a
campaign in flight. Writes a per-iteration CSV, a summary table to stdout, and three figures.

    python tools/ga_report.py output/progress.txt --outdir report/

WHY DIVERSITY IS TRACKED ALONGSIDE FITNESS
------------------------------------------
Rising fitness on its own does not distinguish a search that found a motif from one that
collapsed onto a few near-identical survivors and stopped exploring. The two look the same on a
fitness plot and completely different on positional entropy: a real motif shows entropy falling
at the positions that matter while others stay high, whereas premature convergence pulls every
position down together. Reading them side by side is the point of this report.

WHAT THE PROGRESS FILE CONTAINS, AND WHAT IT DOES NOT
----------------------------------------------------
Each evaluated sequence appears once per iteration with its fitness, the raw pre-orientation
score, the orientation factor, and the number of times it has been simulated. The fitness of a
rerun sequence is the running average over its replicas, so the spread between individual
replicas is not recoverable from this file. That spread is the per-candidate noise, and it is
the quantity that decides whether a fitness difference between two candidates means anything.
To see it, read the individual observables.json files instead.
"""

import argparse
import collections
import csv
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mdtools


ITER_RE = re.compile(r"^Iteration\s+(\d+)\s*:")
# "\tSEQUENCE\tfitness\tddF=...\tpenalty=...\tn_run=..."
ROW_RE = re.compile(
    r"^\t([A-Z]+)\t(-?[\d.eE+naN]+)\tddF=(-?[\d.eE+naN]+)\tpenalty=(-?[\d.eE+naN]+)\tn_run=(\d+)"
)


def _float(text):
    try:
        return float(text)
    except ValueError:
        return float("nan")


def parse_progress(path):
    """
    Return (settings, iterations) where iterations is an ordered {index: [record, ...]}.

    Records are dicts with sequence, fitness, raw, orientation and n_run. Rows that do not parse
    are counted and reported rather than silently dropped: a changed output format should be
    visible, not inferred from a suspiciously empty plot.
    """
    settings = {}
    iterations = collections.OrderedDict()
    current = None
    unparsed = 0

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            m = ITER_RE.match(line)
            if m:
                current = int(m.group(1))
                iterations.setdefault(current, [])
                continue
            if line.startswith("\t") and current is not None:
                m = ROW_RE.match(line)
                if m:
                    iterations[current].append({
                        "sequence": m.group(1),
                        "fitness": _float(m.group(2)),
                        "raw": _float(m.group(3)),
                        "orientation": _float(m.group(4)),
                        "n_run": int(m.group(5)),
                    })
                elif line.strip() and not line.strip().endswith(":"):
                    unparsed += 1
                continue
            if line.startswith("\t") and current is None:
                parts = [p for p in line.split("\t") if p.strip()]
                if len(parts) >= 2:
                    settings[parts[0].strip().rstrip(":")] = parts[1].strip()

    iterations = collections.OrderedDict(
        (k, v) for k, v in sorted(iterations.items()) if v
    )
    return settings, iterations, unparsed


def summarise(iterations, n_parents):
    """Per-iteration metrics. n_parents sets how many top sequences count as the parent pool."""
    rows = []
    best_so_far = -float("inf")
    last_improved = None

    for index, records in iterations.items():
        fitness = np.array([r["fitness"] for r in records], dtype=float)
        finite = fitness[np.isfinite(fitness)]
        if finite.size == 0:
            continue
        order = np.argsort(-fitness)
        parents = fitness[order[:min(n_parents, len(records))]]
        sequences = [r["sequence"] for r in records]

        entropy = mdtools.positional_entropy(sequences)
        best = float(np.max(finite))
        if best > best_so_far:
            best_so_far = best
            last_improved = index

        groups = collections.Counter()
        charges = []
        for s in sequences:
            for g, c in mdtools.group_counts(s).items():
                groups[g] += c
            charges.append(mdtools.net_charge(s))

        row = {
            "iteration": index,
            "n_evaluated": len(records),
            "n_unique": len(set(sequences)),
            "fitness_best": best,
            "fitness_parents_mean": float(parents.mean()),
            "fitness_mean": float(finite.mean()),
            "fitness_sd": float(finite.std(ddof=1)) if finite.size > 1 else 0.0,
            "orientation_mean": float(np.nanmean([r["orientation"] for r in records])),
            "n_rerun": sum(1 for r in records if r["n_run"] > 1),
            "hamming_mean": mdtools.mean_pairwise_hamming(sequences),
            "entropy_mean": float(entropy.mean()) if entropy.size else float("nan"),
            "entropy_min": float(entropy.min()) if entropy.size else float("nan"),
            "charge_mean": float(np.mean(charges)) if charges else float("nan"),
        }
        total_residues = sum(groups.values()) or 1
        for g in mdtools.RESIDUE_GROUPS:
            row["frac_" + g] = groups[g] / total_residues
        rows.append(row)

    return rows, last_improved


def top_sequences(iterations, n):
    """The n highest-fitness sequences across the whole campaign, deduplicated."""
    best = {}
    for records in iterations.values():
        for r in records:
            prev = best.get(r["sequence"])
            if prev is None or r["fitness"] > prev["fitness"]:
                best[r["sequence"]] = r
    ranked = sorted(best.values(), key=lambda r: -r["fitness"])
    return ranked[:n]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def figure_performance(rows, outpath, last_improved):
    plt = mdtools.apply_style()
    it = [r["iteration"] for r in rows]

    # Extra right margin on every panel, so direct labels sit outside the data rather than on
    # top of the last point.
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))

    # Fitness. Direct labels at the right end carry identity; no legend box needed.
    ax = axes[0][0]
    series = [
        ("best", "fitness_best", mdtools.CATEGORICAL[0]),
        ("parents", "fitness_parents_mean", mdtools.CATEGORICAL[1]),
        ("mean", "fitness_mean", mdtools.CATEGORICAL[2]),
    ]
    for label, key, colour in series:
        y = [r[key] for r in rows]
        ax.plot(it, y, color=colour)
        mdtools.label_last(ax, it, y, label, colour)
    ax.set_title("Fitness")
    ax.set_xlabel("iteration")
    mdtools.declutter(ax)
    if last_improved is not None and it:
        stalled = it[-1] - last_improved
        if stalled >= 5:
            ax.axvline(last_improved, color=mdtools.GRID, lw=0.8)
            ax.annotate("flat for %d" % stalled, xy=(last_improved, ax.get_ylim()[0]),
                        xytext=(4, 4), textcoords="offset points",
                        fontsize=7, color=mdtools.INK_MUTED)

    # Diversity. Separate panels rather than a second y axis: different units.
    ax = axes[0][1]
    y_eval = [r["n_evaluated"] for r in rows]
    y_uniq = [r["n_unique"] for r in rows]
    ax.plot(it, y_eval, color=mdtools.CATEGORICAL[3])
    ax.plot(it, y_uniq, color=mdtools.CATEGORICAL[0])
    mdtools.label_last(ax, it, y_eval, "evaluated", mdtools.CATEGORICAL[3])
    mdtools.label_last(ax, it, y_uniq, "unique", mdtools.CATEGORICAL[0])
    ax.set_title("Population turnover")
    ax.set_xlabel("iteration")
    # From zero on purpose. Autoscaling crops to the top few percent and turns a difference of
    # three duplicate sequences into a dramatic-looking collapse; the ratio of unique to
    # evaluated only reads correctly against zero.
    ax.set_ylim(0, max(y_eval) * 1.08)
    mdtools.declutter(ax)

    ax = axes[1][0]
    y = [r["hamming_mean"] for r in rows]
    ax.plot(it, y, color=mdtools.CATEGORICAL[0])
    ax.set_title("Mean pairwise distance")
    ax.set_xlabel("iteration")
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1.0])
    mdtools.declutter(ax, ymax_ticks=3)

    ax = axes[1][1]
    y_mean = [r["entropy_mean"] for r in rows]
    y_min = [r["entropy_min"] for r in rows]
    ax.plot(it, y_mean, color=mdtools.CATEGORICAL[0])
    ax.plot(it, y_min, color=mdtools.CATEGORICAL[2])
    mdtools.label_last(ax, it, y_mean, "mean", mdtools.CATEGORICAL[0])
    mdtools.label_last(ax, it, y_min, "most converged", mdtools.CATEGORICAL[2])
    ax.set_title("Positional entropy (bits)")
    ax.set_xlabel("iteration")
    ax.set_ylim(bottom=0)
    mdtools.declutter(ax)

    # Direct labels sit outside the axes, so reserve room for them rather than letting
    # bbox_inches="tight" crop the longest one.
    fig.tight_layout(w_pad=4.0, h_pad=2.5)
    fig.subplots_adjust(right=0.88)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def figure_composition(rows, outpath):
    """Small multiples, one panel per residue group: six series would exceed the palette."""
    plt = mdtools.apply_style()
    it = [r["iteration"] for r in rows]
    groups = list(mdtools.RESIDUE_GROUPS)

    fig, axes = plt.subplots(2, 3, figsize=(9.5, 4.6), sharex=True, sharey=True)
    for ax, group in zip(axes.ravel(), groups):
        ax.plot(it, [r["frac_" + group] for r in rows], color=mdtools.SEQUENTIAL)
        ax.set_title(group)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1.0])
        mdtools.declutter(ax, ymax_ticks=3, xmax_ticks=4)
    for ax in axes[-1]:
        ax.set_xlabel("iteration")
    axes[0][0].set_ylabel("fraction of residues")
    axes[1][0].set_ylabel("fraction of residues")
    fig.suptitle("Residue group composition of the population",
                 x=0.005, ha="left", fontsize=9.5, color=mdtools.INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2.5, h_pad=2.0)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def figure_positional(iterations, outpath, top_n):
    """
    Positional residue-group frequency among the best sequences of the campaign.

    This is the analogue of the consensus logo both source papers report. One panel per group,
    frequency encoded by height rather than by hue, so no colour has to carry six identities.
    """
    plt = mdtools.apply_style()
    best = top_sequences(iterations, top_n)
    sequences = [r["sequence"] for r in best]
    freq, length = mdtools.positional_group_frequency(sequences)
    if not freq:
        return False

    groups = list(mdtools.RESIDUE_GROUPS)
    positions = np.arange(1, length + 1)
    fig, axes = plt.subplots(len(groups), 1, figsize=(max(6, length * 0.34), 5.6), sharex=True)
    for i, (ax, group) in enumerate(zip(axes, groups)):
        # 4px rounded data-ends anchored to the baseline; 2px gap between adjacent bars.
        ax.bar(positions, freq[group], color=mdtools.SEQUENTIAL, width=0.62)
        ax.set_ylabel(group, rotation=0, ha="right", va="center",
                      fontsize=8, color=mdtools.INK)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["", ""] if i else ["0", "1"])
        # Only a baseline: a y spine per panel would be six vertical rules down the figure.
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[-1].set_xlabel("sequence position")
    axes[-1].set_xticks(positions[::2] if length > 20 else positions)
    fig.suptitle("Residue group frequency by position, top %d sequences" % len(sequences),
                 x=0.005, ha="left", fontsize=9.5, color=mdtools.INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=0.6)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def print_summary(settings, rows, iterations, last_improved, top_n, unparsed):
    if not rows:
        print("No iterations parsed. Is this a progress file written by the GA?")
        return

    first, last = rows[0], rows[-1]
    total_eval = sum(r["n_evaluated"] for r in rows)
    all_sequences = {r["sequence"] for recs in iterations.values() for r in recs}

    print("CAMPAIGN")
    for key in ("Iterations", "Population Size", "Number of parents", "Scoring function",
                "Selection method", "Constraint function", "Genes"):
        if key in settings:
            print("  %-22s %s" % (key, settings[key]))
    print("  %-22s %d" % ("iterations present", len(rows)))
    print("  %-22s %d (%d distinct sequences)" % ("evaluations", total_eval, len(all_sequences)))
    if unparsed:
        print("  %-22s %d rows did not match the expected format" % ("WARNING", unparsed))
    print()

    print("PROGRESS")
    header = ("iter", "best", "parents", "mean", "sd", "uniq", "hamming", "entropy", "charge")
    print("  " + "  ".join("%8s" % h for h in header))
    step = max(1, len(rows) // 20)
    for r in rows[::step] + ([rows[-1]] if len(rows) > 1 and (len(rows) - 1) % step else []):
        print("  " + "  ".join([
            "%8d" % r["iteration"], "%8.3f" % r["fitness_best"],
            "%8.3f" % r["fitness_parents_mean"], "%8.3f" % r["fitness_mean"],
            "%8.3f" % r["fitness_sd"], "%8d" % r["n_unique"],
            "%8.3f" % r["hamming_mean"], "%8.3f" % r["entropy_mean"],
            "%8.2f" % r["charge_mean"],
        ]))
    print()

    print("CONVERGENCE")
    print("  best fitness           %.3f -> %.3f" % (first["fitness_best"], last["fitness_best"]))
    print("  mean fitness           %.3f -> %.3f" % (first["fitness_mean"], last["fitness_mean"]))
    print("  mean pairwise distance %.3f -> %.3f" % (first["hamming_mean"], last["hamming_mean"]))
    print("  positional entropy     %.3f -> %.3f bits" % (first["entropy_mean"], last["entropy_mean"]))
    if last_improved is not None:
        stalled = last["iteration"] - last_improved
        print("  last improvement       iteration %d (%d ago)" % (last_improved, stalled))
    print()

    # The reading is stated rather than left to the reader of a table, because the failure mode
    # this report exists to catch looks like success on the fitness panel alone.
    collapsed = last["hamming_mean"] < 0.15 and last["entropy_mean"] < 0.5
    improving = last["fitness_best"] > first["fitness_best"]
    if collapsed and improving:
        print("  Reading: fitness improved but the population has collapsed. Diversity this low")
        print("  means later iterations are resampling one sequence family, so further iterations")
        print("  will cost simulations and return little. Check whether the entropy fell at a few")
        print("  positions, which is a motif, or across all of them, which is premature")
        print("  convergence.")
    elif collapsed:
        print("  Reading: the population has collapsed without improving fitness. Raise the")
        print("  mutation rate or the parent pool before spending more compute.")
    elif not improving:
        print("  Reading: no net improvement in best fitness. If the observable itself does not")
        print("  separate the calibration set, no optimiser setting fixes this.")
    print()

    print("TOP %d SEQUENCES" % top_n)
    print("  %-34s %9s %9s %7s %6s %5s" % ("sequence", "fitness", "raw", "orient", "charge", "runs"))
    for r in top_sequences(iterations, top_n):
        print("  %-34s %9.3f %9.3f %7.3f %6d %5d" % (
            r["sequence"][:34], r["fitness"], r["raw"], r["orientation"],
            mdtools.net_charge(r["sequence"]), r["n_run"]))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("progress", help="the GA progress file")
    parser.add_argument("--outdir", default="report", help="where figures and the CSV go")
    parser.add_argument("--top", type=int, default=15, help="how many top sequences to report")
    parser.add_argument("--parents", type=int, default=None,
                        help="parent pool size; read from the progress file when omitted")
    parser.add_argument("--no-figures", action="store_true", help="tables only")
    args = parser.parse_args()

    settings, iterations, unparsed = parse_progress(args.progress)
    if not iterations:
        sys.exit("No iteration blocks found in %s" % args.progress)

    n_parents = args.parents
    if n_parents is None:
        try:
            n_parents = int(settings.get("Number of parents", 36))
        except ValueError:
            n_parents = 36

    rows, last_improved = summarise(iterations, n_parents)
    os.makedirs(args.outdir, exist_ok=True)

    write_csv(rows, os.path.join(args.outdir, "progress.csv"))
    with open(os.path.join(args.outdir, "summary.json"), "w") as handle:
        json.dump({"settings": settings, "iterations": rows,
                   "top": top_sequences(iterations, args.top)}, handle, indent=2)

    print_summary(settings, rows, iterations, last_improved, args.top, unparsed)

    if not args.no_figures:
        figure_performance(rows, os.path.join(args.outdir, "performance.png"), last_improved)
        figure_composition(rows, os.path.join(args.outdir, "composition.png"))
        figure_positional(iterations, os.path.join(args.outdir, "positional.png"), args.top)
        print()
        print("Wrote progress.csv, summary.json, performance.png, composition.png, "
              "positional.png to %s" % args.outdir)


if __name__ == "__main__":
    main()
