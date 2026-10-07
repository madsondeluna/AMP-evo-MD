"""
Shared helpers for the analysis scripts: file parsing, block averaging, and plot style.

Nothing here runs a simulation or depends on GROMACS being installed. The heavy extraction is
done by gmx itself (gmx energy, gmx density), which writes .xvg; these functions read what it
wrote. That keeps the analysis runnable on a workstation against files copied off the cluster.
"""

import collections
import math
import os
import re

import numpy as np


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------
# Fixed categorical order, never cycled. Four slots, validated for colour-vision deficiency
# separation, chroma, lightness band and contrast against a light surface. A fifth series does
# not get a generated colour: facet into small multiples instead, which is what the composition
# and positional figures do.
CATEGORICAL = ["#0072B2", "#009E73", "#D55E00", "#8C5AA8"]

# Single hue, light to dark, for magnitude. Never a rainbow.
SEQUENTIAL = "#0072B2"

INK = "#1a1a1a"
INK_MUTED = "#6b6b6b"
GRID = "#d9d9d6"
SURFACE = "#fcfcfb"


def apply_style():
    """
    Minimal chart style: no gridlines, two spines, ink-token text, thin marks.

    Gridlines are off entirely. They exist to let a reader recover a value from a mark, which a
    tick label and a direct annotation do better and with less ink. Everything that is not data
    is either removed or pushed back: the top and right spines go, the remaining two sit in the
    muted grid colour rather than black, and tick marks are drawn outward and short.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK_MUTED,
        "axes.titlesize": 9.5,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.titlecolor": INK,
        "axes.labelsize": 8.5,
        "axes.labelpad": 6,
        "axes.linewidth": 0.7,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.size": 0,
        "ytick.minor.size": 0,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.handlelength": 1.4,
        "legend.borderpad": 0,
        "legend.labelspacing": 0.45,
        "lines.linewidth": 1.5,
        "lines.solid_capstyle": "round",
        "lines.markersize": 4,
        "font.size": 8.5,
        "figure.dpi": 150,
    })
    return plt


def declutter(ax, ymax_ticks=4, xmax_ticks=6):
    """
    Thin out tick labels and lift the spines off the data.

    Matplotlib's default tick density is set for reading values off the axis, which is not what
    these charts are for: the reader wants the shape, and the exact numbers are in the CSV. Fewer
    ticks and detached spines leave the data as the only thing with weight on the page.
    """
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(nbins=ymax_ticks, prune=None))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=xmax_ticks, integer=True))
    for side in ("left", "bottom"):
        ax.spines[side].set_position(("outward", 6))
    return ax


def label_last(ax, x, y, text, colour, dx=6):
    """
    Direct-label a series at its right end, instead of relying on the legend alone.

    Identity then survives being read in greyscale or by someone who cannot separate two of the
    hues, which is the whole reason not to let colour carry meaning by itself.
    """
    if len(x) == 0:
        return
    ax.annotate(text, xy=(x[-1], y[-1]), xytext=(dx, 0), textcoords="offset points",
                va="center", ha="left", fontsize=7.5, color=colour)


# ---------------------------------------------------------------------------
# File parsing
# ---------------------------------------------------------------------------
def read_xvg(path):
    """
    Read a GROMACS .xvg into (x, {label: y}).

    Series labels come from the '@ s<N> legend "..."' lines when present, otherwise the columns
    are named col1, col2 and so on. Comment and command lines are skipped.
    """
    legends = {}
    xlabel = "x"
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if line.startswith("@"):
                m = re.match(r'@\s+s(\d+)\s+legend\s+"(.*)"', line)
                if m:
                    legends[int(m.group(1))] = m.group(2)
                m = re.match(r'@\s+xaxis\s+label\s+"(.*)"', line)
                if m:
                    xlabel = m.group(1)
                continue
            parts = line.split()
            if not parts:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError("no numeric rows in %s" % path)
    # Use the most common row width and drop rows that do not match, rather than truncating every
    # column to the shortest row. Reading an .xvg while gmx is still flushing leaves a final
    # partial line, and taking its width would silently reduce the file to a single column: the
    # caller would then see no series at all and report a result with checks quietly missing.
    widths = collections.Counter(len(r) for r in rows)
    width = widths.most_common(1)[0][0]
    data = np.array([r[:width] for r in rows if len(r) >= width])
    x = data[:, 0]
    series = collections.OrderedDict()
    for i in range(1, width):
        series[legends.get(i - 1, "col%d" % i)] = data[:, i]
    return x, series, xlabel


def read_gro(path):
    """Return (box_x, box_y, box_z, atoms) with atoms as (resname, atomname, x, y, z)."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    n_atoms = int(lines[1].strip())
    atoms = []
    for line in lines[2:2 + n_atoms]:
        try:
            atoms.append((line[5:10].strip(), line[10:15].strip(),
                          float(line[20:28]), float(line[28:36]), float(line[36:44])))
        except ValueError:
            continue
    box = lines[2 + n_atoms].split()
    return float(box[0]), float(box[1]), float(box[2]), atoms


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def block_average(values, n_blocks=5):
    """
    Mean and block-averaged standard error.

    A plain standard error over an MD time series understates the uncertainty badly, because
    consecutive frames are correlated. Splitting into blocks and taking the standard error of the
    block means is the cheap correction: it is still an underestimate if the blocks are shorter
    than the correlation time, but it is far closer than treating every frame as independent.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return float("nan"), float("nan")
    if values.size < n_blocks * 2:
        return float(values.mean()), float(values.std(ddof=1) / math.sqrt(values.size)) if values.size > 1 else 0.0
    blocks = np.array_split(values, n_blocks)
    means = np.array([b.mean() for b in blocks])
    return float(values.mean()), float(means.std(ddof=1) / math.sqrt(n_blocks))


def drift(x, y):
    """Least-squares slope of y against x, as a convergence check. Units are y per unit x."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2:
        return float("nan")
    return float(np.polyfit(x, y, 1)[0])


def tail(x, y, fraction):
    """Keep the last `fraction` of a time series, for averaging after equilibration."""
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    start = x[0] + (1 - fraction) * (x[-1] - x[0])
    mask = x >= start
    return x[mask], y[mask]


# ---------------------------------------------------------------------------
# Sequence chemistry
# ---------------------------------------------------------------------------
# Residue groups follow the scheme used by Biriukov and Vacha, so the composition figures here
# can be read against theirs. The grouping is coarse on purpose: at coarse-grained resolution
# broad physicochemical classes are resolved far more robustly than individual residue
# preferences, and reporting per-residue frequencies invites over-reading of noise.
RESIDUE_GROUPS = collections.OrderedDict([
    ("aromatic",  set("FWY")),
    ("hydrophobic", set("ALMIV")),
    ("polar", set("QSTN")),
    ("positive", set("KR")),
    ("negative", set("DE")),
    ("other", set("GPCH")),
])

GROUP_OF = {}
for _group, _members in RESIDUE_GROUPS.items():
    for _aa in _members:
        GROUP_OF[_aa] = _group


def net_charge(sequence):
    """Net charge at pH 7, histidine neutral, matching the GA's own convention."""
    return (sequence.count("K") + sequence.count("R")) - (sequence.count("D") + sequence.count("E"))


def group_counts(sequence):
    counts = collections.OrderedDict((g, 0) for g in RESIDUE_GROUPS)
    for aa in sequence:
        g = GROUP_OF.get(aa)
        if g:
            counts[g] += 1
    return counts


def hamming(a, b):
    """Mismatches between two equal-length sequences; None if the lengths differ."""
    if len(a) != len(b):
        return None
    return sum(1 for x, y in zip(a, b) if x != y)


def mean_pairwise_hamming(sequences, max_pairs=20000, rng=None):
    """
    Mean pairwise Hamming distance, normalised by length.

    Exhaustive for small populations; for larger ones a random sample of pairs, because the pair
    count grows quadratically and the estimate converges long before the full set is needed.
    Sequences of differing length are skipped rather than padded: padding would invent identity
    at the ends and bias the measure downward.
    """
    seqs = [s for s in sequences if s]
    n = len(seqs)
    if n < 2:
        return float("nan")
    pairs = n * (n - 1) // 2
    total = 0.0
    used = 0
    if pairs <= max_pairs:
        for i in range(n):
            for j in range(i + 1, n):
                d = hamming(seqs[i], seqs[j])
                if d is not None:
                    total += d / len(seqs[i])
                    used += 1
    else:
        import random
        rng = rng or random.Random(0)
        for _ in range(max_pairs):
            i, j = rng.randrange(n), rng.randrange(n)
            if i == j:
                continue
            d = hamming(seqs[i], seqs[j])
            if d is not None:
                total += d / len(seqs[i])
                used += 1
    return total / used if used else float("nan")


def positional_entropy(sequences):
    """
    Shannon entropy in bits at each sequence position.

    Zero means every sequence in the population carries the same residue there; the maximum is
    log2 of the alphabet size, reached when the position is still uniformly random. Watching this
    fall is the most direct view of where the search has converged and where it is still
    exploring, and it separates a population that has genuinely found a motif from one that has
    merely collapsed onto a few near-identical survivors.
    """
    seqs = [s for s in sequences if s]
    if not seqs:
        return np.array([])
    length = min(len(s) for s in seqs)
    out = np.zeros(length)
    for pos in range(length):
        counts = collections.Counter(s[pos] for s in seqs)
        total = sum(counts.values())
        h = 0.0
        for c in counts.values():
            p = c / total
            h -= p * math.log2(p)
        out[pos] = h
    return out


def positional_group_frequency(sequences):
    """Frequency of each residue group at each position, as {group: array(length)}."""
    seqs = [s for s in sequences if s]
    if not seqs:
        return {}, 0
    length = min(len(s) for s in seqs)
    freq = {g: np.zeros(length) for g in RESIDUE_GROUPS}
    for pos in range(length):
        for s in seqs:
            g = GROUP_OF.get(s[pos])
            if g:
                freq[g][pos] += 1
    for g in freq:
        freq[g] /= len(seqs)
    return freq, length
