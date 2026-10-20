#!/usr/bin/env python3
"""
Nature-style plotting theme for the AMP-evo-MD analysis figures.

This module is the single source of visual style for every analysis figure in the pipeline
(membrane QC, GA monitoring, per-candidate detail, cohort comparison, amino-acid analyses). It
is deliberately self-contained: it depends only on matplotlib and numpy, so the figures render
the same way on a workstation, on a cluster login node, or in a headless re-analysis, without
any external design package being installed.

Design decisions, fixed once here so they are not re-taken per figure:

  - Palette: the Nature Publishing Group (NPG) qualitative set. The two membrane roles have
    reserved slots that never change meaning: bacterial (the target) is red, host is dark blue.
  - Continuous ramps are derived from the same NPG hues (NPG_SEQ sequential, NPG_DIV diverging)
    rather than viridis, so a heatmap and a line plot in the same figure read as one palette.
  - Typography: Helvetica/Arial at ~7 pt, the Nature single-column body size.
  - Spines: left and bottom only, thin; ticks point outward; no grid; white background.
  - No bar charts anywhere. Magnitudes are drawn as Cleveland dot plots or lollipops, and
    two-group comparisons as dumbbells, which carry the same information with far less ink and
    read the difference (the selectivity gap) directly.

The helpers at the bottom (dumbbell, lollipop, cleveland) are the non-bar primitives every
report uses so the no-bar rule is enforced by construction rather than by review.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless by default; a workstation session can switch backend before import
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap

# --- colours -----------------------------------------------------------------------------------
# NPG qualitative palette (ggsci npg_nrc). Order is fixed: series are assigned in this order and
# never cycled past slot 8, so the same observable keeps the same colour across figures.
NPG = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488",
       "#F39B7F", "#8491B4", "#91D1C2", "#7E6148"]

# Reserved role colours. These two are semantic and must not be reused for anything else.
TARGET = NPG[0]   # bacterial membrane, the target to damage (red)
HOST = NPG[3]     # host membrane, to be spared (dark blue)

# Neutral inks used for axes, reference lines and faint connectors.
INK = "#222222"
GREY = "#9aa0a6"
FAINT = "#c9ccd1"

# Continuous ramps built from the NPG hues so every colour in a figure stays in the palette.
# NPG_SEQ: low-to-high magnitude (pale to deep NPG blue). NPG_DIV: signed quantities centred on
# white, host-blue negative to bacterial-red positive, matching the role colours.
NPG_SEQ = LinearSegmentedColormap.from_list(
    "npg_seq", ["#f5f8fb", "#c3d0e0", "#7f97b8", "#3C5488", "#243654"])
NPG_DIV = LinearSegmentedColormap.from_list(
    "npg_div", ["#3C5488", "#8fa6c4", "#eef0f2", "#f0a58f", "#E64B35"])

# Residue-group colours, aligned to the NPG set so sequence maps match the rest of the palette.
GROUP_COLORS = {
    "positive": NPG[0], "aromatic": NPG[3], "hydrophobic": NPG[2],
    "polar": NPG[1], "negative": NPG[5], "other": NPG[7],
}


# --- rc profile --------------------------------------------------------------------------------
def rc(base=7.0):
    """
    Return the matplotlib rcParams dict for the Nature theme.

    Kept as a function returning a dict (rather than applied on import) so a caller can scale the
    base font size for a poster or a thumbnail without editing the module.
    """
    return {
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "savefig.dpi": 300, "figure.dpi": 150, "savefig.bbox": "tight",
        "figure.constrained_layout.use": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": base, "axes.labelsize": base, "axes.titlesize": base + 1,
        "xtick.labelsize": base - 1, "ytick.labelsize": base - 1, "legend.fontsize": base - 1,
        "axes.linewidth": 0.75, "axes.edgecolor": INK, "axes.labelcolor": INK,
        "axes.titlecolor": "#000000", "axes.titlelocation": "left",
        "axes.titleweight": "bold", "axes.titlepad": 6,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.75, "ytick.major.width": 0.75,
        "xtick.color": INK, "ytick.color": INK,
        "lines.linewidth": 1.0, "lines.markersize": 4, "lines.solid_capstyle": "round",
        "legend.frameon": False, "text.color": INK,
        "axes.prop_cycle": cycler(color=NPG),
    }


def use(base=7.0):
    """Apply the Nature theme globally and return the pyplot module for convenience."""
    plt.rcParams.update(rc(base))
    return plt


# Figure widths, in inches, matching Nature column specifications.
WIDTH_SINGLE = 3.50   # ~89 mm single column
WIDTH_ONEHALF = 4.80  # ~122 mm
WIDTH_DOUBLE = 7.20   # ~183 mm double column


# --- figure helpers ----------------------------------------------------------------------------
def panel(ax, letter):
    """Bold lowercase panel letter at the top-left corner, Nature convention."""
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="right", color="black")


def finish(ax, title=None, xlabel=None, ylabel=None, letter=None):
    """Apply titles/labels and the panel letter in one call, keeping report code uncluttered."""
    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if letter is not None:
        panel(ax, letter)


def colorbar(fig, ax, mappable, label=None):
    """A compact colourbar sized for a single-column panel."""
    cb = fig.colorbar(mappable, ax=ax, fraction=0.045, pad=0.02)
    if label:
        cb.set_label(label, fontsize=6)
    cb.ax.tick_params(labelsize=5)
    return cb


def save(fig, path, dpi=300):
    """Save and close a figure; returns the path so callers can collect what they wrote."""
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


# --- non-bar primitives ------------------------------------------------------------------------
def dumbbell(ax, labels, left, right, left_color=HOST, right_color=TARGET,
             left_label="host", right_label="bacterial", sort_by_gap=True):
    """
    Horizontal dumbbell plot: one row per label, two dots joined by a faint bar.

    Replaces a grouped bar chart for two-group comparisons. The connecting bar makes the gap
    (for the membrane plots, the selectivity gap) the most salient feature, which is the whole
    point of comparing target against host.
    """
    left = np.asarray(left, float)
    right = np.asarray(right, float)
    order = np.argsort(right - left) if sort_by_gap else np.arange(len(labels))
    y = np.arange(len(labels))
    for yi, i in zip(y, order):
        ax.plot([left[i], right[i]], [yi, yi], color=FAINT, lw=1.6, zorder=1)
    ax.plot(right[order], y, "o", color=right_color, ms=6, label=right_label, zorder=2)
    ax.plot(left[order], y, "o", color=left_color, ms=6, label=left_label, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[i] for i in order])
    return order


def lollipop(ax, labels, values, baseline=0.0, pos_color=TARGET, neg_color=HOST,
             sort=True, orient="h"):
    """
    Lollipop / Cleveland dot plot: a faint stem from a baseline to a dot.

    Replaces a bar chart for ranked magnitudes or signed effects. Colour splits on the sign
    relative to the baseline, so favoured and rejected residues separate at a glance.
    """
    values = np.asarray(values, float)
    order = np.argsort(values) if sort else np.arange(len(labels))
    pos = np.arange(len(labels))
    for p, i in zip(pos, order):
        c = pos_color if values[i] >= baseline else neg_color
        if orient == "h":
            ax.plot([baseline, values[i]], [p, p], color=FAINT, lw=0.9, zorder=1)
            ax.plot(values[i], p, "o", color=c, ms=5, zorder=2)
        else:
            ax.plot([p, p], [baseline, values[i]], color=FAINT, lw=0.9, zorder=1)
            ax.plot(p, values[i], "o", color=c, ms=5, zorder=2)
    if orient == "h":
        ax.set_yticks(pos)
        ax.set_yticklabels([labels[i] for i in order])
        ax.axvline(baseline, color=INK, lw=0.75)
    else:
        ax.set_xticks(pos)
        ax.set_xticklabels([labels[i] for i in order])
        ax.axhline(baseline, color=INK, lw=0.75)
    return order


# --- sequence descriptors ----------------------------------------------------------------------
# Eisenberg consensus hydrophobicity scale, used only for the hydrophobic-moment descriptor.
_EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74,
    "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19,
    "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}


def hydrophobic_moment(sequence, angle_deg=100.0):
    """
    Eisenberg hydrophobic moment for an ideal alpha-helix (100 deg per residue).

    A sequence-only descriptor of amphipathicity, independent of the simulation. Returned per
    residue (divided by length) so peptides of different length are comparable.
    """
    if not sequence:
        return float("nan")
    ang = np.deg2rad(angle_deg)
    sx = sum(_EISENBERG.get(a, 0.0) * np.cos(i * ang) for i, a in enumerate(sequence))
    sy = sum(_EISENBERG.get(a, 0.0) * np.sin(i * ang) for i, a in enumerate(sequence))
    return float(np.hypot(sx, sy) / len(sequence))


def mean_hydrophobicity(sequence):
    """Mean Eisenberg hydrophobicity over the sequence."""
    if not sequence:
        return float("nan")
    return float(np.mean([_EISENBERG.get(a, 0.0) for a in sequence]))
