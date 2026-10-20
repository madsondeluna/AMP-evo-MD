#!/usr/bin/env python3
"""
Render the AMP-evo-MD analysis figures from the saved dataset.

Every figure here reads only the aggregated dataset written by amp_dataset.py
(cohort.csv, sequences.csv, metadata.json) and, for the per-candidate trajectory plots, the
per-candidate series under analysis/<id>/. Nothing in this module re-reads a raw observables.json,
re-parses the progress file, or touches a trajectory: redrawing a figure never recomputes a
measurement. Run amp_dataset.py first; run amp_traj.py to produce the per-candidate series for
the trajectory plots.

Figure groups:
  cohort     C1-C8   compare the whole cohort (100+ AMPs) on the saved observables
  aa         A1-A7   amino-acid selectivity analyses (bacterial target vs host)
  ga         E4-E8   genetic-algorithm monitoring (needs the progress file)
  membrane   E1-E3   membrane QC and reference values
  candidate  P1-P13  per-candidate detail (needs analysis/<id>/ series from amp_traj.py)

All figures use the vendored Nature theme in ampstyle.py: NPG palette, NPG-derived continuous
ramps, no grid, no bar charts (magnitudes are dots/lollipops, comparisons are dumbbells).

USAGE
-----
    python tools/amp_report.py --dataset <campaign>/analysis --group cohort aa --outdir figures/
    python tools/amp_report.py --dataset <campaign>/analysis --group all --sheet
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ampstyle as A
import amp_dataset
import mdtools

plt = A.use()

# Registry of figures: (group, code, title, function). Each function takes (ctx) and draws one
# figure, returning the saved path. Decorator @figure records it so the CLI can select by group.
_FIGURES = []


def figure(group, code, title):
    def deco(fn):
        _FIGURES.append((group, code, title, fn))
        return fn
    return deco


class Context:
    """Holds the loaded dataset and the output directory, passed to every figure function."""
    def __init__(self, dataset_dir, outdir):
        self.dataset_dir = dataset_dir
        self.outdir = outdir
        self.cohort = amp_dataset.load_cohort(os.path.join(dataset_dir, "cohort.csv"))
        seq_path = os.path.join(dataset_dir, "sequences.csv")
        self.seqs = amp_dataset.load_cohort(seq_path) if os.path.isfile(seq_path) else {}
        import json
        meta_path = os.path.join(dataset_dir, "metadata.json")
        self.meta = json.load(open(meta_path)) if os.path.isfile(meta_path) else {}
        self.alphabet = self.meta.get("active_alphabet", list("AEFKLMQRSWY"))
        self.candidate = None  # set by the CLI for the per-candidate (P) group

    def _load_optional(self, name):
        path = os.path.join(self.dataset_dir, name)
        return amp_dataset.load_cohort(path) if os.path.isfile(path) else {}

    def timeline(self):
        if not hasattr(self, "_timeline"):
            self._timeline = self._load_optional("ga_timeline.csv")
        return self._timeline

    def positional(self):
        if not hasattr(self, "_positional"):
            self._positional = self._load_optional("positional.csv")
        return self._positional

    def candidate_dir(self):
        """Per-candidate series directory written by amp_traj.py: <dataset>/<id>/."""
        return os.path.join(self.dataset_dir, str(self.candidate)) if self.candidate else None

    def campaign_candidate_dir(self):
        """The original GROMACS candidate directory, for series kept by the run (distance xvg)."""
        root = self.meta.get("campaign_dir")
        return os.path.join(root, str(self.candidate)) if root and self.candidate else None

    def series(self, name):
        """A per-candidate series CSV from amp_traj, as {column: array}, or {} if absent."""
        d = self.candidate_dir()
        path = os.path.join(d, name) if d else None
        return amp_dataset.load_cohort(path) if path and os.path.isfile(path) else {}

    def col(self, name):
        """A cohort column as a float array, or an array of nan if absent."""
        v = self.cohort.get(name)
        if v is None:
            n = len(next(iter(self.cohort.values()))) if self.cohort else 0
            return np.full(n, np.nan)
        return v

    def counts_matrix(self):
        """Per-residue count matrix (n_candidates x n_alphabet) from sequences.csv."""
        return np.column_stack([self.seqs.get("n_%s" % a,
                               np.zeros(len(self.seqs.get("id", [])))) for a in self.alphabet])

    def new(self, w=A.WIDTH_SINGLE, h=3.0):
        return plt.subplots(figsize=(w, h))

    def out(self, name):
        return os.path.join(self.outdir, name)


def _finite(*arrays):
    """Mask keeping rows where every passed array is finite."""
    m = np.ones(len(arrays[0]), bool)
    for a in arrays:
        m &= np.isfinite(a)
    return m


def _corr(x, y):
    m = _finite(x, y)
    if m.sum() < 3 or np.std(x[m]) < 1e-12 or np.std(y[m]) < 1e-12:
        return 0.0
    return float(np.corrcoef(x[m], y[m])[0, 1])


def _zscore(M):
    return (M - np.nanmean(M, 0)) / (np.nanstd(M, 0) + 1e-9)


# ===================================================================================== COHORT ==
def _cohort_matrix(ctx):
    """The observable matrix used by the multivariate cohort figures, with column labels."""
    cols = [("damage", "damage"), ("host_damage", "host dmg"), ("selectivity", "selectivity"),
            ("elip_target", "|E_lip|"), ("dmem_target", "d_mem"), ("orient_target", "orient"),
            ("charge", "charge"), ("hydrophobic_moment", "H-mom")]
    labels, series = [], []
    for key, lab in cols:
        v = ctx.col(key).astype(float)
        if key == "elip_target":
            v = -v / 100.0  # make magnitude positive and on a comparable scale
        if np.isfinite(v).sum() >= 3:
            labels.append(lab)
            series.append(v)
    return np.column_stack(series), labels


@figure("cohort", "C1", "Candidate by observable")
def c1_heatmap(ctx):
    M, labels = _cohort_matrix(ctx)
    fit = ctx.col("fitness")
    order = np.argsort(-np.where(np.isfinite(fit), fit, -np.inf))
    Z = _zscore(M)[order]
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.4)
    im = ax.imshow(Z, aspect="auto", cmap=A.NPG_DIV, vmin=-2.5, vmax=2.5,
                   extent=[0, len(labels), len(order), 0])
    ax.set_xticks(np.arange(len(labels)) + 0.5)
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=5.5)
    ax.set_yticks([5, len(order) - 3]); ax.set_yticklabels(["best", "worst"])
    A.colorbar(fig, ax, im, "z")
    A.finish(ax, "Candidate by observable", None, None, "a")
    return A.save(fig, ctx.out("C1_candidate_observable.png"))


@figure("cohort", "C2", "Observable distributions")
def c2_dist(ctx):
    specs = [("damage", A.TARGET, "damage"), ("host_damage", A.HOST, "host dmg"),
             ("selectivity", A.NPG[2], "selectivity"), ("dmem_target", A.NPG[4], "d_mem")]
    data, cols, labs = [], [], []
    for key, c, lab in specs:
        v = ctx.col(key); v = v[np.isfinite(v)]
        if v.size:
            data.append(v); cols.append(c); labs.append(lab)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    bp = ax.boxplot(data, widths=0.5, showfliers=False, patch_artist=True)
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.2); patch.set_edgecolor(c)
    for med in bp["medians"]:
        med.set_color(A.INK)
    rng = np.random.default_rng(0)
    for i, (d, c) in enumerate(zip(data, cols), 1):
        ax.plot(i + 0.1 * rng.standard_normal(d.size), d, "o", ms=2, color=c, alpha=0.5)
    ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs)
    A.finish(ax, "Observable distributions", None, "value", "b")
    return A.save(fig, ctx.out("C2_distributions.png"))


@figure("cohort", "C3", "Damage vs selectivity")
def c3_pareto(ctx):
    dmg, sel, fit = ctx.col("damage"), ctx.col("selectivity"), ctx.col("fitness")
    m = _finite(dmg, sel)
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.0)
    mx, my = np.nanmedian(dmg[m]), np.nanmedian(sel[m])
    ax.add_patch(plt.Rectangle((mx, my), np.nanmax(dmg[m]) + 0.4 - mx, np.nanmax(sel[m]) + 0.4 - my,
                 facecolor=A.TARGET, alpha=0.06, zorder=0))
    ax.axvline(mx, color=A.GREY, lw=0.6, ls="--"); ax.axhline(my, color=A.GREY, lw=0.6, ls="--")
    c = fit[m] if np.isfinite(fit[m]).any() else None
    sc = ax.scatter(dmg[m], sel[m], c=c, cmap=A.NPG_SEQ, s=16, edgecolor="white", linewidth=0.3)
    if c is not None:
        A.colorbar(fig, ax, sc, "fitness")
    A.finish(ax, "Damage vs selectivity", r"damage (kJ mol$^{-1}$ nm$^{-2}$)", "selectivity", "c")
    return A.save(fig, ctx.out("C3_pareto.png"))


@figure("cohort", "C4", "Shortlist")
def c4_short(ctx):
    sel = ctx.col("selectivity"); ids = ctx.cohort.get("id")
    m = _finite(sel); idx = np.where(m)[0]
    order = idx[np.argsort(sel[idx])[-15:]]
    err = ctx.col("dg_err_target")
    y = np.arange(len(order))
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.2)
    for i, o in enumerate(order):
        e = err[o] if np.isfinite(err[o]) else 0.0
        ax.plot([sel[o] - e, sel[o] + e], [y[i], y[i]], color=A.FAINT, lw=0.9)
        ax.plot(sel[o], y[i], "o", color=A.TARGET, ms=4)
    labels = [str(ids[o])[:10] if ids is not None else str(o) for o in order]
    ax.set_yticks(y[::2]); ax.set_yticklabels(labels[::2], fontsize=5.5)
    A.finish(ax, "Shortlist (top 15 by selectivity)", "selectivity", None, "d")
    return A.save(fig, ctx.out("C4_shortlist.png"))


@figure("cohort", "C5", "PCA biplot")
def c5_pca(ctx):
    M, labels = _cohort_matrix(ctx)
    m = np.all(np.isfinite(M), 1)
    Z = _zscore(M[m]); fit = ctx.col("fitness")[m]
    U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    pc = U[:, :2] * S[:2]; vr = (S ** 2 / np.sum(S ** 2)) * 100
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.2)
    c = fit if np.isfinite(fit).any() else None
    sc = ax.scatter(pc[:, 0], pc[:, 1], c=c, cmap=A.NPG_SEQ, s=16, edgecolor="white", linewidth=0.3)
    scl = 0.9 * np.abs(pc).max()
    for k, n in enumerate(labels):
        ax.annotate("", (Vt[0, k] * scl, Vt[1, k] * scl), (0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=A.GREY, lw=0.9))
        ax.text(Vt[0, k] * scl * 1.07, Vt[1, k] * scl * 1.07, n, fontsize=5, ha="center")
    if c is not None:
        A.colorbar(fig, ax, sc, "fitness")
    A.finish(ax, "PCA biplot", "PC1 (%.0f%%)" % vr[0], "PC2 (%.0f%%)" % vr[1], "e")
    return A.save(fig, ctx.out("C5_pca.png"))


@figure("cohort", "C6", "Observable correlations")
def c6_corr(ctx):
    M, labels = _cohort_matrix(ctx)
    m = np.all(np.isfinite(M), 1)
    R = np.corrcoef(M[m], rowvar=False)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.2)
    im = ax.imshow(R, cmap=A.NPG_DIV, vmin=-1, vmax=1)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=5.5)
    A.colorbar(fig, ax, im, "r")
    A.finish(ax, "Observable correlations", None, None, "f")
    return A.save(fig, ctx.out("C6_correlations.png"))


@figure("cohort", "C7", "Charge vs selectivity")
def c7_charge(ctx):
    q, sel, hm = ctx.col("charge"), ctx.col("selectivity"), ctx.col("hydrophobic_moment")
    m = _finite(q, sel)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    c = hm[m] if np.isfinite(hm[m]).any() else None
    sc = ax.scatter(q[m], sel[m], c=c, cmap=A.NPG_SEQ, s=16, edgecolor="white", linewidth=0.3)
    if m.sum() >= 2:
        b, a = np.polyfit(q[m], sel[m], 1); xs = np.array([q[m].min(), q[m].max()])
        ax.plot(xs, a + b * xs, color=A.GREY, lw=1.0, ls="--")
    if c is not None:
        A.colorbar(fig, ax, sc, "H-mom")
    A.finish(ax, "Charge vs selectivity", "net charge (+)", "selectivity", "g")
    return A.save(fig, ctx.out("C7_charge_selectivity.png"))


@figure("cohort", "C8", "Top sequences by residue group")
def c8_seq(ctx):
    from matplotlib.colors import ListedColormap, BoundaryNorm
    ids = ctx.cohort.get("id"); seqs = ctx.seqs.get("sequence")
    sel = ctx.col("selectivity")
    if seqs is None:
        seqs = ids
    groups = list(mdtools.RESIDUE_GROUPS.keys())
    order = np.argsort(-np.where(np.isfinite(sel), sel, -np.inf))[:20]
    L = max(len(str(seqs[o])) for o in order)
    grid = np.full((len(order), L), -1)
    for r, o in enumerate(order):
        for c, aa in enumerate(str(seqs[o])):
            g = mdtools.GROUP_OF.get(aa)
            if g in groups:
                grid[r, c] = groups.index(g)
    cmap = ListedColormap([A.GROUP_COLORS.get(g, A.GREY) for g in groups])
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.2)
    ax.imshow(np.ma.masked_less(grid, 0), aspect="auto", cmap=cmap,
              norm=BoundaryNorm(range(len(groups) + 1), len(groups)), extent=[1, L + 1, len(order), 0])
    ax.set_xticks([1, 6, 12, 18, 24]); ax.set_yticks([0.5, len(order) - 0.5])
    ax.set_yticklabels(["#1", "#%d" % len(order)])
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=A.GROUP_COLORS.get(g, A.GREY), label=g) for g in groups],
              loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=5.5)
    A.finish(ax, "Top 20 sequences by residue group", "position", None, "h")
    return A.save(fig, ctx.out("C8_top_sequences.png"))


# ========================================================================================= AA ==
def _aa_corr(ctx, target_col):
    """corr(residue count, target column) for every residue in the active alphabet."""
    C = ctx.counts_matrix(); y = ctx.seqs.get(target_col)
    return np.array([_corr(C[:, k], y) for k in range(len(ctx.alphabet))])


@figure("aa", "A1", "Residue effect: bacterial vs host")
def a1_bachost(ctx):
    rb = _aa_corr(ctx, "damage"); rh = _aa_corr(ctx, "host_damage")
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.2)
    ax.add_patch(plt.Rectangle((0, min(rh.min(), -0.05)), 1.1, 0 - min(rh.min(), -0.05),
                 facecolor=A.TARGET, alpha=0.08, zorder=0))
    ax.axvline(0, color=A.GREY, lw=0.6, ls="--"); ax.axhline(0, color=A.GREY, lw=0.6, ls="--")
    ax.scatter(rb, rh, s=20, color=A.HOST, edgecolor="white", linewidth=0.4)
    for k, a in enumerate(ctx.alphabet):
        ax.annotate(a, (rb[k], rh[k]), textcoords="offset points", xytext=(3, 2), fontsize=6)
    A.finish(ax, "Residue effect: bacterial vs host",
             "corr(count, bacterial damage)", "corr(count, host damage)", "a")
    ax.annotate("selective", (rb.max(), rh.min()), xytext=(-2, 2), textcoords="offset points",
                ha="right", va="bottom", fontsize=6, color=A.TARGET)
    return A.save(fig, ctx.out("A1_residue_bac_host.png"))


@figure("aa", "A2", "Residues favoured vs rejected")
def a2_enrich(ctx):
    r_sel = _aa_corr(ctx, "selectivity")
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.2)
    A.lollipop(ax, ctx.alphabet, r_sel, baseline=0.0)
    A.finish(ax, "Residues favoured vs rejected", "corr(count, selectivity)", "residue", "b")
    return A.save(fig, ctx.out("A2_favoured_rejected.png"))


@figure("aa", "A3", "Amino-acid frequency")
def a3_freq(ctx):
    C = ctx.counts_matrix()
    freq = C.sum(0) / C.sum() * 100 if C.sum() else np.zeros(len(ctx.alphabet))
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.2)
    A.lollipop(ax, ctx.alphabet, freq, baseline=0.0, pos_color=A.HOST, neg_color=A.HOST)
    ax.axvline(100.0 / len(ctx.alphabet), color=A.GREY, lw=0.8, ls="--")
    A.finish(ax, "Amino-acid frequency", "frequency (%)", "residue", "c")
    return A.save(fig, ctx.out("A3_frequency.png"))


@figure("aa", "A4", "Amino-acid populations over generations")
def a4_overgen(ctx):
    C = ctx.counts_matrix(); gen = ctx.seqs.get("generation")
    if gen is None or not np.isfinite(gen).any():
        return None
    gen = gen.astype(float); G = int(np.nanmax(gen)) + 1
    H = np.zeros((len(ctx.alphabet), G))
    for g in range(G):
        mask = gen == g
        if mask.any():
            H[:, g] = (C[mask].sum(0) / C[mask].sum() * 100) if C[mask].sum() else 0
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.0)
    im = ax.imshow(H, aspect="auto", cmap=A.NPG_SEQ, extent=[0, G, len(ctx.alphabet), 0])
    ax.set_yticks(np.arange(len(ctx.alphabet)) + 0.5); ax.set_yticklabels(ctx.alphabet)
    step = max(1, G // 6)
    ax.set_xticks(np.arange(0, G, step) + 0.5); ax.set_xticklabels(range(0, G, step))
    A.colorbar(fig, ax, im, "freq (%)")
    A.finish(ax, "Amino-acid populations over generations", "generation", "residue", "d")
    return A.save(fig, ctx.out("A4_over_generations.png"))


@figure("aa", "A5", "Selectivity by residue group")
def a5_group(ctx):
    C = ctx.counts_matrix()
    groups = ["positive", "aromatic", "hydrophobic", "polar", "negative"]
    gcnt = {g: np.zeros(C.shape[0]) for g in groups}
    for k, a in enumerate(ctx.alphabet):
        g = mdtools.GROUP_OF.get(a)
        if g in gcnt:
            gcnt[g] += C[:, k]
    dmg, hst = ctx.seqs.get("damage"), ctx.seqs.get("host_damage")
    rb = np.array([_corr(gcnt[g], dmg) for g in groups])
    rh = np.array([_corr(gcnt[g], hst) for g in groups])
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.0)
    ax.axvline(0, color=A.GREY, lw=0.6, ls="--")
    A.dumbbell(ax, groups, rh, rb)
    ax.legend(loc="lower right")
    A.finish(ax, "Selectivity by residue group", "corr(group count, damage)", "residue group", "e")
    return A.save(fig, ctx.out("A5_by_group.png"))


@figure("aa", "A6", "Positional discrimination")
def a6_position(ctx):
    seqs = ctx.seqs.get("sequence"); sel = ctx.seqs.get("selectivity")
    if seqs is None:
        return None
    L = min(len(str(s)) for s in seqs)
    score = np.zeros(L)
    for p in range(L):
        has = np.array([1.0 if str(seqs[i])[p] in ("K", "R") else 0.0 for i in range(len(seqs))])
        score[p] = _corr(has, sel)
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 2.6)
    ax.axvspan(0.5, 6.5, color="#eef1f5", zorder=0); ax.axvspan(L - 5.5, L + 0.5, color="#eef1f5", zorder=0)
    ax.axhline(0, color=A.INK, lw=0.75)
    for xi, v in zip(range(1, L + 1), score):
        ax.plot([xi, xi], [0, v], color=A.FAINT, lw=1.0, zorder=1)
        ax.plot(xi, v, "o", color=A.TARGET if v > 0 else A.HOST, ms=4, zorder=2)
    ax.set_xticks([1, 6, 12, 18, 24]); ax.set_xlim(0.5, L + 0.5)
    A.finish(ax, "Positional discrimination (cationic residues)",
             "sequence position", "corr(cationic here, selectivity)", "f")
    return A.save(fig, ctx.out("A6_positional.png"))


@figure("aa", "A7", "Residue pair co-occurrence")
def a7_cooc(ctx):
    C = ctx.counts_matrix(); sel = ctx.seqs.get("selectivity")
    pres = (C > 0).astype(float); m = _finite(sel)
    if m.sum() < 4:
        return None
    top = sel >= np.nanquantile(sel[m], 0.75); rest = m & ~top
    n = len(ctx.alphabet); E = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            both = pres[:, i] * pres[:, j]
            pt = both[top & m].mean() + 1e-3; pr = both[rest].mean() + 1e-3
            E[i, j] = np.log2(pt / pr)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.2)
    mx = np.nanmax(np.abs(E)) or 1.0
    im = ax.imshow(E, cmap=A.NPG_DIV, vmin=-mx, vmax=mx)
    ax.set_xticks(range(n)); ax.set_xticklabels(ctx.alphabet, fontsize=6)
    ax.set_yticks(range(n)); ax.set_yticklabels(ctx.alphabet, fontsize=6)
    A.colorbar(fig, ax, im, r"log$_2$ enrichment")
    A.finish(ax, "Residue pair co-occurrence (top-selective)", None, None, "g")
    return A.save(fig, ctx.out("A7_cooccurrence.png"))


# ========================================================================================= GA ==
@figure("ga", "E4", "GA fitness")
def e4_fitness(ctx):
    tl = ctx.timeline()
    if not tl:
        return None
    g = tl["generation"]
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.plot(g, tl["fitness_best"], color=A.NPG[3], label="best")
    ax.plot(g, tl["fitness_parents"], color=A.NPG[2], label="parents")
    ax.plot(g, tl["fitness_mean"], color=A.NPG[0], label="mean")
    ax.legend(loc="lower right")
    A.finish(ax, "GA fitness", "iteration", "fitness", "a")
    return A.save(fig, ctx.out("E4_fitness.png"))


@figure("ga", "E5", "Population turnover")
def e5_turnover(ctx):
    tl = ctx.timeline()
    if not tl:
        return None
    g = tl["generation"]
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.plot(g, tl["n_evaluated"], color=A.NPG[5], label="evaluated")
    ax.plot(g, tl["n_unique"], color=A.NPG[3], label="unique")
    ax.legend(loc="center right")
    A.finish(ax, "Population turnover", "iteration", "sequences", "b")
    return A.save(fig, ctx.out("E5_turnover.png"))


@figure("ga", "E6", "Diversity")
def e6_diversity(ctx):
    tl = ctx.timeline()
    if not tl:
        return None
    g = tl["generation"]; ent = tl["entropy_mean"]
    ent_norm = ent / np.nanmax(ent) if np.isfinite(ent).any() and np.nanmax(ent) > 0 else ent
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.plot(g, tl["pairwise_hamming"], color=A.NPG[3], label="pairwise dist.")
    ax.plot(g, ent_norm, color=A.NPG[0], label="entropy (norm.)")
    ax.legend(loc="upper right")
    A.finish(ax, "Diversity", "iteration", None, "c")
    return A.save(fig, ctx.out("E6_diversity.png"))


@figure("ga", "E7", "Composition by group")
def e7_composition(ctx):
    tl = ctx.timeline()
    if not tl:
        return None
    g = tl["generation"]; groups = list(mdtools.RESIDUE_GROUPS.keys())
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.0)
    for k, grp in enumerate(groups):
        col = tl.get("frac_%s" % grp)
        if col is not None:
            ax.plot(g, col, color=A.NPG[k % 8], label=grp)
    ax.legend(loc="center right", fontsize=5, ncol=1)
    A.finish(ax, "Composition by group", "iteration", "fraction", "d")
    return A.save(fig, ctx.out("E7_composition.png"))


@figure("ga", "E8", "Group frequency by position")
def e8_positional(ctx):
    pos = ctx.positional()
    if not pos:
        return None
    groups = list(dict.fromkeys(pos["group"]))
    positions = np.unique(pos["position"].astype(int))
    H = np.zeros((len(groups), len(positions)))
    gidx = {g: i for i, g in enumerate(groups)}
    pidx = {p: i for i, p in enumerate(positions)}
    for grp, p, f in zip(pos["group"], pos["position"], pos["frequency"]):
        H[gidx[grp], pidx[int(p)]] = f
    fig, ax = ctx.new(A.WIDTH_ONEHALF, 3.0)
    im = ax.imshow(H, aspect="auto", cmap=A.NPG_SEQ,
                   extent=[positions.min(), positions.max() + 1, len(groups), 0])
    ax.set_yticks(np.arange(len(groups)) + 0.5); ax.set_yticklabels(groups, fontsize=6)
    A.colorbar(fig, ax, im, "freq")
    A.finish(ax, "Group frequency by position", "position", None, "e")
    return A.save(fig, ctx.out("E8_positional.png"))


# =================================================================================== MEMBRANE ==
def _reference_values(ctx):
    """
    Reference membrane values (surften_free, area) per role, recovered from the saved cohort.

    These are campaign constants: every candidate on a given membrane shares the same reference
    area and free-membrane surface tension. surften_free is recovered as surften + delta_gamma
    (delta_gamma = gamma_free - surften), so no separate calibration import is needed and the
    figure stays a pure function of the saved dataset.
    """
    out = {}
    for role in ("target", "host"):
        area = ctx.col("area_%s" % role)
        surften = ctx.col("surften_%s" % role)
        dg = ctx.col("dg_%s" % role)
        if np.isfinite(area).any() and np.isfinite(surften).any() and np.isfinite(dg).any():
            out[role] = {"surften_free": float(np.nanmedian(surften + dg)),
                         "area": float(np.nanmedian(area))}
    return out


def _equil_series(ctx):
    """
    The saved equilibration series per role, as {role: {column: array}}.

    amp_dataset writes equil_<role>.csv (time_ns, surften, area) from the membrane build. Roles
    are read from metadata when present, otherwise by filename convention.
    """
    out = {}
    roles = ctx.meta.get("equilibration") or {role: "equil_%s.csv" % role
                                              for role in ("target", "host")}
    for role, fname in roles.items():
        data = ctx._load_optional(fname)
        if data and "time_ns" in data:
            out[role] = data
    return out


_ROLE_COLOR = {"target": A.TARGET, "host": A.HOST}
_ROLE_LABEL = {"target": "bacterial", "host": "host"}


@figure("membrane", "E1", "Equilibration: surface tension")
def e1_equil_gamma(ctx):
    series = _equil_series(ctx)
    series = {r: d for r, d in series.items() if "surften" in d}
    if not series:
        return None  # needs the build's equilibration series; skipped when absent
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.axhline(0, color=A.GREY, lw=0.6)
    for role, d in series.items():
        ax.plot(d["time_ns"], d["surften"], color=_ROLE_COLOR.get(role, A.HOST), lw=0.8,
                label=_ROLE_LABEL.get(role, role))
    if len(series) > 1:
        ax.legend(loc="upper right")
    A.finish(ax, "Equilibration: surface tension", "time (ns)", r"#Surf*SurfTen (bar nm)", "a")
    return A.save(fig, ctx.out("E1_equil_surften.png"))


@figure("membrane", "E2", "Equilibration: area")
def e2_equil_area(ctx):
    series = _equil_series(ctx)
    series = {r: d for r, d in series.items() if "area" in d}
    if not series:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    for role, d in series.items():
        ax.plot(d["time_ns"], d["area"], color=_ROLE_COLOR.get(role, A.HOST), lw=0.8,
                label=_ROLE_LABEL.get(role, role))
    if len(series) > 1:
        ax.legend(loc="upper right")
    A.finish(ax, "Equilibration: box area", "time (ns)", r"area (nm$^2$)", "b")
    return A.save(fig, ctx.out("E2_equil_area.png"))


@figure("membrane", "E3", "Membrane reference values")
def e3_reference(ctx):
    ref = _reference_values(ctx)
    if not ref or "target" not in ref or "host" not in ref:
        return None
    # Two metrics on very different scales, so each gets its own axis. Each is a one-row dumbbell:
    # the bacterial and host values joined, so the difference between membranes reads directly.
    panels = [("free surface tension", "surften_free", r"bar nm"),
              ("area per lipid", "area", r"nm$^2$")]
    fig, axes = plt.subplots(1, 2, figsize=(A.WIDTH_ONEHALF, 2.0))
    for k, (ax, (title, key, unit)) in enumerate(zip(axes, panels)):
        tv, hv = ref["target"][key], ref["host"][key]
        ax.plot([hv, tv], [0, 0], color=A.FAINT, lw=1.6, zorder=1)
        ax.plot(tv, 0, "o", color=A.TARGET, ms=7, label="bacterial", zorder=2)
        ax.plot(hv, 0, "o", color=A.HOST, ms=7, label="host", zorder=2)
        ax.set_yticks([]); ax.set_ylim(-1, 1); ax.set_xlabel(unit)
        ax.set_title(title)
        pad = 0.15 * (abs(tv - hv) or 1)
        ax.set_xlim(min(tv, hv) - pad - 0.05 * abs(min(tv, hv)), max(tv, hv) + pad)
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2)
    A.panel(axes[0], "c")
    fig.suptitle("Membrane reference values", x=0.01, ha="left", fontsize=8, fontweight="bold")
    return A.save(fig, ctx.out("E3_reference_values.png"))


# ================================================================================== CANDIDATE ==
def _target_membrane_type(ctx):
    for mtype, role in ctx.meta.get("membrane_role_map", {}).items():
        if role == "target":
            return mtype
    return None


def _cohort_row(ctx):
    ids = ctx.cohort.get("id")
    if ids is None:
        return None
    idx = np.where(ids.astype(str) == str(ctx.candidate))[0]
    return int(idx[0]) if idx.size else None


@figure("candidate", "P1", "Density profile")
def p1_density(ctx):
    s = ctx.series("density.csv")
    if not s or "z" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    colors = {"headgroups": A.NPG[3], "tails": A.NPG[2], "water": A.NPG[5], "peptide": A.NPG[0]}
    for name, c in colors.items():
        if name in s:
            ax.plot(s["z"], s[name], color=c, label=name)
    ax.legend(loc="upper right", fontsize=5)
    A.finish(ax, "Density profile along z", "z (nm)", "density (a.u.)", "a")
    return A.save(fig, ctx.out("P1_density.png"))


@figure("candidate", "P2", "Area-per-lipid map")
def p2_apl(ctx):
    d = ctx.candidate_dir(); path = os.path.join(d, "apl_map.txt") if d else None
    if not path or not os.path.isfile(path):
        return None
    grid = np.loadtxt(path)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    im = ax.imshow(grid, origin="lower", cmap=A.NPG_SEQ, aspect="equal")
    A.colorbar(fig, ax, im, r"APL (nm$^2$)")
    A.finish(ax, "Area-per-lipid map", "x (nm)", "y (nm)", "b")
    return A.save(fig, ctx.out("P2_apl_map.png"))


@figure("candidate", "P3", "Thickness map")
def p3_thickness(ctx):
    d = ctx.candidate_dir(); path = os.path.join(d, "thickness_map.txt") if d else None
    if not path or not os.path.isfile(path):
        return None
    grid = np.loadtxt(path)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    im = ax.imshow(grid, origin="lower", cmap=A.NPG_SEQ, aspect="equal")
    A.colorbar(fig, ax, im, "nm")
    A.finish(ax, "Thickness map", "x (nm)", "y (nm)", "c")
    return A.save(fig, ctx.out("P3_thickness_map.png"))


@figure("candidate", "P4", "Tail order parameter")
def p4_order(ctx):
    s = ctx.series("order.csv")
    if not s or "bead" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    if "far" in s:
        ax.plot(s["bead"], s["far"], "o-", color=A.NPG[3], label="far")
    if "under" in s:
        ax.plot(s["bead"], s["under"], "s-", color=A.TARGET, label="under peptide")
    ax.legend(loc="upper right")
    A.finish(ax, "Tail order parameter", "tail bead", "order", "d")
    return A.save(fig, ctx.out("P4_order.png"))


@figure("candidate", "P5", "Lateral pressure profile")
def p5_stress(ctx):
    s = ctx.series("stress.csv")
    if not s or "z" not in s or "pi" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.axhline(0, color=A.GREY, lw=0.6)
    ax.plot(s["z"], s["pi"], color=A.NPG[4])
    A.finish(ax, "Lateral pressure profile", "z (nm)", r"$\pi$(z) (bar)", "e")
    return A.save(fig, ctx.out("P5_stress.png"))


@figure("candidate", "P6", "Insertion depth over time")
def p6_depth(ctx):
    # Prefer the amp_traj series; fall back to the distance xvg the run already kept.
    s = ctx.series("depth.csv")
    t = z = None
    if s and "time" in s and "d_mem" in s:
        t, z = s["time"], s["d_mem"]
    else:
        mtype = _target_membrane_type(ctx); cdir = ctx.campaign_candidate_dir()
        if mtype and cdir:
            xvg = os.path.join(cdir, "production_%s_distance.xvg" % mtype)
            if os.path.isfile(xvg):
                x, series, _ = mdtools.read_xvg(xvg)
                t, z = x / 1000.0, next(iter(series.values()))  # ps -> ns
    if t is None:
        return None
    row = _cohort_row(ctx); b = ctx.col("bnm_target")[row] if row is not None else np.nan
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    if np.isfinite(b):
        ax.axhspan(0, 2 * b, color="#eef1f5", zorder=0); ax.axhline(b, color=A.TARGET, lw=1.0, ls="--")
    ax.plot(t, z, color=A.HOST)
    A.finish(ax, "Insertion depth over time", "time (ns)", r"d$_{mem}$ (nm)", "f")
    return A.save(fig, ctx.out("P6_depth.png"))


@figure("candidate", "P7", "Peptide tilt")
def p7_tilt(ctx):
    s = ctx.series("tilt.csv")
    if not s or "time" not in s or "tilt" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.plot(s["time"], s["tilt"], color=A.NPG[2])
    A.finish(ax, "Peptide tilt", "time (ns)", "tilt (deg)", "g")
    return A.save(fig, ctx.out("P7_tilt.png"))


@figure("candidate", "P8", "Per-residue contact map")
def p8_contacts(ctx):
    d = ctx.candidate_dir(); path = os.path.join(d, "contacts.txt") if d else None
    if not path or not os.path.isfile(path):
        return None
    occ = np.loadtxt(path)
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    im = ax.imshow(occ, origin="lower", aspect="auto", cmap=A.NPG_SEQ)
    A.colorbar(fig, ax, im, "contact")
    A.finish(ax, "Per-residue contact map", "frame", "residue", "h")
    return A.save(fig, ctx.out("P8_contacts.png"))


@figure("candidate", "P9", "Peptide COM vs midplane")
def p9_comz(ctx):
    s = ctx.series("comz.csv")
    if not s or "time" not in s or "z" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    ax.axhline(0, color=A.GREY, lw=0.6)
    ax.plot(s["time"], s["z"], color=A.TARGET)
    A.finish(ax, "Peptide COM vs midplane", "time (ns)", r"z - z$_{mid}$ (nm)", "i")
    return A.save(fig, ctx.out("P9_comz.png"))


@figure("candidate", "P10", "Interaction energy")
def p10_energy(ctx):
    s = ctx.series("energy.csv")
    if not s or "time" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    if "e_lipids" in s:
        ax.plot(s["time"], s["e_lipids"], color=A.NPG[3], lw=0.5, alpha=0.4)
    if "e_solvent" in s:
        ax.plot(s["time"], s["e_solvent"], color=A.NPG[0], lw=0.5, alpha=0.4)
    def run(y):
        return np.cumsum(y) / np.arange(1, len(y) + 1)
    if "e_lipids" in s:
        ax.plot(s["time"], run(s["e_lipids"]), color=A.NPG[3], label="pep-lipids")
    if "e_solvent" in s:
        ax.plot(s["time"], run(s["e_solvent"]), color=A.NPG[0], label="pep-solvent")
    ax.legend(loc="center right")
    A.finish(ax, "Interaction energy", "time (ns)", "E (kJ/mol)", "j")
    return A.save(fig, ctx.out("P10_energy.png"))


@figure("candidate", "P11", "Convergence")
def p11_conv(ctx):
    s = ctx.series("convergence.csv")
    if not s or "time" not in s or "cum" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    if "se" in s:
        ax.fill_between(s["time"], s["cum"] - s["se"], s["cum"] + s["se"],
                        color=A.NPG[3], alpha=0.18, lw=0)
    ax.plot(s["time"], s["cum"], color=A.NPG[3])
    A.finish(ax, r"$\delta\gamma$ convergence", "time (ns)", "running mean", "k")
    return A.save(fig, ctx.out("P11_convergence.png"))


@figure("candidate", "P12", "Target vs host (candidate)")
def p12_tgthost(ctx):
    row = _cohort_row(ctx)
    if row is None:
        return None
    def val(col):
        v = ctx.col(col); return v[row] if row < len(v) else np.nan
    metrics, tgt, hst = [], [], []
    pairs = [("dg", "dg_target", "dg_host"), ("orient", "orient_target", "orient_host"),
             ("|E_lip|/100", "elip_target", "elip_host")]
    for lab, tcol, hcol in pairs:
        tv, hv = val(tcol), val(hcol)
        if lab.startswith("|E"):
            tv, hv = abs(tv) / 100.0, abs(hv) / 100.0
        if np.isfinite(tv) and np.isfinite(hv):
            metrics.append(lab); tgt.append(tv); hst.append(hv)
    if not metrics:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 2.6)
    A.dumbbell(ax, metrics, hst, tgt, sort_by_gap=False)
    ax.legend(loc="lower right")
    A.finish(ax, "Target vs host (candidate)", "value", None, "l")
    return A.save(fig, ctx.out("P12_target_host.png"))


@figure("candidate", "P13", "Binding PMF")
def p13_pmf(ctx):
    s = ctx.series("pmf.csv")
    if not s or "xi" not in s:
        return None
    fig, ax = ctx.new(A.WIDTH_SINGLE, 3.0)
    for role, col, c in [("bacterial", "target", A.TARGET), ("host", "host", A.HOST)]:
        if col in s:
            ax.plot(s["xi"], s[col], color=c, label=role)
            if col + "_err" in s:
                ax.fill_between(s["xi"], s[col] - s[col + "_err"], s[col] + s[col + "_err"],
                                color=c, alpha=0.15, lw=0)
    ax.legend(loc="lower right")
    A.finish(ax, "Binding PMF", "distance (nm)", "PMF (kJ/mol)", "m")
    return A.save(fig, ctx.out("P13_pmf.png"))


# ======================================================================================= CLI ==
def _contact_sheet(paths, out, title, cols=4):
    import matplotlib.image as mpimg
    paths = [p for p in paths if p]
    rows = int(np.ceil(len(paths) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.0, rows * 2.6),
                             constrained_layout=False)
    axes = np.atleast_1d(axes).ravel()
    for k, ax in enumerate(axes):
        ax.axis("off")
        if k < len(paths):
            ax.imshow(mpimg.imread(paths[k]))
            ax.set_title(os.path.basename(paths[k]).split(".")[0], fontsize=8, loc="left")
    fig.suptitle(title, x=0.01, ha="left", fontsize=13, fontweight="bold")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01, wspace=0.05, hspace=0.25)
    fig.savefig(out, dpi=130); plt.close(fig)
    return out


def main(argv=None):
    groups = sorted({g for g, *_ in _FIGURES})
    parser = argparse.ArgumentParser(description="Render AMP-evo-MD analysis figures from the saved dataset.")
    parser.add_argument("--dataset", required=True, help="analysis dir written by amp_dataset.py")
    parser.add_argument("--outdir", default=None, help="where to write figures (default <dataset>/figures)")
    parser.add_argument("--group", nargs="+", default=["all"], choices=groups + ["all"])
    parser.add_argument("--candidate", default=None,
                        help="candidate id (sequence) for the per-candidate (P) figures")
    parser.add_argument("--sheet", action="store_true", help="also write a contact sheet per group")
    args = parser.parse_args(argv)

    outdir = args.outdir or os.path.join(args.dataset, "figures")
    os.makedirs(outdir, exist_ok=True)
    ctx = Context(args.dataset, outdir)
    ctx.candidate = args.candidate
    want = set(groups) if "all" in args.group else set(args.group)
    if "candidate" in want and not ctx.candidate:
        sys.stderr.write("Note: the per-candidate (P) figures need --candidate <id>; skipping them.\n")
        want.discard("candidate")

    written = {}
    for group, code, title, fn in _FIGURES:
        if group not in want:
            continue
        try:
            path = fn(ctx)
        except Exception as exc:  # a single figure failing must not lose the rest
            sys.stderr.write("  [skip] %s %s: %s\n" % (code, title, exc))
            continue
        if path:
            written.setdefault(group, []).append(path)
            print("  %s %-34s %s" % (code, title, os.path.basename(path)))

    if args.sheet:
        for group, paths in written.items():
            sheet = _contact_sheet(paths, os.path.join(outdir, "sheet_%s.png" % group),
                                   "AMP-evo-MD: %s figures" % group)
            print("  sheet  %s" % os.path.basename(sheet))
    total = sum(len(v) for v in written.values())
    print("Rendered %d figures to %s" % (total, outdir))


if __name__ == "__main__":
    main()
