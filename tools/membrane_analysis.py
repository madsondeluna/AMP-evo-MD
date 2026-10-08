#!/usr/bin/env python3
"""
Analyse a peptide-free membrane equilibration and compare membranes against each other.

Reads the .xvg files that gmx writes, so nothing here needs GROMACS installed or a trajectory in
memory. tools/protocol.sh produces the inputs; this script turns them into the equilibration
verdict, the figures, and the four reference constants that user/membranes.py requires.

    python tools/membrane_analysis.py single --name bacterial --dir prep/bacterial
    python tools/membrane_analysis.py compare --dir prep --names bacterial host

WHAT IS BEING CHECKED
---------------------
Three things decide whether a membrane is ready to be used as a reference, and all three can
fail quietly:

  Area has stopped drifting. The constant-area production box must be set to the composition's
  equilibrium area. If the semi-isotropic stage was stopped while the area was still relaxing,
  the production run measures a pre-stressed bilayer and every peptide is then evaluated against
  the wrong baseline.

  Free-membrane surface tension is zero within error. This is the direct test of the above. The
  Evo-MD POPC reference sits at 1.93 +/- 0.99 bar nm. A value far from zero means the box area and
  the equilibrium area disagree, and the fix is to re-plan the lipid count, not to proceed.

  The leaflets are symmetric and the membrane has not drifted across the periodic boundary. An
  area mismatch between leaflets appears as spontaneous curvature and contaminates the surface
  tension; a membrane split across the boundary silently corrupts any thickness measurement.
"""

import argparse
import collections
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mdtools
import membrane


# gmx energy writes these with its own labels; match loosely so a different gmx version that
# renames a term is still picked up rather than silently yielding an empty series.
BOX_X_KEYS = ("Box-X", "Box-X ")
BOX_Y_KEYS = ("Box-Y",)
SURFTEN_KEYS = ("#Surf*SurfTen", "Surf*SurfTen", "SurfTen")

# The fraction of a trajectory discarded before averaging, taken from the pipeline settings so
# this tool, tools/membrane.py and the production runs all use one window. delta_gamma is the
# difference between a reference and a per-candidate value, and averaging the two over different
# windows would leave any residual drift inside the observable.
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import user.usersettings as _UserSettings
    DEFAULT_EQUIL_FRACTION = 1.0 - _UserSettings.ENERGY_init_fraction
except ImportError as _exc:
    # Runnable against copied-off files without the package importable. Narrowed from a bare
    # except, which also caught a syntax error or a failed assertion inside usersettings.py and
    # then silently averaged this reference over a different window than the candidates use,
    # putting the residual drift straight into delta_gamma with nothing saying so.
    DEFAULT_EQUIL_FRACTION = 0.05
    sys.stderr.write(
        "note: user.usersettings is not importable (%s), so the equilibration fraction falls "
        "back to %.3f. If the pipeline uses a different ENERGY_init_fraction, this reference "
        "and the candidates will be averaged over different windows; pass --equil-fraction to "
        "match them.\n" % (_exc, DEFAULT_EQUIL_FRACTION))


def pick_series(series, keys):
    for key in keys:
        for label, values in series.items():
            if label.strip().startswith(key):
                return values
    return None


def analyse_one(name, directory, equil_fraction=None):
    """Measure one membrane from the files protocol.sh leaves in its directory."""
    if equil_fraction is None:
        equil_fraction = DEFAULT_EQUIL_FRACTION
    result = {"name": name, "dir": directory}

    energy_xvg = os.path.join(directory, "%s_energy.xvg" % name)
    if not os.path.exists(energy_xvg):
        raise SystemExit(
            "Missing %s. Produce it with:\n"
            "    printf 'Box-X\\nBox-Y\\n#Surf*SurfTen\\n' | gmx energy -f %s.edr -o %s_energy.xvg"
            % (energy_xvg, os.path.join(directory, name), os.path.join(directory, name))
        )

    time, series, _ = mdtools.read_xvg(energy_xvg)
    box_x = pick_series(series, BOX_X_KEYS)
    box_y = pick_series(series, BOX_Y_KEYS)
    gamma = pick_series(series, SURFTEN_KEYS)
    if gamma is None:
        raise SystemExit("No surface tension column in %s (found: %s)"
                         % (energy_xvg, ", ".join(series)))

    area = box_x * box_y if box_x is not None and box_y is not None else None

    t_eq, gamma_eq = mdtools.tail(time, gamma, equil_fraction)
    gamma_mean, gamma_err = mdtools.block_average(gamma_eq)
    # Block-averaged, which is the right estimator for a correlated series, but NOT the one that
    # ends up in the MembraneDef: tools/membrane.py takes gmx energy's own autocorrelation-based
    # Err.Est., and the production runs take the same. Mixing the two would make delta_gamma_err
    # a quadrature sum of two different estimators. This value is a diagnostic for the
    # equilibration verdict below; the canonical reference constant comes from
    # `tools/membrane.py calibrate`.
    # bar*nm: the unit of the #Surf*SurfTen column, which is the surface tension times
    # the number of interfaces. Not bar.
    result["surften"] = {"mean": gamma_mean, "err": gamma_err, "units": "bar*nm"}
    result["surften_drift"] = mdtools.drift(t_eq, gamma_eq)
    result["time"] = {"start": float(time[0]), "end": float(time[-1]),
                      "analysed_from": float(t_eq[0])}

    if area is not None:
        t_eq_a, area_eq = mdtools.tail(time, area, equil_fraction)
        area_mean, area_err = mdtools.block_average(area_eq)
        # Time average over the equilibration, which is a drift diagnostic. It is NOT the area
        # the production run is held at: that is box.x * box.y of the final structure, which is
        # what tools/membrane.py calibrate reads and writes into the MembraneDef.
        result["area"] = {"mean": area_mean, "err": area_err, "units": "nm2",
                          "note": "time average, diagnostic only; the held area is box.x*box.y"}
        # Drift expressed over the analysed window, which is what matters: a slope that is
        # negligible per picosecond can still be a percent of the area over 100 ns.
        slope = mdtools.drift(t_eq_a, area_eq)
        span = float(t_eq_a[-1] - t_eq_a[0])
        result["area_drift_pct_over_window"] = 100.0 * slope * span / area_mean if area_mean else float("nan")

    gro = os.path.join(directory, "%s_eq.gro" % name)
    if os.path.exists(gro):
        bx, by, bz, atoms = mdtools.read_gro(gro)
        result["box"] = {"x": bx, "y": by, "z": bz}

        # Counting goes through membrane.count_lipids, not a local bead match. Cardiolipin
        # carries PO41 and PO42 rather than PO4 and sterols carry neither, so a direct name
        # match undercounts the first and drops the second, inflating the area per lipid and
        # the composition line this writes into the report.
        flat = [(res, atom, z) for res, atom, _x, _y, z in atoms]
        counts, _zs = membrane.count_lipids(flat)
        result["lipids"] = dict(counts)
        n_lipids = sum(counts.values())
        if n_lipids and area is not None:
            result["area_per_lipid_nm2"] = result["area"]["mean"] / (n_lipids / 2.0)

        # Thickness and the leaflet balance are measured on phosphates only, which is the right
        # reference for insertion depth; sterols are deliberately absent from this count, so the
        # two numbers do not have to agree with n_lipids.
        try:
            thickness, n_upper, n_lower = membrane.measure_thickness(flat)
            result["thickness_nm"] = thickness
            result["leaflets"] = {"upper": n_upper, "lower": n_lower}
        except SystemExit:
            # No phosphate beads at all, which a pure sterol or ether-lipid membrane could
            # produce. Report the rest rather than aborting the whole analysis.
            result["thickness_nm"] = float("nan")

    result["_series"] = {"time": time, "gamma": gamma,
                         "area": area if area is not None else None}
    return result


def verdict(result):
    """Plain-language readings, because each failure here is silent in the numbers alone."""
    notes = []
    g = result["surften"]
    if abs(g["mean"]) > 10 and abs(g["mean"]) > 3 * max(g["err"], 1e-9):
        notes.append(
            "Surface tension is %.1f +/- %.1f bar nm, not zero. The constant-area box does not "
            "match this composition's equilibrium area, so peptides would be evaluated against "
            "a pre-stressed bilayer. Re-plan the lipid count and rebuild."
            % (g["mean"], g["err"]))
    else:
        notes.append("Surface tension is %.1f +/- %.1f bar nm, zero within error. Ready to use as a "
                     "reference." % (g["mean"], g["err"]))

    if "area_drift_pct_over_window" in result:
        d = result["area_drift_pct_over_window"]
        if abs(d) > 1.0:
            notes.append("Area still drifting: %.2f %% over the analysed window. Equilibrate "
                         "longer before trusting the reference area." % d)
        else:
            notes.append("Area drift %.2f %% over the analysed window, flat enough." % d)

    if "leaflets" in result:
        u, l = result["leaflets"]["upper"], result["leaflets"]["lower"]
        if u != l:
            notes.append("Leaflets asymmetric (%d vs %d phosphate beads). Intentional only for "
                         "an outer-membrane model; otherwise this shows up as spontaneous "
                         "curvature and contaminates the surface tension. Note that sterols are "
                         "not counted here, so a sterol-rich membrane can be balanced overall "
                         "while these two differ." % (u, l))
    return notes


def figure_single(result, outpath):
    plt = mdtools.apply_style()
    s = result["_series"]
    panels = 2 if s["area"] is not None else 1
    fig, axes = plt.subplots(panels, 1, figsize=(7.5, 2.6 * panels), sharex=True)
    axes = np.atleast_1d(axes)

    start = result["time"]["analysed_from"]

    ax = axes[0]
    ax.plot(s["time"], s["gamma"], color=mdtools.CATEGORICAL[0], lw=0.8)
    ax.axhline(0, color=mdtools.GRID, lw=0.8)
    ax.axvspan(s["time"][0], start, color=mdtools.GRID, alpha=0.3, lw=0)
    ax.set_ylabel("#Surf*SurfTen (bar nm)")
    ax.set_title("%s: free membrane" % result["name"])
    ax.annotate("discarded", xy=(start, ax.get_ylim()[1]), xytext=(-4, -10),
                textcoords="offset points", ha="right", fontsize=7,
                color=mdtools.INK_MUTED)
    mdtools.declutter(ax)

    if s["area"] is not None:
        ax = axes[1]
        ax.plot(s["time"], s["area"], color=mdtools.CATEGORICAL[0], lw=0.8)
        ax.axvspan(s["time"][0], start, color=mdtools.GRID, alpha=0.3, lw=0)
        ax.set_ylabel("area (nm2)")
        mdtools.declutter(ax)
    axes[-1].set_xlabel("time (ps)")

    fig.tight_layout(h_pad=2.0)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def figure_compare(results, outpath):
    """
    One panel per quantity, membranes as the categorical series.

    Capped at four membranes: a fifth would need a generated colour, and the honest alternative
    is to facet rather than invent a hue.
    """
    plt = mdtools.apply_style()
    if len(results) > len(mdtools.CATEGORICAL):
        raise SystemExit("compare handles at most %d membranes at once; run it in groups."
                         % len(mdtools.CATEGORICAL))

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))

    ax = axes[0]
    for colour, r in zip(mdtools.CATEGORICAL, results):
        s = r["_series"]
        ax.plot(s["time"], s["gamma"], color=colour, lw=0.8)
        mdtools.label_last(ax, s["time"], s["gamma"], r["name"], colour)
    ax.axhline(0, color=mdtools.GRID, lw=0.8)
    ax.set_title("#Surf*SurfTen (bar nm)")
    ax.set_xlabel("time (ps)")
    mdtools.declutter(ax)

    # Dot with interval, not a bar. The question this panel answers is whether the tension is
    # zero within error, so the reference line at zero and the span of the interval are the
    # content; a bar would encode the value as area and bury the interval that matters more.
    ax = axes[1]
    names = [r["name"] for r in results]
    values = [r["surften"]["mean"] for r in results]
    errors = [r["surften"]["err"] for r in results]
    ypos = np.arange(len(results))[::-1]
    ax.axvline(0, color=mdtools.INK_MUTED, lw=0.8)
    for y, v, e, colour in zip(ypos, values, errors, mdtools.CATEGORICAL):
        ax.plot([v - e, v + e], [y, y], color=colour, lw=1.6, solid_capstyle="round")
        ax.plot([v], [y], "o", color=colour, markersize=6,
                markeredgecolor=mdtools.SURFACE, markeredgewidth=1.5)
        ax.annotate("%.1f +/- %.1f" % (v, e), xy=(v + e, y), xytext=(8, 0),
                    textcoords="offset points", va="center", fontsize=7.5,
                    color=mdtools.INK)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names)
    ax.set_ylim(-0.6, len(results) - 0.4)
    ax.set_title("Equilibrium #Surf*SurfTen (bar nm)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    span = max(abs(v) + e for v, e in zip(values, errors)) * 2.2
    ax.set_xlim(-span, span)
    # Only the x axis is decluttered here: the y axis carries the membrane names, and letting a
    # locator near them drops rows.
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.spines["bottom"].set_position(("outward", 6))

    ax = axes[2]
    apl = [r.get("area_per_lipid_nm2", float("nan")) for r in results]
    ax.bar(names, apl, color=mdtools.CATEGORICAL[:len(results)], width=0.5)
    ax.set_title("Area per lipid (nm2)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelleft=False)
    for i, v in enumerate(apl):
        if np.isfinite(v):
            ax.annotate("%.3f" % v, xy=(i, v), xytext=(0, 6), textcoords="offset points",
                        ha="center", fontsize=8, color=mdtools.INK)

    fig.tight_layout(w_pad=3.0)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def print_table(results):
    header = ("membrane", "area nm2", "APL nm2", "thickness", "surf bar nm", "err", "lipids")
    print("  " + "  ".join("%11s" % h for h in header))
    for r in results:
        print("  " + "  ".join([
            "%11s" % r["name"][:11],
            "%11.3f" % r.get("area", {}).get("mean", float("nan")),
            "%11.4f" % r.get("area_per_lipid_nm2", float("nan")),
            "%11.4f" % r.get("thickness_nm", float("nan")),
            "%11.2f" % r["surften"]["mean"],
            "%11.2f" % r["surften"]["err"],
            "%11d" % sum(r.get("lipids", {}).values()),
        ]))


def strip_series(result):
    out = {k: v for k, v in result.items() if k != "_series"}
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    s = sub.add_parser("single", help="analyse one membrane")
    s.add_argument("--name", required=True)
    s.add_argument("--dir", required=True)
    s.add_argument("--equil-fraction", type=float, default=DEFAULT_EQUIL_FRACTION,
                   help="fraction of the trajectory kept for averaging (default %.2f, taken from "
                        "UserSettings.ENERGY_init_fraction so this agrees with the window the "
                        "production runs and tools/membrane.py use)" % DEFAULT_EQUIL_FRACTION)
    s.add_argument("--outdir", default=None)

    c = sub.add_parser("compare", help="compare several membranes")
    c.add_argument("--dir", required=True, help="parent directory holding one subdir per membrane")
    c.add_argument("--names", nargs="+", required=True)
    c.add_argument("--equil-fraction", type=float, default=DEFAULT_EQUIL_FRACTION)
    c.add_argument("--outdir", default=None)

    args = parser.parse_args()

    if args.command == "single":
        outdir = args.outdir or args.dir
        os.makedirs(outdir, exist_ok=True)
        result = analyse_one(args.name, args.dir, args.equil_fraction)
        print("MEMBRANE %s" % result["name"])
        print_table([result])
        print()
        for note in verdict(result):
            print("  " + note)
        figure_single(result, os.path.join(outdir, "%s_equilibration.png" % args.name))
        with open(os.path.join(outdir, "%s_analysis.json" % args.name), "w") as handle:
            json.dump(strip_series(result), handle, indent=2)
        print()
        print("Wrote %s_equilibration.png and %s_analysis.json to %s"
              % (args.name, args.name, outdir))

    else:
        outdir = args.outdir or args.dir
        os.makedirs(outdir, exist_ok=True)
        results = [analyse_one(n, os.path.join(args.dir, n), args.equil_fraction)
                   for n in args.names]
        print("MEMBRANE COMPARISON")
        print_table(results)
        print()
        for r in results:
            print("  %s:" % r["name"])
            for note in verdict(r):
                print("    " + note)
        areas = [r.get("area", {}).get("mean") for r in results]
        if all(a for a in areas):
            spread = 100 * (max(areas) - min(areas)) / min(areas)
            print()
            print("  Area spread across membranes: %.1f %%." % spread)
            if spread > 10:
                print("  Over 10 %%. Delta gamma is per unit area, so a spread this large means the")
                print("  per-peptide perturbation is diluted noticeably differently in each")
                print("  membrane and the comparison between them carries a correction nobody")
                print("  wants to defend. Adjust the lipid counts so the equilibrium areas match.")
        figure_compare(results, os.path.join(outdir, "membrane_comparison.png"))
        with open(os.path.join(outdir, "membrane_comparison.json"), "w") as handle:
            json.dump([strip_series(r) for r in results], handle, indent=2)
        print()
        print("Wrote membrane_comparison.png and membrane_comparison.json to %s" % outdir)


if __name__ == "__main__":
    main()
