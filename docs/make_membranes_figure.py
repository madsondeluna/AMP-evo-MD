# Build docs/membranes.png: the three lipids shown atom-by-atom and as Martini 3 CG beads, side by
# side and in the same colours so the mapping is readable, plus the two membranes' topology and
# charge balance. Needs only matplotlib. Regenerate with: python docs/make_membranes_figure.py
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

# ---- shared palette (all-atom atoms and CG beads use the same colours) ----------------
C_POS   = "#4C72B0"   # cationic head (NC3, NH3)  -> +1
C_NEG   = "#C44E52"   # phosphate (PO4)           -> -1
C_NEUT  = "#8172B3"   # neutral polar head (GL0)  ->  0
C_GLY   = "#55A868"   # glycerol / ester linkers (GL1, GL2)
C_SAT   = "#B8B8B8"   # saturated tail bead (C1)
C_UNSAT = "#DD8452"   # unsaturated tail bead (C4h / D, the double bond)

HEADNAME  = {"POPC": "choline", "POPE": "ethanolamine", "POPG": "glycerol"}
HEADBEAD  = {"POPC": ("NC3", +1), "POPE": ("NH3", +1), "POPG": ("GL0", 0)}
NETCHARGE = {"POPC": 0, "POPE": 0, "POPG": -1}

# ---- all-atom hairpin: every heavy atom as a circle, coloured by the CG bead it maps to -------
def draw_aa(ax, lip):
    A, B = [], []
    def add(x, y, c):
        A.append((x, y, c)); return len(A) - 1
    def link(i, j):
        B.append((i, j))

    def tail(x0, ncarb, ytop, ystep=0.27, unsat=None):
        idx = []
        for i in range(ncarb):
            xo = 0.11 * (-1 if i % 2 else 1)
            col = C_UNSAT if (unsat and unsat[0] <= i <= unsat[1]) else C_SAT
            k = add(x0 + xo, ytop - i * ystep, col)
            if idx:
                link(idx[-1], k)
            idx.append(k)
        return idx

    # two tails: A = sn-2 oleoyl 18:1 (kink, orange beads), B = sn-1 palmitoyl 16:0 (saturated)
    tA = tail(-0.5, 15, 4.35, unsat=(7, 8))
    tB = tail(0.5, 14, 4.35)

    # glycerol + ester carbonyls (green)
    gC = add(0.0, 4.95, C_GLY)
    gA = add(-0.5, 4.72, C_GLY)
    gB = add(0.5, 4.72, C_GLY)
    link(gC, gA); link(gC, gB); link(gA, tA[0]); link(gB, tB[0])

    # phosphate (red): P and two O
    p = add(0.0, 5.55, C_NEG)
    link(gC, p)
    link(p, add(-0.24, 5.78, C_NEG))
    link(p, add(0.24, 5.78, C_NEG))

    # head
    if lip in ("POPE", "POPC"):
        hc = C_POS
        c1 = add(0.0, 6.15, hc); c2 = add(0.0, 6.55, hc); n = add(0.0, 7.00, hc)
        link(p, c1); link(c1, c2); link(c2, n)
        if lip == "POPC":                       # choline trimethyl
            for dx in (-0.30, 0.30, 0.0):
                link(n, add(dx, 7.28, hc))
    else:                                        # POPG neutral glycerol head (purple)
        hc = C_NEUT
        c1 = add(0.0, 6.15, hc)
        c2 = add(-0.26, 6.55, hc); c3 = add(0.26, 6.55, hc)
        link(p, c1); link(c1, c2); link(c1, c3)
        link(c2, add(-0.46, 6.82, hc)); link(c3, add(0.46, 6.82, hc))

    for i, j in B:
        ax.plot([A[i][0], A[j][0]], [A[i][1], A[j][1]], color="#8a8a8a", lw=1.4, zorder=1)
    for x, y, c in A:
        ax.add_patch(Circle((x, y), 0.115, facecolor=c, edgecolor="#4a4a4a", lw=0.5, zorder=2))
    ax.set_xlim(-1.15, 1.15); ax.set_ylim(0.1, 7.7)
    ax.set_aspect("equal"); ax.axis("off")

# ---- Martini 3 CG bead model ----------------------------------------------------------
def cg_beads(headbead, headq, headcolour):
    return [
        (headbead, 0.0, 6.0, headcolour, headq),
        ("PO4",    0.0, 5.1, C_NEG, -1),
        ("GL1",   -0.45, 4.2, C_GLY, 0),
        ("GL2",    0.45, 4.2, C_GLY, 0),
        ("C1A",   -0.45, 3.3, C_SAT, 0),
        ("D2A",   -0.45, 2.4, C_UNSAT, 0),
        ("C3A",   -0.45, 1.5, C_SAT, 0),
        ("C4A",   -0.45, 0.6, C_SAT, 0),
        ("C1B",    0.45, 3.3, C_SAT, 0),
        ("C2B",    0.45, 2.4, C_SAT, 0),
        ("C3B",    0.45, 1.5, C_SAT, 0),
        ("C4B",    0.45, 0.6, C_SAT, 0),
    ]

BONDS = [(0,1),(1,2),(1,3),(2,4),(4,5),(5,6),(6,7),(3,8),(8,9),(9,10),(10,11)]

def draw_cg(ax, lip):
    headbead, headq = HEADBEAD[lip]
    headcolour = C_POS if headq > 0 else C_NEUT
    beads = cg_beads(headbead, headq, headcolour)
    for i, j in BONDS:
        ax.plot([beads[i][1], beads[j][1]], [beads[i][2], beads[j][2]],
                color="#9a9a9a", lw=1.3, zorder=1)
    for label, x, y, col, q in beads:
        ax.add_patch(Circle((x, y), 0.33, facecolor=col, edgecolor="#3a3a3a", lw=0.8, zorder=2))
        ax.text(x, y, label, ha="center", va="center", fontsize=6.3,
                color="white", fontweight="bold", zorder=3)
        if q > 0:
            ax.text(x+0.34, y+0.3, "+", ha="center", va="center", fontsize=12,
                    color=C_POS, fontweight="bold", zorder=4)
        elif q < 0:
            ax.text(x+0.34, y+0.3, "−", ha="center", va="center", fontsize=12,
                    color=C_NEG, fontweight="bold", zorder=4)
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(0.0, 6.7)
    ax.set_aspect("equal"); ax.axis("off")

# ---- figure ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "DejaVu Sans"
fig = plt.figure(figsize=(13.6, 11.8))
fig.patch.set_facecolor("white")

fig.text(0.5, 0.977, "Membranes: lipids, coarse-graining and charge",
         ha="center", fontsize=18, fontweight="bold", color="#1a1a1a")
fig.text(0.5, 0.957,
         "Martini 3, non-polarisable water. The three lipids share tails, glycerol and phosphate; "
         "only the first head bead differs, and that sets the charge.",
         ha="center", fontsize=10.5, color="#555555")

lipids = ["POPE", "POPG", "POPC"]
left, colw, gap = 0.050, 0.308, 0.006
y_panel, h_panel, w_aa = 0.515, 0.370, 0.150

for k, lip in enumerate(lipids):
    x0 = left + k * (colw + gap)
    anion = " (anionic, net −1)" if NETCHARGE[lip] < 0 else " (zwitterionic, net 0)"
    fig.text(x0 + colw/2, 0.920, lip + anion, ha="center", fontsize=13, fontweight="bold",
             color=("#C44E52" if NETCHARGE[lip] < 0 else "#333333"))

    ax_aa = fig.add_axes([x0, y_panel, w_aa, h_panel])
    draw_aa(ax_aa, lip)
    fig.text(x0 + w_aa/2, 0.897, "all-atom", ha="center", fontsize=8.5, color="#999999")

    x_cg = x0 + w_aa + 0.012
    w_cg = colw - w_aa - 0.012
    ax_cg = fig.add_axes([x_cg, y_panel, w_cg, h_panel])
    draw_cg(ax_cg, lip)
    fig.text(x_cg + w_cg/2, 0.897, "Martini 3 CG, 12 beads", ha="center",
             fontsize=8.5, color="#999999")

    hb, hq = HEADBEAD[lip]
    fig.text(x0 + colw/2, 0.497,
             "head: %s  →  %s(%s) + PO4(−1)" %
             (HEADNAME[lip], hb, ("+1" if hq > 0 else "0")),
             ha="center", fontsize=9, color="#444444")

# bead legend
lx, ly = 0.050, 0.440
fig.text(lx, ly + 0.022, "Shared colours (atoms and beads)", fontsize=10.5,
         fontweight="bold", color="#1a1a1a")
items = [("cationic head (+1)", C_POS), ("neutral head (0)", C_NEUT),
         ("phosphate (−1)", C_NEG), ("glycerol / ester", C_GLY),
         ("saturated tail", C_SAT), ("unsaturated tail (double bond)", C_UNSAT)]
for i, (lab, col) in enumerate(items):
    xx = lx + (i % 3) * 0.235
    yy = ly - (i // 3) * 0.024
    fig.patches.append(Circle((xx, yy), 0.007, transform=fig.transFigure,
                             facecolor=col, edgecolor="#3a3a3a", lw=0.6))
    fig.text(xx + 0.014, yy, lab, va="center", fontsize=9, color="#444444")

# ---- membranes panel ------------------------------------------------------------------
def membrane_panel(ax, title, subtitle, lines, charge_lines, charge_colour):
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.01, 0.01), 0.98, 0.98,
                 boxstyle="round,pad=0.01,rounding_size=0.03",
                 facecolor="#fbfbfa", edgecolor="#d9d9d4", lw=1.2, transform=ax.transAxes))
    ax.text(0.05, 0.915, title, fontsize=12.5, fontweight="bold",
            transform=ax.transAxes, color="#1a1a1a")
    ax.text(0.05, 0.845, subtitle, fontsize=9, transform=ax.transAxes, color="#777")
    y = 0.72
    for lab, val in lines:
        ax.text(0.07, y, lab, fontsize=10, transform=ax.transAxes, color="#444")
        ax.text(0.93, y, val, fontsize=10, transform=ax.transAxes, color="#222",
                ha="right", fontfamily="DejaVu Sans Mono")
        y -= 0.10
    y -= 0.035
    for cl in charge_lines:
        ax.text(0.05, y, cl, fontsize=9.5, transform=ax.transAxes,
                color=charge_colour, fontweight="bold")
        y -= 0.082

fig.text(0.5, 0.380, "The two membranes: composition and charge balance",
         ha="center", fontsize=14, fontweight="bold", color="#1a1a1a")

axb = fig.add_axes([0.050, 0.048, 0.44, 0.312])
membrane_panel(
    axb, "Bacterial target  —  POPE:POPG 3:1",
    "Gram-negative inner-membrane mimic. 198 lipids, ~5285 beads, 0.15 M NaCl.",
    [("POPE (zwitterionic)", "140"), ("POPG (anionic, −1 each)", "58"),
     ("W  (Martini water)", "2847"), ("Na⁺  /  Cl⁻", "60 / 2")],
    ["Membrane charge −58 e  (58 × POPG).",
     "Excess Na⁺ neutralises it → system neutral.",
     "Anionic surface."], "#C44E52")

axh = fig.add_axes([0.515, 0.048, 0.44, 0.312])
membrane_panel(
    axh, "Host  —  POPC (neutral)",
    "Mammalian-membrane mimic. Cholesterol and sphingomyelin deferred.",
    [("POPC (zwitterionic)", "194"), ("W  (Martini water)", "~2716"),
     ("Na⁺  /  Cl⁻", "equal (salt only)")],
    ["Membrane charge 0  (zwitterionic POPC).",
     "Na⁺ = Cl⁻ : the 0.15 M salt only.",
     "Neutral surface."], "#2E6FAE")

fig.text(0.5, 0.028,
         "Selectivity contrast: the anionic bacterial surface (POPG) versus the neutral host (POPC) "
         "is what a cationic antimicrobial peptide discriminates on.",
         ha="center", fontsize=10, color="#555555", style="italic")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "membranes.png")
fig.savefig(out, dpi=200, facecolor="white")
print("wrote", out)
