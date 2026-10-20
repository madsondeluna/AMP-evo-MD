# Build docs/membranes.png: the three lipids shown atom by atom and as Martini 3 CG beads, in the
# same colours, plus the two membranes' composition and charge balance. Needs only matplotlib.
# Regenerate with: python docs/make_membranes_figure.py
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.lines import Line2D

# ---- shared palette (all-atom atoms and CG beads use the same colours) ----------------
C_POS   = "#4C72B0"   # cationic head (NC3, NH3)  -> +1
C_NEG   = "#C44E52"   # phosphate (PO4)           -> -1
C_NEUT  = "#8172B3"   # neutral polar head (GL0)  ->  0
C_GLY   = "#55A868"   # glycerol / ester linkers (GL1, GL2)
C_SAT   = "#B8B8B8"   # saturated tail bead (C1)
C_UNSAT = "#DD8452"   # unsaturated tail bead (C4h / D, the double bond)
INK     = "#1a1a1a"
MUTE    = "#6b6b6b"
FAINT   = "#9a9a9a"

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

    tA = tail(-0.5, 15, 4.35, unsat=(7, 8))   # sn-2 oleoyl 18:1 (kink, orange beads)
    tB = tail(0.5, 14, 4.35)                   # sn-1 palmitoyl 16:0 (saturated)

    gC = add(0.0, 4.95, C_GLY)
    gA = add(-0.5, 4.72, C_GLY)
    gB = add(0.5, 4.72, C_GLY)
    link(gC, gA); link(gC, gB); link(gA, tA[0]); link(gB, tB[0])

    p = add(0.0, 5.55, C_NEG)
    link(gC, p)
    link(p, add(-0.24, 5.78, C_NEG))
    link(p, add(0.24, 5.78, C_NEG))

    if lip in ("POPE", "POPC"):
        hc = C_POS
        c1 = add(0.0, 6.15, hc); c2 = add(0.0, 6.55, hc); n = add(0.0, 7.00, hc)
        link(p, c1); link(c1, c2); link(c2, n)
        if lip == "POPC":
            for dx in (-0.30, 0.30, 0.0):
                link(n, add(dx, 7.28, hc))
    else:
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
            ax.text(x+0.33, y+0.31, "+", ha="center", va="center", fontsize=11,
                    color=C_POS, fontweight="bold", zorder=4)
        elif q < 0:
            ax.text(x+0.33, y+0.31, "−", ha="center", va="center", fontsize=11,
                    color=C_NEG, fontweight="bold", zorder=4)
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(0.1, 6.7)
    ax.set_aspect("equal"); ax.axis("off")

# ---- figure ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "DejaVu Sans"
fig = plt.figure(figsize=(13.2, 10.6))
fig.patch.set_facecolor("white")
LX = 0.055            # common left margin
RX = 0.945            # common right margin

def rule(y, x0=LX, x1=RX, c="#e4e4df", lw=1.1):
    fig.add_artist(Line2D([x0, x1], [y, y], color=c, lw=lw, transform=fig.transFigure))

def section(y, label, note):
    fig.text(LX, y, label, fontsize=12.5, fontweight="bold", color=INK)
    fig.text(LX + 0.115, y, note, fontsize=10.5, color=MUTE, va="baseline")
    rule(y - 0.013)

# title block
fig.text(LX, 0.968, "Membranes: lipids, coarse-graining and charge",
         fontsize=18, fontweight="bold", color=INK)
fig.text(LX, 0.944,
         "Martini 3, non-polarisable water. The three lipids share tails, glycerol and phosphate; "
         "only the first head bead differs, and that is what sets the charge.",
         fontsize=10.5, color=MUTE)

# ---- section A: lipids -----------------------------------------------------------------
section(0.905, "Lipids", "all-atom (left) and Martini 3 coarse-grained (right), in shared colours")

lipids = ["POPE", "POPG", "POPC"]
colw, cgap = 0.287, 0.018
y_panel, h_panel, w_aa, inner = 0.470, 0.370, 0.148, 0.006
for k, lip in enumerate(lipids):
    x0 = LX + k * (colw + cgap)
    xc = x0 + colw / 2.0
    net = NETCHARGE[lip]
    fig.text(xc, 0.864, lip, ha="center", fontsize=14.5, fontweight="bold", color=INK)
    fig.text(xc, 0.846, ("anionic  ·  net −1" if net < 0 else "zwitterionic  ·  net 0"),
             ha="center", fontsize=9.5, color=(C_NEG if net < 0 else "#5a6472"), fontweight="bold")

    ax_aa = fig.add_axes([x0, y_panel, w_aa, h_panel])
    draw_aa(ax_aa, lip)
    x_cg = x0 + w_aa + inner
    w_cg = colw - w_aa - inner
    ax_cg = fig.add_axes([x_cg, y_panel, w_cg, h_panel])
    draw_cg(ax_cg, lip)

    hb, hq = HEADBEAD[lip]
    fig.text(xc, 0.450,
             "%s  →  %s (%s),  PO4 (−1)" %
             (HEADNAME[lip], hb, ("+1" if hq > 0 else "0")),
             ha="center", fontsize=9.5, color="#555")

# colour key (single aligned row)
rule(0.420)
fig.text(LX, 0.400, "Colours", fontsize=10.5, fontweight="bold", color=INK)
items = [("cationic head +1", C_POS), ("neutral head 0", C_NEUT), ("phosphate −1", C_NEG),
         ("glycerol / ester", C_GLY), ("saturated tail", C_SAT), ("double-bond tail", C_UNSAT)]
kx = LX + 0.085
for i, (lab, col) in enumerate(items):
    xx = kx + i * 0.147
    fig.patches.append(Circle((xx, 0.4005), 0.0066, transform=fig.transFigure,
                             facecolor=col, edgecolor="#3a3a3a", lw=0.6))
    fig.text(xx + 0.013, 0.400, lab, va="center", fontsize=9, color="#444")

# ---- section B: membranes --------------------------------------------------------------
section(0.360, "Membranes", "composition and charge balance")

def membrane_panel(ax, title, subtitle, lines, charge_lines, accent):
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.012, 0.012), 0.976, 0.976,
                 boxstyle="round,pad=0.008,rounding_size=0.025",
                 facecolor="#fbfbfa", edgecolor="#dedeea", lw=1.1, transform=ax.transAxes))
    ax.add_patch(FancyBboxPatch((0.012, 0.012), 0.012, 0.976,
                 boxstyle="square,pad=0", facecolor=accent, edgecolor="none",
                 transform=ax.transAxes))
    ax.text(0.055, 0.90, title, fontsize=12.5, fontweight="bold", transform=ax.transAxes, color=INK)
    ax.text(0.055, 0.815, subtitle, fontsize=9, transform=ax.transAxes, color=MUTE)
    y = 0.685
    for lab, val in lines:
        ax.text(0.075, y, lab, fontsize=10, transform=ax.transAxes, color="#444")
        ax.text(0.925, y, val, fontsize=10, transform=ax.transAxes, color="#222",
                ha="right", fontfamily="DejaVu Sans Mono")
        y -= 0.105
    y -= 0.03
    for cl in charge_lines:
        ax.text(0.055, y, cl, fontsize=9.5, transform=ax.transAxes, color=accent, fontweight="bold")
        y -= 0.085

pw, pgap = (RX - LX - 0.03) / 2.0, 0.03
axb = fig.add_axes([LX, 0.050, pw, 0.295])
membrane_panel(
    axb, "Bacterial target · POPE:POPG 3:1",
    "Gram-negative inner-membrane mimic. 198 lipids, ~5285 beads, 0.15 M NaCl.",
    [("POPE (zwitterionic)", "140"), ("POPG (anionic, −1 each)", "58"),
     ("W  (Martini water)", "2847"), ("Na⁺  /  Cl⁻", "60 / 2")],
    ["Membrane charge −58 e  (58 × POPG).",
     "Excess Na⁺ neutralises it → system neutral.",
     "Anionic surface."], C_NEG)

axh = fig.add_axes([LX + pw + pgap, 0.050, pw, 0.295])
membrane_panel(
    axh, "Host · POPC (neutral)",
    "Mammalian-membrane mimic. Cholesterol and sphingomyelin deferred.",
    [("POPC (zwitterionic)", "194"), ("W  (Martini water)", "~2716"),
     ("Na⁺  /  Cl⁻", "equal (salt only)")],
    ["Membrane charge 0  (zwitterionic POPC).",
     "Na⁺ = Cl⁻ : the 0.15 M salt only.",
     "Neutral surface."], "#2E6FAE")

fig.text(0.5, 0.022,
         "Selectivity contrast: the anionic bacterial surface (POPG) versus the neutral host (POPC) "
         "is what a cationic antimicrobial peptide discriminates on.",
         ha="center", fontsize=9.5, color=MUTE, style="italic")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "membranes.png")
fig.savefig(out, dpi=200, facecolor="white")
print("wrote", out)
