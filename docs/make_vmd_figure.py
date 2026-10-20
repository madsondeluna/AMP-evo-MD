# Build docs/vmd_view.png: a usage/dataflow diagram for tools/vmd_view.sh -- what it takes, what it
# does, and what you end up looking at. Same visual language as make_membranes_figure.py.
# Regenerate with: python docs/make_vmd_figure.py
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

INK, MUTE, FAINT = "#1a1a1a", "#6b6b6b", "#9a9a9a"
C_LIP   = "#55A868"   # lipids (coloured by residue in VMD; green stands for "the bilayer")
C_PEP   = "#E1A100"   # peptide (VMD ColorID 4, yellow)
C_ION   = "#4C72B0"   # ions
C_WAT   = "#9ec9e0"   # water surface (iceblue, transparent)
PANEL_FC, PANEL_EC = "#fbfbfa", "#dedeea"

plt.rcParams["font.family"] = "DejaVu Sans"
fig = plt.figure(figsize=(12.8, 6.7))
fig.patch.set_facecolor("white")
LX, RX = 0.04, 0.96

def rule(y, x0=LX, x1=RX, c="#e4e4df", lw=1.1):
    from matplotlib.lines import Line2D
    fig.add_artist(Line2D([x0, x1], [y, y], color=c, lw=lw, transform=fig.transFigure))

fig.text(LX, 0.945, "Visualising a membrane in VMD",
         fontsize=17, fontweight="bold", color=INK)
fig.text(LX, 0.902,
         "tools/vmd_view.sh  ·  Martini-friendly representations, auto-detected from residue names. "
         "Works on a bare membrane or on a peptide in a membrane.",
         fontsize=10.5, color=MUTE)
rule(0.875)

# ---- three stage panels ---------------------------------------------------------------
py, ph = 0.300, 0.500          # panel bottom, height
pA = [0.040, py, 0.250, ph]
pB = [0.368, py, 0.280, ph]
pC = [0.726, py, 0.234, ph]

def panel(rect, tag, title, accent):
    ax = fig.add_axes(rect); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch((0.01, 0.01), 0.98, 0.98,
                 boxstyle="round,pad=0.006,rounding_size=0.03",
                 facecolor=PANEL_FC, edgecolor=PANEL_EC, lw=1.1, transform=ax.transAxes))
    ax.text(0.06, 0.885, tag, fontsize=9, fontweight="bold", color=accent, transform=ax.transAxes)
    ax.text(0.06, 0.80, title, fontsize=12.5, fontweight="bold", color=INK, transform=ax.transAxes)
    return ax

# panel A: inputs
axA = panel(pA, "1 · INPUTS", "What you give it", "#2E6FAE")
axA.text(0.06, 0.63, "structure", fontsize=10.5, fontweight="bold", color="#333", transform=axA.transAxes)
axA.text(0.06, 0.555, "equil.gro   (coordinates)", fontsize=10, color="#222",
         family="DejaVu Sans Mono", transform=axA.transAxes)
axA.text(0.06, 0.475, "VMD can't open a .tpr; the bonds\ncome from cg_bonds (see RUN)",
         fontsize=8.8, color=MUTE, transform=axA.transAxes, va="top")
axA.text(0.06, 0.30, "trajectory", fontsize=10.5, fontweight="bold", color="#333", transform=axA.transAxes)
axA.text(0.06, 0.225, "equil.xtc  (optional)", fontsize=10, color="#222",
         family="DejaVu Sans Mono", transform=axA.transAxes)
axA.text(0.06, 0.145, "the movie over time; candidates\nneed WRITE_XTC = True",
         fontsize=8.8, color=MUTE, transform=axA.transAxes, va="top")

# panel B: what it does
axB = panel(pB, "2 · RUN", "tools/vmd_view.sh", "#1a7a5e")
cmd = [("tools/vmd_view.sh \\", "#222"),
       ("  <structure> [trajectory]", "#222")]
yy = 0.655
for t, c in cmd:
    axB.text(0.055, yy, t, fontsize=9.5, color=c, family="DejaVu Sans Mono", transform=axB.transAxes)
    yy -= 0.072
notes = [
    "one rep per lipid type + peptide,",
    "ions and water (each toggleable)",
    "sticks via cg_bonds-v5.tcl (CG_BONDS)",
    "trajectory centred, smoothed (5)",
    "VIEW_HIDE_WATER=1 · VIEW_SHOW_BOX=1",
]
yy = 0.44
for n in notes:
    axB.text(0.055, yy, "•  " + n, fontsize=9, color="#444", transform=axB.transAxes)
    yy -= 0.082

# panel C: what you see
axC = panel(pC, "3 · WHAT YOU SEE", "Representations", "#9a3b8f")
rows = [
    (C_LIP, "Lipids", "one rep per type, by colour"),
    (C_PEP, "Peptide", "yellow, slightly larger"),
    (C_ION, "Ions", "small beads, by name"),
    (C_WAT, "Water", "transparent QuickSurf slab"),
]
yy = 0.64
for col, lab, desc in rows:
    axC.add_patch(Circle((0.10, yy), 0.045, facecolor=col, edgecolor="#3a3a3a", lw=0.7,
                         transform=axC.transAxes))
    axC.text(0.20, yy + 0.018, lab, fontsize=10.5, fontweight="bold", color="#222",
             transform=axC.transAxes, va="center")
    axC.text(0.20, yy - 0.055, desc, fontsize=8.6, color=MUTE, transform=axC.transAxes, va="center")
    yy -= 0.155

# ---- arrows between panels ------------------------------------------------------------
ymid = py + ph / 2.0
for x0, x1 in [(pA[0] + pA[2], pB[0]), (pB[0] + pB[2], pC[0])]:
    fig.add_artist(FancyArrowPatch((x0 + 0.004, ymid), (x1 - 0.004, ymid),
                   transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
                   lw=2.0, color="#9aa0a8"))

# ---- try-it / render callouts ---------------------------------------------------------
def callout(y, rect_x, rect_w, tag, tagcol, text, stacked=False, tagx=0.205):
    ax = fig.add_axes([rect_x, y, rect_w, 0.115]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch((0.004, 0.08), 0.992, 0.84,
                 boxstyle="round,pad=0.004,rounding_size=0.10",
                 facecolor="#f4f6f8", edgecolor="#d7dde3", lw=1.0, transform=ax.transAxes))
    if stacked:
        ax.text(0.04, 0.70, tag, fontsize=9, fontweight="bold", color=tagcol,
                transform=ax.transAxes, va="center")
        ax.text(0.04, 0.30, text, fontsize=9.6, color="#222", family="DejaVu Sans Mono",
                transform=ax.transAxes, va="center")
    else:
        ax.text(0.03, 0.55, tag, fontsize=9, fontweight="bold", color=tagcol,
                transform=ax.transAxes, va="center")
        ax.text(tagx, 0.55, text, fontsize=9.8, color="#222", family="DejaVu Sans Mono",
                transform=ax.transAxes, va="center")

callout(0.145, 0.040, 0.608, "TRY IT NOW",
        "#1a7a5e", "tools/vmd_view.sh user/base/insert_peptide/membrane_bacterial.gro", tagx=0.165)
callout(0.145, 0.672, 0.288, "STILL IMAGE  (from the VMD console)",
        "#9a3b8f", "render TachyonInternal view.tga", stacked=True)

fig.text(0.5, 0.045,
         "The reference membrane ships with the repository. VMD draws no bonds for a bare .gro, so "
         "make cg_bonds-v5.tcl available (CG_BONDS=...) to render the Martini sticks.",
         ha="center", fontsize=9.3, color=MUTE, style="italic")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vmd_view.png")
fig.savefig(out, dpi=200, facecolor="white")
print("wrote", out)
