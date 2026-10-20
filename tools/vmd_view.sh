#!/usr/bin/env bash
#
# Open a coarse-grained membrane in VMD with Martini-friendly representations, optionally with a
# peptide and a trajectory. Lipids are drawn as beads coloured by type, the peptide is
# highlighted, water as a transparent surface, ions small, and the periodic box drawn, all in a chalky
# matte style (VMD's AOChalky material with ambient occlusion). It is meant to make this
# pipeline's outputs easy to inspect, in the spirit of the Martini tutorials' visualisation setup,
# with this repository's file names.
#
# Usage:
#   tools/vmd_view.sh user/base/insert_peptide/membrane_bacterial.gro            # the committed reference
#   tools/vmd_view.sh prep/bacterial/equil.gro                                   # a membrane you just built
#   tools/vmd_view.sh prep/bacterial/min.tpr prep/bacterial/equil.xtc            # with a trajectory
#   tools/vmd_view.sh output/<SEQ>/production_bacterial.tpr \
#                     output/<SEQ>/production_bacterial.xtc                      # a peptide in a membrane
#
# The first line runs straight after a clone or pull: the equilibrated bacterial reference membrane
# ships with the repository, so there is nothing to build first.
#
# Notes:
#   - A .tpr as the structure carries the bonds, so the beads are drawn connected. A .gro has no
#     connectivity, so beads show without sticks, which is still perfectly readable for a membrane.
#   - Water is hidden by default. Edit the water representation block below to show it.
#   - Set VMD=/path/to/vmd if vmd is not on PATH (for example a module's binary).
#
set -euo pipefail

STRUCT="${1:?usage: vmd_view.sh <structure .gro|.tpr> [trajectory .xtc|.trr]}"
TRAJ="${2:-}"
VMD_BIN="${VMD:-vmd}"

command -v "${VMD_BIN}" >/dev/null 2>&1 || {
    echo "vmd not found on PATH. Set VMD=/path/to/vmd, or load the VMD module." >&2
    exit 1
}
[ -f "${STRUCT}" ] || { echo "structure not found: ${STRUCT}" >&2; exit 1; }
[ -z "${TRAJ}" ] || [ -f "${TRAJ}" ] || { echo "trajectory not found: ${TRAJ}" >&2; exit 1; }

# The structure and trajectory are passed through the environment so the TCL can load them with
# mol new / mol addfile, which puts the trajectory into the same molecule as the structure. That
# is more reliable than listing both on the vmd command line, where the trajectory can land in a
# second molecule instead.
export VIEW_STRUCT VIEW_TRAJ
VIEW_STRUCT="$(cd "$(dirname "${STRUCT}")" && pwd)/$(basename "${STRUCT}")"
VIEW_TRAJ=""
[ -z "${TRAJ}" ] || VIEW_TRAJ="$(cd "$(dirname "${TRAJ}")" && pwd)/$(basename "${TRAJ}")"

TCL="$(mktemp "${TMPDIR:-/tmp}/vmd_cg.XXXXXX.tcl")"
trap 'rm -f "${TCL}"' EXIT

cat > "${TCL}" <<'TCL'
# ---------------------------------------------------------------------------
# Load structure, then trajectory (if any) into the same molecule.
# ---------------------------------------------------------------------------
set struct $env(VIEW_STRUCT)
set traj   $env(VIEW_TRAJ)

mol new $struct waitfor all
if {$traj ne ""} {
    mol addfile $traj waitfor all
}

# Coarse-grained beads carry no van der Waals radius, so VDW spheres would be tiny. A uniform
# ~2 A radius makes a Martini bilayer read as space-filling without being faithful to the small
# differences between bead classes, which do not matter for inspection.
set everything [atomselect top all]
$everything set radius 2.0
$everything delete

# ---------------------------------------------------------------------------
# Selections. Lipids, water and ions are matched by residue name; the peptide is whatever is
# left, so this works for any composition and for a membrane with or without a peptide.
# ---------------------------------------------------------------------------
set lipids "resname POPE POPG POPC POPS CDL0 CDL2 DPSM CHOL"
set water  "resname W WF PW"
set ions   "resname NA CL ION NA+ CL- SOD CLA"
set pept   "not ($lipids or $water or $ions)"

mol delrep 0 top

# lipids: beads coloured by residue name, so POPE and POPG are distinguishable
mol representation VDW 1.0 20.0
mol color ResName
mol selection "$lipids"
mol material AOChalky
mol addrep top

# peptide: highlighted in yellow, slightly larger, so it stands out against the bilayer
mol representation VDW 1.1 20.0
mol color ColorID 4
mol selection "$pept"
mol material AOChalky
mol addrep top

# ions: small beads coloured by name
mol representation VDW 0.6 16.0
mol color Name
mol selection "$ions"
mol material AOChalky
mol addrep top

# water: a transparent surface, the Martini-tutorial look. QuickSurf wraps the water beads in a
# smooth isosurface so the solvent reads as a faint slab around the bilayer rather than thousands
# of beads; a transparent material and a light blue (iceblue) keep it from hiding the membrane.
# For a glassier finish swap Transparent for Glass3 or Glass1; for a smoother or coarser surface
# change the QuickSurf grid spacing, the third number. Set VIEW_HIDE_WATER=1 to drop it entirely.
if {![info exists env(VIEW_HIDE_WATER)] || $env(VIEW_HIDE_WATER) eq ""} {
    mol representation QuickSurf 1.0 0.5 1.0 1.0
    mol color ColorID 15
    mol selection "$water"
    mol material Transparent
    mol addrep top
}

# ---------------------------------------------------------------------------
# Display and periodic box.
# ---------------------------------------------------------------------------
display projection Orthographic
display depthcue off
color Display Background white
axes location Off

# Chalky, matte look. The AOChalky material on the representations above gives the clay/chalk
# surface; ambient occlusion is the soft self-shadowing in the creases that makes it read as
# matte rather than plastic, and a directional shadow adds depth. These render live in the
# interactive window on a recent VMD, and fully in a ray-traced snapshot. For a clean still,
# use File > Render with the Tachyon or TachyonLOptiXInternal renderer, or from the console:
#   render TachyonInternal view.tga
display shadows on
display ambientocclusion on
display aoambient 0.9
display aodirect 0.3
display antialias on
display rendermode GLSL

# pbctools ships with VMD. Wrap molecules into the box so the bilayer is continuous, then draw
# the box. Both are wrapped in catch so a VMD build without pbctools still opens the system.
if {![catch {package require pbctools}]} {
    if {$traj ne ""} {
        catch {pbc wrap -center bb -compound res -all}
    }
    catch {pbc box -center bb -color gray -width 1}
}

display resetview
TCL

"${VMD_BIN}" -e "${TCL}"
