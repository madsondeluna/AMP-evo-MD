#!/usr/bin/env bash
#
# Open a coarse-grained membrane in VMD with Martini-friendly representations, optionally with a
# peptide and a trajectory. Lipids are drawn as VDW beads coloured by type (sticks are optional, via
# cg_bonds; see below), the peptide highlighted, water as a transparent surface, and ions small, all
# in a chalky matte style (VMD's AOChalky material with ambient occlusion). It is meant to make this
# pipeline's outputs easy to inspect, in the spirit of the Martini tutorials' visualisation setup,
# with this repository's file names.
#
# The structure is always a coordinate file VMD can read (.gro); VMD cannot open a .tpr directly.
# To draw Martini sticks, the bonds are injected from the topology by cg_bonds-v5.tcl (see below).
#
# Usage:
#   tools/vmd_view.sh prep/bacterial/equil.gro                                   # bonds from equil.tpr beside it
#   tools/vmd_view.sh prep/bacterial/equil.gro prep/bacterial/equil.xtc          # ...with the trajectory (movie)
#   tools/vmd_view.sh user/base/insert_peptide/membrane_bacterial.gro            # committed reference (.top beside it)
#   tools/vmd_view.sh output/<SEQ>/production_bacterial.gro \
#                     output/<SEQ>/production_bacterial.xtc                      # a candidate: peptide in a membrane
#
# Notes:
#   - Martini bonds. VMD draws no bonds for a coarse-grained .gro, so the lipids look like beads
#     until cg_bonds-v5.tcl adds the bonds from the topology. Make it available with
#     CG_BONDS=/path/to/cg_bonds-v5.tcl, or drop it in tools/vmd/. The bonds are read from a .tpr
#     (via gmx dump, most reliable) or a .top, auto-detected beside the structure, or VIEW_BONDS=...
#     Without it the lipids fall back to beads, which is still perfectly readable for a membrane.
#   - To watch the trajectory (the "movie"), pass the .gro and its .xtc, e.g. equil.gro equil.xtc.
#     The bilayer is kept centred over the run and every representation is smoothed (5-frame window).
#   - One representation is created per lipid type present, plus peptide, ions and water, each
#     toggleable in Graphics > Representations.
#   - Lipids are drawn as VDW beads by default (no external script needed). Sticks need Martini
#     bonds, which are optional and off unless cg_bonds is set up (see above). Force any style with
#     VIEW_LIPID_REP, e.g. VIEW_LIPID_REP="Bonds 1.0 12.0".
#   - Water is shown as a transparent surface by default; set VIEW_HIDE_WATER=1 to hide it.
#   - The periodic box outline is not drawn by default; set VIEW_SHOW_BOX=1 to draw it.
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

# ---------------------------------------------------------------------------
# Martini bonds. VMD cannot recover coarse-grained bonds from a .gro (beads sit too far apart for
# its distance search) and does not read a .tpr, so a bare membrane shows only as beads. The Martini
# cg_bonds-v5.tcl script injects the bonds from the topology; point CG_BONDS at it, or drop it in
# tools/vmd/, and the lipids are then drawn as sticks. The bonds are taken from a .tpr (via gmx dump,
# the robust route) or from a .top, auto-detected beside the structure; override with VIEW_BONDS.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CG_BONDS="${CG_BONDS:-}"
if [ -z "${CG_BONDS}" ] && [ -f "${SCRIPT_DIR}/vmd/cg_bonds-v5.tcl" ]; then
    CG_BONDS="${SCRIPT_DIR}/vmd/cg_bonds-v5.tcl"
fi

VIEW_BONDS="${VIEW_BONDS:-}"
if [ -z "${VIEW_BONDS}" ]; then
    _base="${STRUCT%.*}"
    if   [ -f "${_base}.tpr" ];                        then VIEW_BONDS="${_base}.tpr"
    elif [ -f "$(dirname "${STRUCT}")/equil.tpr" ];    then VIEW_BONDS="$(dirname "${STRUCT}")/equil.tpr"
    elif [ -f "${_base}.top" ];                        then VIEW_BONDS="${_base}.top"
    fi
fi

export VIEW_CG_BONDS VIEW_BONDS_ABS VIEW_BONDS_KIND VIEW_GMX
VIEW_CG_BONDS=""
[ -z "${CG_BONDS}" ] || [ ! -f "${CG_BONDS}" ] || VIEW_CG_BONDS="$(cd "$(dirname "${CG_BONDS}")" && pwd)/$(basename "${CG_BONDS}")"
VIEW_BONDS_ABS=""
VIEW_BONDS_KIND=""
if [ -n "${VIEW_BONDS}" ] && [ -f "${VIEW_BONDS}" ]; then
    VIEW_BONDS_ABS="$(cd "$(dirname "${VIEW_BONDS}")" && pwd)/$(basename "${VIEW_BONDS}")"
    case "${VIEW_BONDS}" in
        *.tpr) VIEW_BONDS_KIND="tpr" ;;
        *)     VIEW_BONDS_KIND="top" ;;
    esac
fi
VIEW_GMX="$(command -v gmx || true)"

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
# Add the Martini bonds with cg_bonds-v5.tcl, if it and a topology were found. Without them a
# coarse-grained .gro has no connectivity, so the lipids can only be drawn as beads (not sticks).
# cg_bonds adds the bonds to the loaded molecule; the bond count checked below then turns on the
# stick (Licorice) representation.
# ---------------------------------------------------------------------------
if {[info exists env(VIEW_CG_BONDS)] && $env(VIEW_CG_BONDS) ne "" && [file exists $env(VIEW_CG_BONDS)]
        && [info exists env(VIEW_BONDS_ABS)] && $env(VIEW_BONDS_ABS) ne ""} {
    if {[catch {source $env(VIEW_CG_BONDS)} _e]} {
        puts "cg_bonds: could not source $env(VIEW_CG_BONDS): $_e"
    } elseif {$env(VIEW_BONDS_KIND) eq "tpr"} {
        set _gmx $env(VIEW_GMX)
        if {$_gmx eq ""} { set _gmx "gmx" }
        puts "cg_bonds: adding Martini bonds from [file tail $env(VIEW_BONDS_ABS)] (gmx dump)..."
        if {[catch {cg_bonds -gmx $_gmx -tpr $env(VIEW_BONDS_ABS)} _e]} {
            puts "cg_bonds (-tpr) failed: $_e"
        }
    } else {
        puts "cg_bonds: adding Martini bonds from [file tail $env(VIEW_BONDS_ABS)] (.top)..."
        if {[catch {cg_bonds -top $env(VIEW_BONDS_ABS)} _e]} {
            puts "cg_bonds (-top) failed: $_e"
        }
    }
}

# ---------------------------------------------------------------------------
# Selections and representations. Each lipid type that is present gets its own representation, so it
# can be toggled on and off and told apart by colour; water, ions and the peptide (protein: anything
# not matched as lipid, water or ion) each get one too. Matching is by residue name, so this works
# for any composition and for a membrane with or without a peptide. Unmatched names (special
# characters like NA+ / CL-) are double-quoted, which VMD's selection parser requires.
# ---------------------------------------------------------------------------
set lipid_order {POPE POPG POPC POPS CDL0 CDL2 DPSM CHOL}
array set lipid_color {POPE 23  POPG 32  POPC 6  POPS 29  CDL0 25  CDL2 13  DPSM 19  CHOL 14}

set lipids "resname [join $lipid_order { }]"
set water  "resname W WF PW"
set ions   {resname NA CL ION "NA+" "CL-" SOD CLA}
set pept   "not ($lipids or $water or $ions)"

# ions split by sign so Na and Cl carry their own colour (names with + / - are double-quoted)
set ion_cations {resname NA "NA+" SOD ION}
set ion_anions  {resname CL "CL-" CLA}
set ion_cation_color 5    ;# Na: tan
set ion_anion_color  9    ;# Cl: pink

mol delrep 0 top

# Lipid drawing style. Override it with VIEW_LIPID_REP to try others live without editing, e.g.
#   VIEW_LIPID_REP="Bonds 1.0 12.0"      sticks without end spheres
#   VIEW_LIPID_REP="CPK 1.0 0.3 12 12"   small balls joined by thin sticks
#   VIEW_LIPID_REP="VDW 1.0 20.0"        space-filling beads (no sticks)
#   VIEW_LIPID_REP="Licorice 1.0 12 12"  the default stick look (rounded joints)
# The default is Licorice, the Martini-tutorial stick look. Licorice, Bonds and CPK all need bonds
# to draw their sticks: a .tpr carries them from the topology, but a .gro does not (VMD's distance
# search finds none for coarse-grained beads), so with no explicit override and no bonds the lipids
# fall back to small VDW beads, which is all a bare .gro can show.
set _all [atomselect top all]
set _nbonds 0
foreach _bl [$_all getbonds] { incr _nbonds [llength $_bl] }
$_all delete
if {[info exists env(VIEW_LIPID_REP)] && $env(VIEW_LIPID_REP) ne ""} {
    set lipid_rep $env(VIEW_LIPID_REP)
} else {
    set lipid_rep [expr {$_nbonds > 0 ? "Licorice 1.0 12.0 12.0" : "VDW 1.0 20.0"}]
}

# one representation per lipid type actually present, each in its own flat colour (unchanged by the
# stick/bead choice) so the types are visually separate and individually toggleable
foreach lip $lipid_order {
    set sel [atomselect top "resname $lip"]
    set n [$sel num]
    $sel delete
    if {$n == 0} { continue }
    mol representation $lipid_rep
    mol color ColorID $lipid_color($lip)
    mol selection "resname $lip"
    mol material AOChalky
    mol addrep top
}

# peptide (protein): highlighted in yellow, slightly larger, so it stands out against the bilayer
mol representation VDW 1.1 20.0
mol color ColorID 4
mol selection $pept
mol material AOChalky
mol addrep top

# ions: small beads, Na and Cl each in their own colour
mol representation VDW 0.6 16.0
mol color ColorID $ion_cation_color
mol selection $ion_cations
mol material AOChalky
mol addrep top

mol representation VDW 0.6 16.0
mol color ColorID $ion_anion_color
mol selection $ion_anions
mol material AOChalky
mol addrep top

# water: a transparent QuickSurf slab, the Martini-tutorial look. QuickSurf wraps the water beads
# in a smooth isosurface so the solvent reads as a faint slab around the bilayer rather than
# thousands of beads; a transparent material and a light blue (iceblue) keep it from hiding the
# membrane. For a glassier finish swap Transparent for Glass3 or Glass1; for a smoother or coarser
# surface change the QuickSurf grid spacing, the third number. Set VIEW_HIDE_WATER=1 to drop it.
if {![info exists env(VIEW_HIDE_WATER)] || $env(VIEW_HIDE_WATER) eq ""} {
    mol representation QuickSurf 1.0 0.5 1.0 1.0
    mol color ColorID 15
    mol selection $water
    mol material Transparent
    mol addrep top
}

# Smooth every representation's motion over a 5-frame window, so the trajectory plays without
# per-frame bead jitter. This is applied to all reps created above and is harmless on a single
# frame (a bare .gro), where there is nothing to smooth.
set nreps [molinfo top get numreps]
for {set i 0} {$i < $nreps} {incr i} {
    mol smoothrep top $i 5
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

# pbctools ships with VMD. Wrap molecules into the box so the bilayer is continuous across the
# periodic boundary (useful mainly when playing a trajectory). The box outline is NOT drawn by
# default; set VIEW_SHOW_BOX=1 to draw it. Both calls are wrapped in catch so a VMD build without
# pbctools still opens the system.
if {![catch {package require pbctools}]} {
    if {$traj ne ""} {
        # keep the bilayer centred across the whole trajectory by centring the box on the lipids'
        # centre of mass before wrapping, so the membrane neither drifts nor jumps across the
        # periodic boundary as the movie plays
        catch {pbc wrap -centersel $lipids -center com -compound res -all}
    }
    if {[info exists env(VIEW_SHOW_BOX)] && $env(VIEW_SHOW_BOX) ne ""} {
        catch {pbc box -center bb -color gray -width 1}
    }
}

display resetview
TCL

"${VMD_BIN}" -e "${TCL}"
