#!/usr/bin/env bash
#
# Open a coarse-grained membrane in VMD with Martini-friendly representations, optionally with a
# peptide and a trajectory. Lipids are coloured by bead region (tails, glycerol, phosphate and the
# head, the head split by charge), so the bacterial and host membranes look alike except at the
# head; the peptide, ions and water each have their own colour, in a muted palette and a chalky
# matte style (VMD's AOChalky material with ambient occlusion). It is meant to make this pipeline's
# outputs easy to inspect, in the spirit of the Martini tutorials' visualisation setup.
#
# The structure is always a coordinate file VMD can read (.gro); VMD cannot open a .tpr directly.
# Martini sticks come from cg_bonds, which is bundled in tools/vmd and injects the bonds (see below).
#
# Usage:
#   tools/vmd_view.sh prep/bacterial/equil.gro                                   # bonds from equil.tpr beside it
#   tools/vmd_view.sh prep/bacterial/equil.gro prep/bacterial/equil.xtc          # ...with the trajectory (movie)
#   tools/vmd_view.sh user/base/insert_peptide/membrane_bacterial.gro            # committed reference (.top beside it)
#   tools/vmd_view.sh output/<SEQ>/production_bacterial.gro \
#                     output/<SEQ>/production_bacterial.xtc                      # a candidate: peptide in a membrane
#
# Notes:
#   - Martini bonds / sticks. VMD draws no bonds for a coarse-grained .gro. cg_bonds (bundled in
#     tools/vmd/) injects them from the topology, auto-detected beside the structure: a .tpr (via
#     gmx dump) or a .top. With bonds the lipids are Licorice sticks; without them they fall back to
#     VDW beads. Point CG_BONDS=/path or VIEW_BONDS=/path to override, or VIEW_LIPID_REP to force a
#     style, e.g. VIEW_LIPID_REP="VDW 1.0 20.0".
#   - To watch the trajectory (the "movie"), pass the .gro and its .xtc, e.g. equil.gro equil.xtc.
#     The bilayer is kept centred over the run and every representation is smoothed (5-frame window).
#   - Colours use a muted palette (seaborn "deep", matching docs/membranes.png) and are tunable live
#     via env vars, each a VMD ColorID (0-32): VIEW_TAIL_COLOR, VIEW_HEADCAT_COLOR, VIEW_HEADNEU_COLOR,
#     VIEW_PO4_COLOR, VIEW_GLYC_COLOR, VIEW_PEPTIDE_COLOR, VIEW_NA_COLOR, VIEW_CL_COLOR, VIEW_WATER_COLOR.
#   - VIEW_COLOR_BY selects the lipid body scheme: "region" (default) draws every lipid body grey
#     with the head/phosphate/glycerol marks on top; "charge" additionally tints the neutral (GL0)
#     lipid, POPG, purple, leaving every other lipid grey. Purple is then the one saturated lipid
#     colour, so POPG stands out and the bacterial (POPE:POPG) and host (POPC) membranes are
#     distinguishable side by side, without a blue body competing with the blue water. The region
#     marks still overlay in both modes.
#   - Every component has its own representation, toggleable in Graphics > Representations.
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

# VMD reads coordinate files (.gro, .pdb, ...), not GROMACS run inputs. A .tpr/.top passed as the
# structure makes VMD load no topology and then segfault once a selection or representation runs, so
# swap in the sibling .gro of the same basename if there is one (the .tpr is still picked up as the
# bond source below), or stop with a clear message.
case "${STRUCT}" in
    *.tpr|*.top|*.itp|*.ndx|*.mdp|*.edr|*.log|*.cpt)
        _alt="${STRUCT%.*}.gro"
        if [ -f "${_alt}" ]; then
            echo "note: VMD cannot open .${STRUCT##*.} as a structure; using ${_alt} instead." >&2
            STRUCT="${_alt}"
        else
            echo "error: VMD cannot open a .${STRUCT##*.} as a structure. Pass the coordinate file:" >&2
            echo "       tools/vmd_view.sh ${STRUCT%.*}.gro ${TRAJ:+$TRAJ}" >&2
            exit 1
        fi
        ;;
esac

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
if [ -z "${CG_BONDS}" ]; then
    for _cand in "${SCRIPT_DIR}/vmd/cg_bonds.tcl" "${SCRIPT_DIR}/vmd/cg_bonds-v5.tcl"; do
        [ -f "${_cand}" ] && { CG_BONDS="${_cand}"; break; }
    done
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
        if {[catch {cg_bonds -gmx "$_gmx dump" -tpr $env(VIEW_BONDS_ABS)} _e]} {
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
# Selections and representations. Lipids are coloured by bead REGION (head, phosphate, glycerol,
# tails), not by lipid type, so the bacterial (POPE:POPG) and host (POPC) membranes render
# identically except at the head bead, where the chemistry and the charge differ. Water, ions and
# the peptide (anything not matched as lipid, water or ion) each get their own representation.
# Matching is by residue or bead name; names with special characters (NA+, CL-) are double-quoted,
# which VMD's selection parser requires.
# ---------------------------------------------------------------------------
# Colours are VMD ColorIDs (0-32), each overridable live from the environment, so the scheme is
# tunable without editing, e.g.:
#   VIEW_PO4_COLOR=1 VIEW_HEADCAT_COLOR=0 VIEW_NA_COLOR=10 tools/vmd_view.sh membrane.gro
proc env_or {name default} {
    if {[info exists ::env($name)] && $::env($name) ne ""} { return $::env($name) }
    return $default
}
# Colour IDs, defaulting to a harmonious muted palette (seaborn "deep", matching docs/membranes.png).
set c_tail      [env_or VIEW_TAIL_COLOR     2]   ;# lipid tails / body      : grey
set c_head_cat  [env_or VIEW_HEADCAT_COLOR 23]   ;# cationic head NC3/NH3    : blue
set c_head_neu  [env_or VIEW_HEADNEU_COLOR 11]   ;# neutral head GL0 (POPG)  : purple
set c_po4       [env_or VIEW_PO4_COLOR      1]   ;# phosphate PO4            : red
set c_glyc      [env_or VIEW_GLYC_COLOR     7]   ;# glycerol GL1/GL2         : green
set pept_color  [env_or VIEW_PEPTIDE_COLOR  3]   ;# peptide                  : orange
set water_color [env_or VIEW_WATER_COLOR   15]   ;# water surface            : light blue (transparent)

# Recolour each ColorID to the muted palette, but only when the user has not pinned it via its env
# var (an explicit override is taken as a raw VMD ColorID and left untouched).
proc harmonise {colorid envname r g b} {
    if {![info exists ::env($envname)] || $::env($envname) eq ""} {
        catch {color change rgb $colorid $r $g $b}
    }
}

set lipids "resname POPE POPG POPC POPS CDL0 CDL2 DPSM CHOL"
set water  "resname W WF PW"
set ions   {resname NA CL ION "NA+" "CL-" SOD CLA}
set pept   "not ($lipids or $water or $ions)"

# ions split by sign so Na and Cl carry their own colour. Defaults are far from the lipid palette
# and from each other (the earlier red3/orange2 pair read as a single red).
set ion_cations {resname NA "NA+" SOD ION}
set ion_anions  {resname CL "CL-" CLA}
set ion_cation_color [env_or VIEW_NA_COLOR 10]   ;# Na : cyan
set ion_anion_color  [env_or VIEW_CL_COLOR  9]   ;# Cl : pink

# Apply the harmonious palette (seaborn "deep") to the default colour IDs.
harmonise $c_tail           VIEW_TAIL_COLOR    0.72 0.72 0.72
harmonise $c_head_cat       VIEW_HEADCAT_COLOR 0.30 0.45 0.69
harmonise $c_head_neu       VIEW_HEADNEU_COLOR 0.51 0.45 0.70
harmonise $c_po4            VIEW_PO4_COLOR     0.77 0.31 0.32
harmonise $c_glyc           VIEW_GLYC_COLOR    0.33 0.66 0.41
harmonise $pept_color       VIEW_PEPTIDE_COLOR 0.87 0.52 0.32
harmonise $ion_cation_color VIEW_NA_COLOR      0.39 0.71 0.80
harmonise $ion_anion_color  VIEW_CL_COLOR      0.85 0.55 0.76
harmonise $water_color      VIEW_WATER_COLOR   0.65 0.81 0.89

mol delrep 0 top

# Lipid body drawing style. Licorice (sticks) needs the Martini bonds, which cg_bonds adds from the
# topology when available (it ships in tools/vmd/); without bonds this falls back to VDW beads.
# Override with VIEW_LIPID_REP, e.g. VIEW_LIPID_REP="VDW 1.0 20.0" or VIEW_LIPID_REP="Bonds 1.0 12".
set _all [atomselect top all]
set _nbonds 0
foreach _bl [$_all getbonds] { incr _nbonds [llength $_bl] }
$_all delete
if {[info exists env(VIEW_LIPID_REP)] && $env(VIEW_LIPID_REP) ne ""} {
    set lipid_rep $env(VIEW_LIPID_REP)
} else {
    set lipid_rep [expr {$_nbonds > 0 ? "Licorice 0.5 12.0 12.0" : "VDW 0.8 20.0"}]
}

# Lipid body. Two schemes, chosen by VIEW_COLOR_BY:
#   region (default) one grey chain, so the Licorice bonds draw continuously and the coloured
#                    region marks below carry all the colour.
#   charge           the neutral (GL0) lipid, POPG, is tinted purple and every other lipid stays
#                    grey. Whole lipids stay in one representation, so the bonds remain continuous,
#                    and the region marks still overlay. Purple is the one saturated lipid colour,
#                    so POPG pops and the bacterial (POPE:POPG) and host (POPC) membranes are told
#                    apart side by side, with no blue body competing with the blue water.
set color_by [env_or VIEW_COLOR_BY region]
if {$color_by eq "charge"} {
    set neu_body  "($lipids) and (same residue as (name GL0))"
    set rest_body "($lipids) and not (same residue as (name GL0))"
    foreach pair [list [list $neu_body $c_head_neu] [list $rest_body $c_tail]] {
        set sel [lindex $pair 0]
        set probe [atomselect top $sel]
        set n [$probe num]
        $probe delete
        if {$n == 0} { continue }
        mol representation $lipid_rep
        mol color ColorID [lindex $pair 1]
        mol selection $sel
        mol material AOChalky
        mol addrep top
    }
} else {
    mol representation $lipid_rep
    mol color ColorID $c_tail
    mol selection $lipids
    mol material AOChalky
    mol addrep top
}

# Region highlights overlaid as coloured spheres on the head, phosphate and glycerol beads. Each is
# added only if present, so this adapts to any composition and stays identical across membranes.
set lipid_marks [list \
    [list "name NC3 NH3" $c_head_cat 1.1] \
    [list "name GL0"     $c_head_neu 1.1] \
    [list "name PO4"     $c_po4      1.0] \
    [list "name GL1 GL2" $c_glyc     0.9]]
foreach mark $lipid_marks {
    set sel "($lipids) and ([lindex $mark 0])"
    set probe [atomselect top $sel]
    set n [$probe num]
    $probe delete
    if {$n == 0} { continue }
    mol representation VDW [lindex $mark 2] 16.0
    mol color ColorID [lindex $mark 1]
    mol selection $sel
    mol material AOChalky
    mol addrep top
}

# peptide (protein): highlighted, slightly larger, so it stands out against the bilayer
mol representation VDW 1.1 20.0
mol color ColorID $pept_color
mol selection $pept
mol material AOChalky
mol addrep top

# ions: small beads, Na and Cl each in their own colour
mol representation VDW 0.8 16.0
mol color ColorID $ion_cation_color
mol selection $ion_cations
mol material AOChalky
mol addrep top

mol representation VDW 0.8 16.0
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
    mol color ColorID $water_color
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
