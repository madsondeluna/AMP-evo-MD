#!/usr/bin/env bash
#
# Membrane preparation protocol, end to end.
#
#   tools/protocol.sh bacterial "POPE:70,POPG:20,CDL2:10" target
#   tools/protocol.sh host      "POPC:60,DPSM:10,CHOL:30" host
#
# Runs: build -> minimise -> equilibrate at zero tension -> re-plan the lipid count ->
# rebuild -> production at constant area -> extract -> analyse -> emit the MembraneDef.
#
# The re-plan step is why this is a script and not a single command. The constant-area box must
# sit at the composition's equilibrium area, and the equilibrium area per lipid is only known
# after equilibrating, so the lipid count is corrected on a second pass. The box never changes;
# the lipid count does.
#
# Requires gmx and insane on PATH. Set GMX to use a different binary.

set -euo pipefail

NAME="${1:?usage: protocol.sh <name> <composition> <role>}"
COMPOSITION="${2:?composition, e.g. POPE:70,POPG:20,CDL2:10}"
ROLE="${3:-target}"

GMX="${GMX:-gmx}"
# Membrane preparation is a handful of long runs rather than many short ones, so here a single
# process with many threads is the right shape, unlike the candidate evaluations. Leaves one
# core free by default; override with NT.
NT="${NT:-$(python3 -c "
import os
try:
    n = len(os.sched_getaffinity(0))
except Exception:
    n = os.cpu_count() or 1
print(max(1, n - 1))
")}"
BOX="${BOX:-8.0}"              # lateral box in nm, both dimensions
MAX_RESIDUES="${MAX_RESIDUES:-30}"
SALT="${SALT:-0.15}"           # mol/L, on top of neutralising counterions
EQUIL_NS="${EQUIL_NS:-200}"    # semi-isotropic, to find the equilibrium area
PROD_NS="${PROD_NS:-400}"      # constant area, the reference run
ROOT="${ROOT:-prep}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Absolute before the cd. ROOT is relative by default, so after 'cd "${DIR}"' a later
# 'cd "${HERE}"' would leave the analysis commands resolving ${DIR} against the repository root
# rather than against the directory the script was invoked from.
mkdir -p "${ROOT}/${NAME}"
DIR="$(cd "${ROOT}/${NAME}" && pwd)"
cd "${DIR}"

MDP_DIR="${HERE}/user/base/production"

say() { printf '\n=== %s ===\n' "$1"; }

# ---------------------------------------------------------------------------
# Everything this script needs, checked before the first step rather than discovered at the one
# that happens to use it. The build is step 2 and the analysis is step 7, so a missing insane or
# a missing matplotlib used to surface after minutes or after the whole equilibration.
missing=""
command -v "${GMX}" >/dev/null 2>&1 || missing="${missing}\n  ${GMX} (GROMACS): not on PATH. Load the module, or set GMX=<path>."
command -v insane  >/dev/null 2>&1 || missing="${missing}\n  insane: not on PATH. It builds the membrane; see 'insane, for building membranes' in the README."
command -v python3 >/dev/null 2>&1 || missing="${missing}\n  python3: not on PATH."

if command -v python3 >/dev/null 2>&1; then
    python3 - <<'CHECKPY' || missing="${missing}\n  python packages: see the message above."
import sys
absent = []
for module in ("numpy", "matplotlib"):
    try:
        __import__(module)
    except ImportError:
        absent.append(module)
if absent:
    sys.stderr.write("  missing python package(s): %s. Activate the virtualenv "
                     "(source ~/amp-env.sh) or pip install them.\n" % ", ".join(absent))
    sys.exit(1)
CHECKPY
fi

if [ -n "${missing}" ]; then
    printf 'Cannot run the membrane protocol. Missing:%b\n' "${missing}" >&2
    exit 1
fi

if [ ! -d "${MDP_DIR}" ]; then
    echo "Cannot find the production .mdp directory: ${MDP_DIR}" >&2
    echo "This script resolves it relative to its own location, so it has to stay inside the" >&2
    echo "repository; run it as tools/protocol.sh rather than from a copy elsewhere." >&2
    exit 1
fi

echo "Using ${NT} thread(s) for GROMACS. Override with NT=<n>."
say "1/7  plan"
# ---------------------------------------------------------------------------
python3 "${HERE}/tools/membrane.py" plan \
    --name "${NAME}" --composition "${COMPOSITION}" \
    --box "${BOX}" --max-residues "${MAX_RESIDUES}" --salt "${SALT}" \
    | tee plan_pass1.txt

echo
if [ "${PASS:-1}" = "2" ]; then
    echo "(PASS=2: the plan above is informational. Production uses the system already in"
    echo " ${DIR}/equil.gro, which you rebuilt from plan_pass2.txt.)"
    echo
fi
echo "The plan above contains the insane command for this composition. Run it, or let this"
echo "script run the commands it printed if you have reviewed them. Stopping here on the first"
echo "pass is deliberate: lipid naming and the insane flags differ between Martini lipid sets,"
echo "and a silently wrong -u/-l split produces an asymmetric bilayer that looks fine until the"
echo "surface tension refuses to reach zero."
echo
if [ "${AUTO:-0}" != "1" ]; then
    echo "Re-run with AUTO=1 once the plan looks right, to execute the remaining steps."
    exit 0
fi

if [ "${PASS:-1}" = "2" ]; then
    # Pass 2 resumes from a system the user rebuilt and re-equilibrated by hand, so rebuilding
    # here would overwrite it with the pass-1 lipid counts and feed those into production. That
    # is the one outcome this whole two-pass procedure exists to avoid.
    if [ ! -f equil.gro ]; then
        echo "PASS=2 expects an equilibrated system in ${DIR}/equil.gro." >&2
        echo "Rebuild with the lipid counts from plan_pass2.txt, re-run the minimisation and" >&2
        echo "the zero-tension equilibration, then run PASS=2 again." >&2
        exit 1
    fi
    echo "PASS=2: resuming from the existing equil.gro, skipping build and equilibration."
    say "6/7  production at constant area (the peptide-free reference run)"
else

# The insane command is extracted from the plan rather than duplicated here, so that the lipid
# counts cannot drift between the two. -A2 assumes the plan prints it over exactly three lines,
# which is the kind of assumption that breaks silently: a truncated command still runs, and an
# insane call that lost its -u leaflet string or its -salt produces a membrane that looks
# plausible and is wrong. So the result is asserted before it is evaluated.
if ! grep -q '^insane' plan_pass1.txt; then
    echo "Could not find the insane command in ${DIR}/plan_pass1.txt." >&2
    echo "tools/membrane.py plan is expected to print it starting at the beginning of a line." >&2
    exit 1
fi

BUILD_CMD=$(grep -A2 '^insane' plan_pass1.txt | tr '\n' ' ' | sed 's/\\ / /g')

for required in " -l " " -u " " -sol " " -x " " -y "; do
    case "${BUILD_CMD}" in
        *"${required}"*) ;;
        *)
            echo "The insane command extracted from plan_pass1.txt is missing '${required# }':" >&2
            echo "  ${BUILD_CMD}" >&2
            echo "It was most likely truncated: the extraction reads the three lines starting" >&2
            echo "at 'insane', and the plan now prints more. Build by hand from plan_pass1.txt" >&2
            echo "or fix the extraction in tools/protocol.sh." >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
say "2/7  build"
# ---------------------------------------------------------------------------
eval "${BUILD_CMD}"

# insane writes a topology for the lipid set its own templates use, which is Martini 2: a
# generic '#include "martini.itp"' that does not exist here, and ions named NA+ and CL-. This
# pipeline is Martini 3, where the force field is split across named itp files and the ions are
# NA and CL. Rewrite both so grompp can read the system, right after the build, so every grompp
# that follows sees a corrected topology. Assumes ff = martini3001 (the configured force field);
# a martini22 run would need the Martini 2 include set instead.
FF_DIR="${HERE}/user/base/insert_peptide"
for itp in martini_v3.0.0.itp martini_v3.0.0_phospholipids_v1.itp \
           martini_v3.0.0_solvents_v2.itp martini_v3.0.0_ions_v1.itp; do
    if [ ! -f "${FF_DIR}/${itp}" ]; then
        echo "Missing force field file: ${FF_DIR}/${itp}" >&2
        exit 1
    fi
    cp "${FF_DIR}/${itp}" .
done
# The placeholder include becomes the Martini 3 set: the core field, the phospholipids (POPE,
# POPG and the rest), the solvents (W), and the ions (NA, CL).
sed -i 's|#include "martini.itp"|#include "martini_v3.0.0.itp"\n#include "martini_v3.0.0_phospholipids_v1.itp"\n#include "martini_v3.0.0_solvents_v2.itp"\n#include "martini_v3.0.0_ions_v1.itp"|' "${NAME}.top"
# NA+ -> NA and CL- -> CL in both the topology and the structure. The replacements are the same
# width as what they replace, so the fixed-column .gro stays aligned.
sed -i -e 's/NA+/NA /g' -e 's/CL-/CL /g' "${NAME}.top" "${NAME}_raw.gro"

# ---------------------------------------------------------------------------
say "3/7  minimise"
# ---------------------------------------------------------------------------
sed "s|\$INCLUDE_CWD|-I$(pwd)|" "${MDP_DIR}/steep_softcore_system.mdp" > min.mdp
"${GMX}" grompp -f min.mdp -c "${NAME}_raw.gro" -p "${NAME}.top" -o min.tpr -maxwarn 2
"${GMX}" mdrun -deffnm min -nt "${NT}"

# ---------------------------------------------------------------------------
say "4/7  equilibrate at zero tension (semi-isotropic) to find the equilibrium area"
# ---------------------------------------------------------------------------
# Semi-isotropic with zero reference tension lets the area relax. This is NOT the production
# ensemble: production fixes the area so that gamma can respond to the peptide.
cat > equil.mdp <<EOF
integrator               = md
dt                       = 0.02
nsteps                   = $(python3 -c "print(int(${EQUIL_NS}*1000/0.02))")
nstxout-compressed       = 50000
nstenergy                = 1000
cutoff-scheme            = Verlet
nstlist                  = 20
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale
tc-grps                  = Solvent non-Solvent
tau_t                    = 1.0 1.0
ref_t                    = 310 310
Pcoupl                   = C-rescale
Pcoupltype               = semiisotropic
tau_p                    = 12.0
compressibility          = 3e-4 3e-4
ref_p                    = 1.0 1.0
gen_vel                  = yes
gen_temp                 = 310
refcoord_scaling         = all
EOF
"${GMX}" grompp -f equil.mdp -c min.gro -p "${NAME}.top" -o equil.tpr -maxwarn 2
"${GMX}" mdrun -deffnm equil -nt "${NT}"

# ---------------------------------------------------------------------------
say "5/7  measure the equilibrium area per lipid and re-plan"
# ---------------------------------------------------------------------------
printf 'Box-X\nBox-Y\n\n' | "${GMX}" energy -f equil.edr -o equil_box.xvg -b $((EQUIL_NS/2))000
MEAN_APL=$(python3 - <<PY
import sys
sys.path.insert(0, "${HERE}/tools")
import mdtools
import membrane

t, series, _ = mdtools.read_xvg("equil_box.xvg")


def pick(prefix):
    # A bare next() over a generator raises StopIteration here, which inside a heredoc reports
    # no file, no term and no list of what was actually found. The same lookup in
    # tools/membrane_analysis.py names the labels it saw, and so does this one now.
    for key, values in series.items():
        if key.strip().startswith(prefix):
            return values
    raise SystemExit(
        "equil_box.xvg has no %s column. Found: %s.\n"
        "That file comes from 'gmx energy' on equil.edr with the box terms selected; if the "
        "term indices differ in this GROMACS version, the selection in step 5 picked the "
        "wrong ones." % (prefix, ", ".join(repr(k.strip()) for k in series) or "nothing"))


bx = pick("Box-X")
by = pick("Box-Y")
area = (bx * by).mean()

_x, _y, _z, atoms = mdtools.read_gro("equil.gro")
# Reuse membrane.count_lipids rather than matching bead names here. An inline copy is how this
# diverged before: it matched PO4 and ROH only, which misses cardiolipin entirely (its beads are
# PO41 and PO42), so the measured value driving the whole pass-2 correction was wrong by about
# 11 percent for the shipped bacterial mixture.
counts, _zs = membrane.count_lipids([(r, a, z) for r, a, _px, _py, z in atoms])
n_lipids = sum(counts.values())
if n_lipids == 0:
    sys.exit("no lipids found in equil.gro; check the head bead names")
print("%.4f" % (area / (n_lipids / 2.0)))
PY
)
echo "measured mean area per lipid: ${MEAN_APL} nm2"

python3 "${HERE}/tools/membrane.py" plan \
    --name "${NAME}" --composition "${COMPOSITION}" \
    --box "${BOX}" --max-residues "${MAX_RESIDUES}" --salt "${SALT}" \
    --apl "MEAN=${MEAN_APL}" | tee plan_pass2.txt

echo
echo "Pass 2 uses the measured area per lipid. Rebuild with its lipid counts, re-run the"
echo "minimisation and the zero-tension equilibration, and then run this script again with"
echo "PASS=2 to continue from the rebuilt equil.gro into production."
exit 0

fi   # end of the pass-1 branch

# ---------------------------------------------------------------------------
cp equil.gro "${NAME}_eq.gro"

# The production mdp is written for peptide-in-membrane systems: its tc-grps names Protein,
# Solvent and Lipids, which the pipeline builds with its own make_ndx. This membrane has no
# peptide, so a Protein group cannot exist and grompp would abort on it. Couple the two groups
# that do exist here instead, leaving every other setting in npt.mdp untouched so the reference
# run and the candidate runs share an identical ensemble, which is the whole point of measuring
# a reference.
sed -e "s|\$INCLUDE_CWD|-I$(pwd)|" \
    -e "s|\$SIM_DT|0.03|" \
    -e "s|\$SIM_NSTEPS|$(python3 -c "print(int(${PROD_NS}*1000/0.03))")|" \
    -e "s|^tc-grps.*|tc-grps                  = Solvent non-Solvent|" \
    -e "s|^tau_t.*|tau_t                    = 1.0 1.0|" \
    -e "s|^ref_t.*|ref_t                    = 310 310|" \
    -e "s|^energygrps.*|energygrps               =|" \
    "${MDP_DIR}/npt.mdp" > prod.mdp

# Solvent and non-Solvent are GROMACS default groups, so no index file is needed. Confirm the
# substitution took, because a silently unmodified mdp fails later and more confusingly.
if grep -q "Protein" prod.mdp; then
    echo "prod.mdp still references a Protein group; the tc-grps substitution did not apply." >&2
    exit 1
fi

"${GMX}" grompp -f prod.mdp -c "${NAME}_eq.gro" -p "${NAME}.top" -o "${NAME}.tpr" -maxwarn 2
"${GMX}" mdrun -deffnm "${NAME}" -nt "${NT}"

# ---------------------------------------------------------------------------
say "7/7  extract, analyse, emit the MembraneDef"
# ---------------------------------------------------------------------------
# -b uses the SAME fraction the candidate runs use, read from usersettings rather than hardcoded.
# This matters: delta_gamma is the difference between this reference and a per-candidate value, so
# averaging the two over different windows leaves any residual drift inside the observable. The
# previous value here was half the trajectory while production takes the last five percent.
PROD_START=$(python3 -c "
import sys
sys.path.insert(0, '${HERE}')
import user.usersettings as U
print(int(${PROD_NS} * 1000 * U.ENERGY_init_fraction))
")
echo "Averaging the reference over the last $(python3 -c "print(${PROD_NS} - int(${PROD_START})//1000)") ns, matching the production window."
printf 'Box-X\nBox-Y\n#Surf*SurfTen\n\n' \
    | "${GMX}" energy -f "${NAME}.edr" -o "${NAME}_energy.xvg" -b "${PROD_START}" \
      > "${NAME}_free.ENERGY" 2>&1

cd "${HERE}"
python3 tools/membrane_analysis.py single --name "${NAME}" --dir "${DIR}"
python3 tools/membrane.py calibrate --name "${NAME}" --role "${ROLE}" \
    --gro "${DIR}/${NAME}_eq.gro" --energy "${DIR}/${NAME}_free.ENERGY"

echo
echo "Paste the MembraneDef above into MEMBRANES in user/membranes.py, and copy"
echo "${DIR}/${NAME}_eq.gro and the matching .top into user/base/insert_peptide/ as"
echo "membrane_${NAME}.gro and membrane_${NAME}.top."
