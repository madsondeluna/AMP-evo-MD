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
echo "Using ${NT} thread(s) for GROMACS. Override with NT=<n>."
say "1/7  plan"
# ---------------------------------------------------------------------------
python3 "${HERE}/tools/membrane.py" plan \
    --name "${NAME}" --composition "${COMPOSITION}" \
    --box "${BOX}" --max-residues "${MAX_RESIDUES}" --salt "${SALT}" \
    | tee plan_pass1.txt

echo
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

BUILD_CMD=$(grep -A2 '^insane' plan_pass1.txt | tr '\n' ' ' | sed 's/\\ / /g')

# ---------------------------------------------------------------------------
say "2/7  build"
# ---------------------------------------------------------------------------
eval "${BUILD_CMD}"

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
import sys, os
sys.path.insert(0, "${HERE}/tools")
import mdtools, collections
t, s, _ = mdtools.read_xvg("equil_box.xvg")
bx = next(v for k, v in s.items() if k.strip().startswith("Box-X"))
by = next(v for k, v in s.items() if k.strip().startswith("Box-Y"))
area = (bx * by).mean()
_x, _y, _z, atoms = mdtools.read_gro("equil.gro")
heads = sum(1 for r, a, *_ in atoms if a in ("PO4", "ROH"))
print("%.4f" % (area / (heads / 2.0)))
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
# -b discards the first half, matching the equil_fraction membrane_analysis.py uses. Without
# it the MembraneDef would carry a whole-trajectory average while the verdict printed next to it
# reports the equilibrated one, and the two would disagree with no indication why.
PROD_HALF=$(python3 -c "print(int(${PROD_NS}*1000/2))")
printf 'Box-X\nBox-Y\n#Surf*SurfTen\n\n' \
    | "${GMX}" energy -f "${NAME}.edr" -o "${NAME}_energy.xvg" -b "${PROD_HALF}" \
      > "${NAME}_free.ENERGY" 2>&1

cd "${HERE}"
python3 tools/membrane_analysis.py single --name "${NAME}" --dir "${DIR}"
python3 tools/membrane.py calibrate --name "${NAME}" --role "${ROLE}" \
    --gro "${DIR}/${NAME}_eq.gro" --energy "${DIR}/${NAME}_free.ENERGY"

echo
echo "Paste the MembraneDef above into MEMBRANES in user/membranes.py, and copy"
echo "${DIR}/${NAME}_eq.gro and the matching .top into user/base/insert_peptide/ as"
echo "membrane_${NAME}.gro and membrane_${NAME}.top."
