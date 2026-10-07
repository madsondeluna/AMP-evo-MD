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
BOX="${BOX:-8.0}"              # lateral box in nm, both dimensions
MAX_RESIDUES="${MAX_RESIDUES:-30}"
SALT="${SALT:-0.15}"           # mol/L, on top of neutralising counterions
EQUIL_NS="${EQUIL_NS:-200}"    # semi-isotropic, to find the equilibrium area
PROD_NS="${PROD_NS:-400}"      # constant area, the reference run
ROOT="${ROOT:-prep}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT}/${NAME}"
mkdir -p "${DIR}"
cd "${DIR}"

MDP_DIR="${HERE}/user/base/production"

say() { printf '\n=== %s ===\n' "$1"; }

# ---------------------------------------------------------------------------
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
"${GMX}" mdrun -deffnm min

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
"${GMX}" mdrun -deffnm equil

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
    --apl "MEAN=${MEAN_APL}" | tee plan_pass2.txt || true

echo
echo "Pass 2 uses the measured area per lipid. Rebuild with its lipid counts, then re-run"
echo "steps 3 and 4, before continuing to production. Set PASS=2 to skip straight to"
echo "production once the rebuilt system has been equilibrated."

if [ "${PASS:-1}" != "2" ]; then
    exit 0
fi

# ---------------------------------------------------------------------------
say "6/7  production at constant area (the peptide-free reference run)"
# ---------------------------------------------------------------------------
cp equil.gro "${NAME}_eq.gro"
sed -e "s|\$INCLUDE_CWD|-I$(pwd)|" \
    -e "s|\$SIM_DT|0.03|" \
    -e "s|\$SIM_NSTEPS|$(python3 -c "print(int(${PROD_NS}*1000/0.03))")|" \
    "${MDP_DIR}/npt.mdp" > prod.mdp
"${GMX}" grompp -f prod.mdp -c "${NAME}_eq.gro" -p "${NAME}.top" -o "${NAME}.tpr" -maxwarn 2
"${GMX}" mdrun -deffnm "${NAME}"

# ---------------------------------------------------------------------------
say "7/7  extract, analyse, emit the MembraneDef"
# ---------------------------------------------------------------------------
printf 'Box-X\nBox-Y\n#Surf*SurfTen\n\n' \
    | "${GMX}" energy -f "${NAME}.edr" -o "${NAME}_energy.xvg" > "${NAME}_free.ENERGY" 2>&1

cd "${HERE}"
python3 tools/membrane_analysis.py single --name "${NAME}" --dir "${DIR}"
python3 tools/membrane.py calibrate --name "${NAME}" --role "${ROLE}" \
    --gro "${DIR}/${NAME}_eq.gro" --energy "${DIR}/${NAME}_free.ENERGY"

echo
echo "Paste the MembraneDef above into MEMBRANES in user/membranes.py, and copy"
echo "${DIR}/${NAME}_eq.gro and the matching .top into user/base/insert_peptide/ as"
echo "membrane_${NAME}.gro and membrane_${NAME}.top."
