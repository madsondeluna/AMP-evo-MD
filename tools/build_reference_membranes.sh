#!/usr/bin/env bash
#
# Canonical reference membranes for AMP-evo-MD.
#
# This is the standard, versioned recipe for the equilibrated membranes the pipeline is built
# on. Running it on any machine reproduces the same protocol: the composition, box, salt and
# equilibration length are pinned here, so the result is the same membrane within statistical
# error rather than an ad-hoc one.
#
# A molecular-dynamics trajectory is not bit-reproducible across machines or GROMACS builds, so
# the equilibrated structure itself is committed to the repository under
# user/base/insert_peptide/ as membrane_<name>.gro and membrane_<name>.top. That committed file
# is the one the per-candidate pipeline reuses; this script is how it was made and how to remake
# it elsewhere.
#
# Each entry runs pass 1 of tools/protocol.sh: build with insane, minimise, equilibrate at zero
# tension for EQUIL_NS, measure the equilibrium area per lipid, and emit plan_pass2.txt. The
# constant-area reference (pass 2 plus production) and copying the structure into place are the
# steps documented in the README, "Step 1: prepare and calibrate the membranes".
#
# Usage:
#   tools/build_reference_membranes.sh          # build every reference membrane below
#   NT=10 tools/build_reference_membranes.sh    # cap the GROMACS thread count
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Pinned parameters. These define the reference membranes; change them here, not at the call
# site, so that what the repository ships and what a rebuild produces cannot drift apart.
export BOX="${BOX:-8.0}"            # nm, lateral, both dimensions
export SALT="${SALT:-0.15}"         # mol/L, physiological, on top of neutralising counterions
export EQUIL_NS="${EQUIL_NS:-500}"  # zero-tension equilibration length, the reference standard
export AUTO=1                       # run the steps, do not stop after printing the plan

# name | composition | role. One row per reference membrane.
run_one() {
    local name="$1" composition="$2" role="$3"
    echo "=== reference membrane: ${name}  (${composition}, role ${role}) ==="
    "${HERE}/tools/protocol.sh" "${name}" "${composition}" "${role}"
}

# bacterial: Gram-negative inner-membrane mimic, anionic through POPG, buildable with the shipped
# Martini 3 lipid set. It is the target the damage objective evolves against. Cardiolipin (CDL2)
# is deferred until its Martini 3 parameters are added; see the README for why that matters.
run_one bacterial "POPE:70,POPG:30" target

# host: a mammalian-membrane mimic for the binding objective's selectivity term. A faithful host
# needs cholesterol and sphingomyelin, whose Martini 3 parameters are not shipped here, so it is
# left commented until they are added rather than built as a misleading neutral POPC bilayer.
# run_one host "POPC:100" host
