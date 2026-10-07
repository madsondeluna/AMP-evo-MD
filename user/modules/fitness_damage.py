"""
Approach B: induced leaflet tension on the target membrane.

What is being maximised is delta_gamma on the bacterial membrane, with a hinge penalty on the
host membrane. delta_gamma = gamma(free) - gamma(with peptide) at constant area, so a positive
value means the peptide relieves tension in the leaflet it occupies, which is the same
statement as the peptide generating leaflet strain.

Why this is a damage observable and not a binding one. The mechanical reading in van Hilten
et al. 2023 is that a peptide which senses leaflet tension also induces it, and the paper's
discussion notes that membrane-binding peptides therefore carry an inherent destabilising
propensity: they lower the energetic cost of the highly curved interface of a toroidal pore.
The paper needs the derivative of this quantity with respect to area, because the curvature
selectivity it is after is a differential property. Damage is not differential, so the value
at a single area is what matters and the second area state is unnecessary. The cost is the
same two production runs; they now differ in composition rather than in tension.

Three limitations that belong in any claim made from this number, not in a footnote:

  Density. One peptide in a 128-lipid patch is about P:L 1:128. AMPs act cooperatively at
  ratios an order of magnitude higher. delta_gamma here is dilute-limit, per-peptide tension
  generation. Extrapolation to lysis is an assumption.

  Threshold. Generating tension is necessary for rupture, not sufficient. Nothing in this
  objective knows where the rupture threshold is. Establishing that requires long unbiased
  runs or biased sampling on the top candidates, outside the loop.

  Ensemble. The constant-area coupling in npt.mdp is load-bearing. Under a constant-tension
  (semi-isotropic) barostat the area would respond instead and gamma would be pinned near its
  reference value, and this observable would vanish. Do not change Pcoupltype without
  replacing the observable.
"""

import user.usersettings as UserSettings
import user.membranes as Membranes


def membrane_entries(observables):
    """
    The membrane measurements in an observables dict, excluding metadata.

    compute_observables stores provenance under a leading-underscore key alongside the per
    membrane entries. Filtering by prefix keeps the objectives working when further metadata is
    added later, which iterating over every key would not.
    """
    return {k: v for k, v in observables.items() if not k.startswith("_")}

def fitness_from_observables(observables):
    """
    Scalarise the observables dict into (fitness, components).

    Pure function of the JSON written by compute_observables, deliberately: it can be rerun
    over a finished campaign without touching GROMACS, which is what makes comparing the two
    objectives on one dataset a post-processing step.
    """
    target = Membranes.the_target()
    host = Membranes.the_host()

    target_obs = observables[target.name]
    dg_target = target_obs["delta_gamma"]["mean"]
    dg_target_err = target_obs["delta_gamma"]["err"]

    components = {
        "objective": "damage",
        "delta_gamma_target": dg_target,
        "delta_gamma_target_err": dg_target_err,
    }

    raw = dg_target

    # Host term as a hinge rather than a symmetric subtraction. Subtracting the host strain
    # outright would reward peptides that are merely inert on the host, including peptides
    # that are inert on everything; the hinge leaves fitness untouched until host strain
    # exceeds the tolerance and only then starts to cost. The tolerance is the knob that turns
    # this from "maximise selectivity" into "stay inside a window", which is the lesson the
    # curvature paper paid for: maximising a single physical quantity returns its degenerate
    # extreme, not a usable peptide.
    if host is not None:
        dg_host = observables[host.name]["delta_gamma"]["mean"]
        excess = max(0.0, dg_host - UserSettings.DAMAGE_host_tolerance)
        penalty = UserSettings.DAMAGE_host_weight * excess
        raw -= penalty
        components["delta_gamma_host"] = dg_host
        components["delta_gamma_host_err"] = observables[host.name]["delta_gamma"]["err"]
        components["host_penalty"] = penalty
        # Reported for interpretation only; the objective is not this difference.
        components["selectivity_gamma"] = dg_target - dg_host

    # Orientation factor, taken as the minimum across membranes as in the original: a peptide
    # that has left the surface of any of the systems is not measuring what this objective
    # means to measure, even if it behaves on the other.
    if UserSettings.APPLY_orientation_factor:
        factor = min(entry["orientation"] for entry in membrane_entries(observables).values())
    else:
        factor = 1.0

    components["orientation_factor"] = factor

    # The orientation factor is applied only to a positive raw score, restoring a guard the
    # upstream code had and this rewrite had dropped. Multiplying a negative raw score by a
    # factor in [0, 1] moves it toward zero, which on a maximised fitness makes it BETTER: a
    # peptide that had left the membrane entirely (factor 0) would score 0 and outrank every
    # attached peptide that merely performed badly. The optimiser would then be actively
    # rewarded for producing peptides that detach, and nothing in the fitness trace would show
    # it. A candidate whose raw score is already negative is already being rejected on the
    # merits, so the orientation penalty has nothing to add.
    if raw > 0:
        fitness = factor * raw
    else:
        fitness = raw
    components["raw"] = raw
    components["fitness"] = fitness

    return fitness, components
