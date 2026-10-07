"""
Approach A: partitioning selectivity, target membrane against host.

What is being maximised is how much more favourably the peptide associates with the bacterial
membrane than with the host mimic, using a cheap proxy for partitioning that is already
present in the production .edr.

The proxy is deliberately configurable (UserSettings.BINDING_proxy) because it should be
chosen from the calibration run rather than assumed. The three options are not equivalent and
none of them is the quantity one actually wants:

  interaction_energy
      Coul-SR + LJ-SR between peptide and lipids. Free in the .edr once energygrps is on.
      It is an enthalpy of interaction: no entropy, no desolvation cost, no accounting for the
      solvent the peptide gave up. It is known to behave poorly when comparing peptides of
      different composition, which is precisely the comparison a sequence optimiser makes.
      Use it knowing that, and check it against the calibration set before trusting a ranking.

  depth
      Centre-of-mass insertion depth. More robust than the interaction energy, since it does
      not pretend to be an energy at all, and correspondingly less informative: it saturates,
      and two peptides at the same depth can partition very differently.

  combined
      Both, z-scored within the generation so that neither unit dominates by magnitude.
      Z-scoring makes the fitness relative to the current population, which interacts badly
      with elitism and with the rerun averaging in the GA: a sequence's score then depends on
      the company it kept. Use only with that understood.

The correct quantity is the partitioning free energy, by thermodynamic integration, and it
does not fit inside the loop. The structure the curvature paper used applies unchanged here:
a cheap observable drives the evolution, and the expensive one re-ranks the survivors.
"""

import numpy as np

import user.usersettings as UserSettings
import user.membranes as Membranes


class BindingProxyUnavailable(RuntimeError):
    """Raised when the configured proxy cannot be computed from the observables at hand."""


def _interaction_energy(entry, membrane_name):
    if entry.get("e_lipids") is None:
        raise BindingProxyUnavailable(
            "interaction energies are missing for membrane %r (%s). The binding objective "
            "needs energygrps enabled in npt.mdp; the damage objective does not, so a run "
            "configured for damage cannot be re-scored for binding after the fact."
            % (membrane_name, entry.get("energy_warning", "no reason recorded"))
        )
    return entry["e_lipids"]["mean"]


def _affinity(entry, membrane_name, population_stats=None):
    """
    A scalar that increases with affinity for this membrane.

    Interaction energies are negative for a bound peptide, so the sign is flipped to keep the
    convention that larger means more strongly associated throughout this module.
    """
    proxy = UserSettings.BINDING_proxy

    if proxy == "interaction_energy":
        return -_interaction_energy(entry, membrane_name)

    if proxy == "depth":
        # Deeper insertion means a smaller centre-of-mass separation, hence the sign flip.
        return -entry["d_mem"]["mean"]

    if proxy == "combined":
        if population_stats is None:
            raise BindingProxyUnavailable(
                "BINDING_proxy='combined' z-scores against the current generation and needs "
                "population statistics, which are not available when scoring a single "
                "candidate in isolation. Use 'interaction_energy' or 'depth' for calibration "
                "runs and single-sequence rescoring."
            )
        e = -_interaction_energy(entry, membrane_name)
        d = -entry["d_mem"]["mean"]
        e_z = _z(e, population_stats["e_mean"], population_stats["e_std"])
        d_z = _z(d, population_stats["d_mean"], population_stats["d_std"])
        return (1.0 - UserSettings.BINDING_depth_weight) * e_z + UserSettings.BINDING_depth_weight * d_z

    raise ValueError("Unknown BINDING_proxy %r; expected interaction_energy, depth or combined"
                     % proxy)


def _z(value, mean, std):
    if std is None or std == 0 or not np.isfinite(std):
        return 0.0
    return (value - mean) / std


def membrane_entries(observables):
    """
    The membrane measurements in an observables dict, excluding metadata.

    compute_observables stores provenance under a leading-underscore key alongside the per
    membrane entries. Filtering by prefix keeps the objectives working when further metadata is
    added later, which iterating over every key would not.
    """
    return {k: v for k, v in observables.items() if not k.startswith("_")}

def fitness_from_observables(observables, population_stats=None):
    """
    Scalarise the observables dict into (fitness, components).

    Same signature contract as fitness_damage.fitness_from_observables, so the two are
    interchangeable at the call site and a finished campaign can be re-scored under either.
    """
    target = Membranes.the_target()
    host = Membranes.the_host()

    affinity_target = _affinity(observables[target.name], target.name, population_stats)

    components = {
        "objective": "binding",
        "proxy": UserSettings.BINDING_proxy,
        "affinity_target": affinity_target,
    }

    raw = affinity_target

    # Selectivity as a weighted difference. Unlike the damage objective this one is naturally
    # differential, and that difference carries a built-in brake the curvature paper's
    # objective lacked: hydrophobic insertion raises affinity for both membranes and largely
    # cancels in the subtraction, so what survives is electrostatic discrimination between an
    # anionic target and a zwitterionic, sterol-rich host.
    #
    # That does not remove the degeneracy, it moves it. The extreme of this objective is a
    # maximally cationic peptide: selective, weakly potent, salt-sensitive. A different
    # failure, still a failure, and still an argument for bounding the objective with a window
    # once the calibration set has told you where the window is.
    if host is not None:
        affinity_host = _affinity(observables[host.name], host.name, population_stats)
        raw -= UserSettings.BINDING_host_weight * affinity_host
        components["affinity_host"] = affinity_host
        components["selectivity"] = affinity_target - affinity_host

    if UserSettings.APPLY_orientation_factor:
        factor = min(entry["orientation"] for entry in membrane_entries(observables).values())
    else:
        factor = 1.0

    components["orientation_factor"] = factor
    fitness = factor * raw
    components["raw"] = raw
    components["fitness"] = fitness

    return fitness, components
