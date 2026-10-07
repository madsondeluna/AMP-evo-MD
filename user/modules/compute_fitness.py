"""
Dispatch the configured objective over the observables produced by compute_observables.

This replaces the original module_compute_fitness_2mem_dict, which both measured and scored.
The split matters: measurement happens once per candidate, scoring is a pure function of the
result, and switching objectives therefore does not mean rerunning anything.

Writes a human-readable 'fitness' file alongside the machine-readable observables.json, and
publishes 'fitness' and the objective's components as simulation output variables so the GA
progress file records which objective produced a given run.
"""

import json
import os

import core.system as System

import user.usersettings as UserSettings
from user.modules import fitness_binding, fitness_damage


OBJECTIVES = {
    "damage": fitness_damage,
    "binding": fitness_binding,
}


def get_objective(name=None):
    name = name or UserSettings.OBJECTIVE
    try:
        return OBJECTIVES[name]
    except KeyError:
        raise ValueError("Unknown OBJECTIVE %r; expected one of %s"
                         % (name, ", ".join(sorted(OBJECTIVES))))


def score_observables(observables, objective=None, population_stats=None):
    """
    Score an observables dict without any simulation context.

    Used by the GA through module_compute_fitness, and directly by calibrate.py and by any
    later re-scoring of a finished campaign.
    """
    module = get_objective(objective)
    # Only the binding objective takes population statistics, and only for its combined proxy.
    if module is fitness_binding:
        return module.fitness_from_observables(observables, population_stats=population_stats)
    return module.fitness_from_observables(observables)


def module_compute_fitness(simulation):
    observables = simulation._output_variables.get("observables")

    if observables is None:
        # Fall back to the file, so the module still works when the pipeline is resumed or the
        # observables were produced by an earlier run.
        obs_path = os.path.join(simulation._temp_dir, "observables.json")
        with open(obs_path, "r") as handle:
            observables = json.load(handle)

    fitness, components = score_observables(observables)

    lines = ["objective=%s" % components["objective"],
             "fitness=%.6f" % fitness,
             "raw=%.6f" % components["raw"],
             "orientation_factor=%.6f" % components["orientation_factor"]]
    for key in sorted(components):
        if key in ("objective", "fitness", "raw", "orientation_factor"):
            continue
        lines.append("%s=%s" % (key, components[key]))

    filename = UserSettings.name_output_score
    System.write_text_file(os.path.join(simulation._temp_dir, filename), lines)

    print("\t".join(lines) + "\n")

    simulation._output_variables["fitness"] = fitness
    simulation._output_variables["fitness_components"] = components
    # Kept under its original name so downstream code and the GA progress file that expect a
    # 'ddF' column still find a value. It is no longer a curvature free energy: it is the raw,
    # pre-orientation score of whichever objective is active.
    simulation._output_variables["ddF"] = components["raw"]

    simulation._share.add_file("fitness", filename, True)
