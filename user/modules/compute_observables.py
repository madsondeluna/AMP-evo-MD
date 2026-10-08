"""
Extract every observable from the production trajectories and write them to observables.json.

This module computes no fitness. It separates measurement from scalarisation, for one
practical reason: the binding objective and the damage objective read the same .edr produced
by the same run, so there is no reason for either of them to own the extraction, and good
reason for neither to. Running both objectives over a finished campaign then becomes a
post-processing step on the JSON files rather than a second set of simulations.

Per membrane the following are recorded:

  surften        the #Surf*SurfTen column of the system with peptide, converted to
                 kJ mol-1 nm-2. Note that GROMACS reports the surface tension multiplied by
                 the number of interfaces, and a bilayer in a periodic box has two, so this is
                 2*gamma rather than gamma. The factor is carried through consistently and
                 matches the upstream convention, so the reference constants and the published
                 values are directly comparable; it is named here rather than divided out
                 because dividing would break that comparability.
  delta_gamma    surften(free) - surften(with peptide), kJ mol-1 nm-2
                 Positive means the peptide relieves tension in the leaflet it occupies.
                 At constant area that is the same statement as the peptide generating
                 leaflet strain, which is the mechanical step upstream of pore formation.
  delta_gamma_work
                 delta_gamma * area, kJ mol-1. Because delta_gamma already carries the factor
                 of two for the two interfaces, this is the change in the total interfacial
                 free energy of the box, which is a legitimate work-equivalent.

                 It is NOT comparable to the ddF of van Hilten et al. Theirs is a trapezoidal
                 integral of delta_gamma over the area DIFFERENCE between two tension states
                 (dA about 7 nm2); this is an absolute interfacial work at one area (about
                 64 nm2), so it is roughly an order of magnitude larger and means something
                 else. delta_gamma itself is intensive and is what the objectives use.
  e_lipids       Coul-SR + LJ-SR between peptide and lipids, kJ mol-1
  e_solvent      Coul-SR + LJ-SR between peptide and solvent, kJ mol-1
  d_mem          centre-of-mass distance between peptide and membrane along z, nm
  orientation    cosine factor from d_mem (Eq. 3 of van Hilten et al. 2023), in [0, 1]

Every quantity carries its standard error where one is available, because the averaging
window is short (see UserSettings.ENERGY_init_fraction) and a value without its error is not
interpretable here.

Sign conventions, which differ from the published figures: delta_gamma is positive for a
peptide that relieves tension, matching the internal convention of the original code, where
the fitness is maximised. The paper plots the corresponding quantity as negative and inverts
the y axis of Fig. 2B to reconcile the two. Do not compare these numbers to the published
ones without accounting for that.
"""

import json
import os

import numpy as np

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands
import core.provenance as Provenance

import user.usersettings as UserSettings
import user.membranes as Membranes


OBSERVABLES_FILENAME = "observables.json"

# Captured once per process rather than per candidate: running `gmx --version` for every
# evaluation would add a subprocess call to a loop that runs thousands of times, and the
# environment cannot change within a process.
_PROVENANCE = Provenance.capture()


def orientation_factor(d_mem, b):
    """
    Cosine window on insertion depth, Eq. 3 of van Hilten et al. 2023.

    Returns 1 when the peptide sits at the membrane surface (d_mem == b) and falls to 0 for a
    peptide that has dissolved into solution (d_mem >= 2b) or crossed into a transmembrane or
    inter-leaflet configuration (d_mem <= 0). Its purpose is to stop the optimiser from
    scoring configurations that are not what the objective means to measure.
    """
    if b is None or b <= 0:
        raise ValueError("orientation_factor needs a positive b (got %r)" % b)
    # Upper bound is 2b - a, where the cosine completes its period, not 2(b - a). The two
    # coincide only because UserSettings.a is 0.0; with a nonzero lower bound the factor would
    # jump discontinuously to zero while the cosine was still well away from it.
    if d_mem < UserSettings.a or d_mem > 2 * b - UserSettings.a:
        return 0.0
    return -0.5 * (np.cos((d_mem - UserSettings.a) * np.pi / (b - UserSettings.a)) - 1)


def module_compute_observables(simulation):

    results = {}

    for membrane_type in UserSettings.membrane_types:
        membrane = Membranes.get(membrane_type)
        membrane.require_calibrated()

        class LocalFiles():
            prod_XTC = simulation._share.get_filename("prod-%s-XTC" % membrane_type)
            prod_TPR = simulation._share.get_filename("prod-%s-TPR" % membrane_type)
            prod_EDR = simulation._share.get_filename("prod-%s-EDR" % membrane_type)
            # Per-membrane index file. Sharing one across compositions would silently read the
            # wrong atom indices; see the note in production.py.
            NDX      = simulation._share.get_filename("prod-%s-NDX" % membrane_type)

            surften_XVG  = UserSettings.name_output_production + "_" + membrane_type + "_surften.xvg"
            distance_XVG = UserSettings.name_output_production + "_" + membrane_type + "_distance.xvg"
            energy_XVG   = UserSettings.name_output_production + "_" + membrane_type + "_energy.xvg"

            SURFTEN_OUT = UserSettings.name_output_production + "_" + membrane_type + ".SURFTEN"
            ENERGY_OUT  = UserSettings.name_output_production + "_" + membrane_type + ".ENERGY"
            DIST_OUT    = UserSettings.name_output_production + "_" + membrane_type + ".DIST"

            stdout_stderr = None
            if CoreSettings.write_stdout_stderr_to_file:
                stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

        def path(filename):
            if filename is None:
                return None
            return os.path.join(simulation._temp_dir, filename)

        def gmx_energy(selection, path_xvg, path_out):
            """Run gmx energy with an interactive selection and return the parsed mean/err."""
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.ENERGY(
                    path(LocalFiles.prod_EDR),
                    path(path_xvg),
                    UserSettings.ENERGY_init_time,
                ),
                input=selection.encode(),
                path_stdout=path(path_out),
                path_stderr=path(path_out),
            )
            parsed = Gromacs.GMX_average_ENERGY_output(path(path_out))
            if parsed is None:
                # Raised as a SimulationError so run.py's retry logic engages. A bare TypeError
                # from unpacking None is not caught there, and a single truncated file would
                # otherwise kill an entire campaign instead of triggering a rerun.
                raise System.SimulationError(
                    "could not parse the statistics table from %s. The gmx energy output is "
                    "empty or truncated, which usually means the run was killed mid-write."
                    % path_out)
            return parsed

        def read_surften():
            mean, err = gmx_energy("#Surf*SurfTen\n",
                                   LocalFiles.surften_XVG, LocalFiles.SURFTEN_OUT)
            f = UserSettings.conversion_factor
            return mean * f, err * f

        def read_interaction_energies():
            """
            Short-range interaction energies between peptide and the lipid and solvent groups.

            Requires energygrps in npt.mdp (enabled in this fork). If the terms are absent the
            .edr simply does not contain them and gmx energy returns nothing usable, so a
            failure here is reported rather than silently defaulted: a zero interaction energy
            is indistinguishable from a peptide that never touched the membrane.
            """
            selection = "".join(term + "\n" for term in UserSettings.ENERGY_group_terms)
            try:
                _ = gmx_energy(selection, LocalFiles.energy_XVG, LocalFiles.ENERGY_OUT)
            except Exception as exc:
                return None, ("gmx energy failed for the interaction terms (%s). "
                              "Check that energygrps is set in npt.mdp." % exc)

            terms = _parse_energy_terms(path(LocalFiles.ENERGY_OUT),
                                        UserSettings.ENERGY_group_terms)
            if terms is None:
                return None, "could not parse interaction terms from the gmx energy output"
            # Every requested term, not merely one of them. A partial parse used to pass here and
            # _sum_terms would then skip the missing keys and return 0.0, which reaches the
            # binding objective as a perfectly plausible zero interaction energy. Half a sum is
            # worse than no sum, because nothing downstream can tell the difference.
            missing = [t for t in UserSettings.ENERGY_group_terms if t not in terms]
            if missing:
                return None, ("gmx energy output is missing %s (found %s). The group names in "
                              "the .edr do not match ENERGY_group_terms; check energygrps in "
                              "npt.mdp and the index groups."
                              % (", ".join(missing), ", ".join(sorted(terms))))
            return terms, None

        def read_distance():
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.DISTANCE(
                    path(LocalFiles.prod_XTC),
                    path(LocalFiles.prod_TPR),
                    path(LocalFiles.NDX),
                    path(LocalFiles.distance_XVG),
                    UserSettings.ENERGY_init_time
                ),
                input="com of group Protein plus com of group Lipids\n".encode(),
                path_stdout=path(LocalFiles.DIST_OUT),
                path_stderr=path(LocalFiles.DIST_OUT)
            )
            z_values = []
            for line in open(path(LocalFiles.distance_XVG), "r"):
                if "@" in line or "#" in line:
                    continue
                parts = line.split()
                if parts:
                    z_values.append(float(parts[-1]))
            if not z_values:
                raise RuntimeError("no distance frames parsed for membrane %r; the trajectory "
                                   "may be missing (WRITE_XTC off?)" % membrane_type)
            z = np.asarray(z_values)
            # Standard error of the mean, uncorrected for autocorrelation. It understates the
            # true uncertainty of a correlated MD time series and is reported as a floor, not
            # as a converged estimate.
            err = float(z.std(ddof=1) / np.sqrt(len(z))) if len(z) > 1 else 0.0
            return float(z.mean()), err

        # --- measurements ------------------------------------------------------------------
        surften, surften_err = read_surften()
        gamma_free, gamma_free_err = membrane.surften_free_converted()

        # Both terms carry the same n_surf factor, so it cancels in the difference as long as
        # the peptide perturbs one leaflet only. That is what the orientation factor enforces:
        # a transmembrane configuration would break the cancellation, and it is scored zero.
        delta_gamma = gamma_free - surften
        delta_gamma_err = float(np.sqrt(gamma_free_err ** 2 + surften_err ** 2))

        d_mem, d_mem_err = read_distance()
        factor = orientation_factor(d_mem, membrane.b_nm)

        energies, energy_warning = read_interaction_energies()

        entry = {
            "role": membrane.role,
            "area_nm2": membrane.area_nm2,
            "b_nm": membrane.b_nm,
            "surften": {"mean": surften, "err": surften_err, "units": "kJ mol-1 nm-2"},
            "surften_free": {"mean": gamma_free, "err": gamma_free_err, "units": "kJ mol-1 nm-2"},
            "delta_gamma": {"mean": delta_gamma, "err": delta_gamma_err, "units": "kJ mol-1 nm-2"},
            "delta_gamma_work": {"mean": delta_gamma * membrane.area_nm2,
                                 "err": delta_gamma_err * membrane.area_nm2,
                                 "units": "kJ mol-1"},
            "d_mem": {"mean": d_mem, "err": d_mem_err, "units": "nm"},
            "orientation": factor,
        }

        if energies is not None:
            e_lip = _sum_terms(energies, ["Coul-SR:Protein-Lipids", "LJ-SR:Protein-Lipids"])
            e_sol = _sum_terms(energies, ["Coul-SR:Protein-Solvent", "LJ-SR:Protein-Solvent"])
            entry["e_lipids"] = {"mean": e_lip[0], "err": e_lip[1], "units": "kJ mol-1"}
            entry["e_solvent"] = {"mean": e_sol[0], "err": e_sol[1], "units": "kJ mol-1"}
            entry["energy_terms"] = {k: {"mean": v[0], "err": v[1]} for k, v in energies.items()}
        else:
            entry["e_lipids"] = None
            entry["e_solvent"] = None
            entry["energy_warning"] = energy_warning

        results[membrane_type] = entry

        simulation._share.add_file("surften-%s-XVG" % membrane_type, LocalFiles.surften_XVG, False)
        simulation._share.add_file("distance-%s-XVG" % membrane_type, LocalFiles.distance_XVG, True)
        simulation._share.add_file("%s-SURFTEN" % membrane_type, LocalFiles.SURFTEN_OUT, True)
        simulation._share.add_file("%s-DIST" % membrane_type, LocalFiles.DIST_OUT, False)

    # --- write ----------------------------------------------------------------------------
    # One JSON per candidate. Across a campaign these files are the dataset: every observable
    # on every sequence, including the ones the active objective ignored. Keeping the
    # discarded observables is the cheapest thing in the pipeline and the hardest to recover
    # later.
    # Provenance travels with the measurements, not only in the campaign header, so a single
    # observables file is self-describing: re-scoring one later, or comparing two candidates
    # from different sessions, does not require trusting a note about which build was loaded.
    results["_provenance"] = {
        "gmx_version": _PROVENANCE.get("gmx_version"),
        "gmx_binary": _PROVENANCE.get("gmx_binary"),
        "martinize2": _PROVENANCE.get("martinize2_version"),
        "modules": _PROVENANCE.get("modules"),
    }

    out_path = os.path.join(simulation._temp_dir, OBSERVABLES_FILENAME)
    with open(out_path, "w") as handle:
        json.dump(results, handle, indent=2, sort_keys=True)

    simulation._output_variables["observables"] = results
    simulation._share.add_file("observables", OBSERVABLES_FILENAME, True)


def _parse_energy_terms(path_output, wanted):
    """
    Pull 'term  average  err.est.  rmsd  tot-drift' rows out of a gmx energy stdout dump.

    gmx prints the statistics table after the interactive selection; the term name can contain
    spaces and colons, so the row is matched by prefix against the requested term names rather
    than split blindly on whitespace.
    """
    try:
        lines = open(path_output, "r").read().splitlines()
    except OSError:
        return None

    found = {}
    for line in lines:
        stripped = line.strip()
        for term in wanted:
            if not stripped.startswith(term):
                continue
            rest = stripped[len(term):].split()
            if len(rest) < 2:
                continue
            try:
                found[term] = (float(rest[0]), float(rest[1]))
            except ValueError:
                continue
    return found or None


def _sum_terms(terms, keys):
    """
    Sum a set of energy terms, combining their errors in quadrature.

    Quadrature assumes the terms are independent, and Coul-SR and LJ-SR for the same group pair
    are not: they come from the same frames and are anticorrelated, so this understates the error
    on the sum. Summing frame by frame before averaging would be correct. Left as is because the
    binding objective uses the mean and not the error, but the reported error should not be
    quoted as if it were rigorous.
    """
    total = 0.0
    var = 0.0
    for key in keys:
        if key not in terms:
            continue
        mean, err = terms[key]
        total += mean
        var += err ** 2
    return total, float(np.sqrt(var))
