"""
Executable wrapper. Initialises the genetic algorithm and defines the fitness function.

Adapted from the Evo-MD run script of van Hilten et al. 2023. Differences from upstream:

  - the pipeline runs insert_peptide_nmem (N compositions) instead of insert_peptide_2mem
    (two tension states of one composition)
  - measurement and scoring are separate modules, so the objective is a setting rather than
    a code path
  - the Constraint class is actually wired to the GA. Upstream it was defined and never
    assigned to ptr_constraint_function, so the amphipathicity and charge constraints were
    dead code and only the default length and proline checks ever ran.
"""

import sys
import os

if not sys.version_info[0] == 3:
    sys.exit("Script requires Python 3.")

import core.settings            as CoreSettings
import core.simulation          as Sim
import core.system              as System
import core.genetic_algorithm   as GA
import core.provenance          as Provenance

import user.usersettings        as UserSettings
import user.membranes           as Membranes

from core.system                import mprint


def _setup_simulation_directory(sequence, path_parent_dir, dirname=None):
    if dirname is None:
        dirname = sequence.upper()
    sequence_dir = os.path.join(path_parent_dir, dirname)
    System.mkdir(sequence_dir)
    System.rm_dir_content(sequence_dir) # Ensure folder is empty in case of a rerun

    System.write_text_file(os.path.join(sequence_dir, UserSettings.filename_seq),
                           [">Test", sequence], add_newline=True)

    # Secondary structure handed to martinize2. The all-H assignment is the central structural
    # assumption of this approach: Martini cannot fold, so the helix is imposed rather than
    # predicted. For AMPs this is a stronger assumption than it was for curvature sensors,
    # since many are disordered in solution and only fold on binding. It is the first thing to
    # state as a limitation and the first thing to revisit if the calibration set misranks.
    System.write_text_file(os.path.join(sequence_dir, UserSettings.filename_ssd),
                           [str(len(sequence)), "H" * (len(sequence))], add_newline=True)
    return sequence_dir


def _compute_fitness_simulation(sequence, dirname=None):
    from user.modules.generate_peptide_PB  import module_generate_peptide
    from user.modules.insert_peptide_nmem  import module_insert_peptide_nmem
    from user.modules.production           import module_production_2mem
    from user.modules.compute_observables  import module_compute_observables
    from user.modules.compute_fitness      import module_compute_fitness

    sequence = "".join(sequence)
    sequence_dir = _setup_simulation_directory(sequence, CoreSettings.output_dir, dirname)

    N_retries = 0
    fitness = -1
    raw = 0.0

    while True:
        try:
            simulation = Sim.Simulation(sequence_dir)

            simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_seq))
            simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_ssd))

            simulation.add_module(Sim.Module(module_generate_peptide))
            simulation.add_module(Sim.Module(module_insert_peptide_nmem))
            # module_production_2mem is named for its origin but loops over membrane_types and
            # is already generic over the number of compositions.
            simulation.add_module(Sim.Module(module_production_2mem))
            simulation.add_module(Sim.Module(module_compute_observables))
            simulation.add_module(Sim.Module(module_compute_fitness))

            simulation.run()

            raw = simulation.get_output_variable("ddF")
            fitness = simulation.get_output_variable("fitness")

        except System.SimulationError:
            N_retries += 1
            if N_retries > CoreSettings.N_simulation_retry:
                print("Simulation stopped due to an error. Simulation retries exceeded. Exiting.")
                sys.exit(1)
            else:
                print("Simulation stopped due to an error. Retrying simulation. %d retries left."
                      % (CoreSettings.N_simulation_retry - N_retries))
            continue
        break

    return [fitness, raw]


def _run_specific_sequences(sequences):
    """
    Run a fixed list of sequences through the pipeline with no genetic algorithm.

    This is how the calibration set is measured: the same modules, the same membranes, the
    same sampling, applied to peptides whose behaviour is already known experimentally. Until
    that run shows the observables separating known actives from known inactives, nothing the
    GA produces can be interpreted.
    """
    # Runs serially when mpi4py is absent, rather than exiting as upstream did. The calibration
    # set is a few dozen peptides, which is slow without MPI but entirely possible, and on a
    # machine without root the MPI stack is the single most likely thing to be missing. The
    # genetic algorithm already degrades to serial the same way; this path should not be the
    # one hard requirement that blocks the whole pipeline.
    MPI = None
    try:
        from mpi4py import MPI
    except ImportError:
        pass

    if MPI is None:
        size, rank = 1, 0
        mprint("mpi4py not available: running %d sequences serially in this process."
               % len(sequences))
    else:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()

    buffer = [[] for _ in range(size)]

    count = 0
    while len(sequences) > 0:
        buffer[count].append(sequences.pop())
        count += 1
        if count >= size:
            count = 0

    queue = buffer[rank]
    for count, sequence in enumerate(queue):
        _compute_fitness_simulation(sequence, "".join(sequence) + str(rank) + "_" + str(count))
    print("Done.")


class Constraint:
    """
    Search-space constraints, composed by Constraint.combined.

    A constraint returns True when the sequence is NOT allowed. The GA discards such sequences
    before simulating them, so these are free compared with anything measured, and they are
    the right place to encode what the force field cannot represent and what the biology
    rules out.
    """

    @staticmethod
    def length(sequence):
        return len(sequence) < UserSettings.min_length or len(sequence) > UserSettings.max_length

    @staticmethod
    def terminal_proline(sequence):
        # martinize2 cannot assign a terminal proline; this is a tooling constraint and should
        # stay active regardless of the biology.
        return sequence[0] == "P" or sequence[-1] == "P"

    @staticmethod
    def charge_window(sequence):
        """Reject sequences outside the configured net charge window."""
        lo, hi = min(UserSettings.charge_range), max(UserSettings.charge_range)
        pos = sum(sequence.count(aa) for aa in ("K", "R"))
        neg = sum(sequence.count(aa) for aa in ("D", "E"))
        return not (lo <= pos - neg <= hi)

    @staticmethod
    def hydrophobic_fraction(sequence, max_fraction=0.65):
        """
        Cap the hydrophobic fraction.

        This is the cheapest available guard against the failure the curvature paper
        documented: left to maximise a membrane-interaction observable, the optimiser converges
        on sequences so hydrophobic that they would aggregate and precipitate out of solution
        rather than act on anything. A composition cap does not model solubility, it only
        keeps the search away from the region where solubility is obviously lost. A real
        aggregation term belongs in the objective, not here.
        """
        hydrophobic = ("A", "C", "F", "I", "L", "M", "V", "W", "Y")
        count = sum(sequence.count(aa) for aa in hydrophobic)
        return (count / float(len(sequence))) > max_fraction

    @staticmethod
    def combined(sequence):
        return (Constraint.length(sequence)
                or Constraint.terminal_proline(sequence)
                or Constraint.charge_window(sequence)
                or Constraint.hydrophobic_fraction(sequence))


def main():
    # Refuse to start if any membrane still carries placeholder reference values. Catching this
    # here costs a second; catching it after a campaign costs the campaign.
    Membranes.require_all_calibrated()

    # Record which GROMACS actually produced this campaign, and stop now if it is not the pinned
    # one. A campaign spanning several sessions can otherwise be produced by more than one build
    # without leaving any trace in the output.
    provenance = Provenance.capture()
    Provenance.check_expected(provenance)

    g = GA.GeneticAlgorithm()
    g.description = "AMP design run, objective=%s" % UserSettings.OBJECTIVE
    for key, value in Provenance.format_block(provenance):
        g.additional_output[key] = value

    g.ptr_scoring_function      = _compute_fitness_simulation
    g.population_size           = 144
    g.iterations                = 100
    g.N_parents                 = 36
    g.ptr_selection_function    = GA.SelectionMethods.best_N
    g.N_elites                  = 1
    g.N_rerun_elites            = 1
    g.N_rerun_threshold         = 4
    # Wired, unlike upstream.
    g.ptr_constraint_function   = Constraint.combined
    g.filepath_progress_output  = os.path.join(CoreSettings.output_dir, CoreSettings.name_output_progress)
    g.filepath_checkpoint       = os.path.join(CoreSettings.output_dir, CoreSettings.name_checkpoint_file)

    g.initialize(path_restart_file=CoreSettings.filepath_checkpoint)
    g.run()
    result = g.get_result()
    # get_result yields (sequence, [fitness, raw]); the second element is a list, not a scalar.
    best_sequence, best_scores = result[0][0], result[0][1]
    print("best\t%s\tfitness=%s\traw=%s"
          % (best_sequence, best_scores[0], best_scores[1] if len(best_scores) > 1 else ""))


if __name__ == "__main__":
    CoreSettings.parse_user_input()

    # 0: genetic algorithm
    # 1: fixed sequence list, no GA (calibration set); see also calibrate.py
    STATE = 0

    # WRITE_XTC does not control whether the trajectory is written: production.py always calls
    # MDRUN_XTC, so the .xtc always exists in the simulation's temp directory, which is where the
    # insertion-depth analysis reads it from. The flag only decides whether that file is copied
    # into the output directory and kept. So it is purely a disk-space question, and neither
    # stage needs a particular value of it to work.
    if STATE == 0:
        if UserSettings.WRITE_XTC:
            mprint("Refusing to start: WRITE_XTC is True, which keeps a trajectory for every "
                   "candidate. A campaign of thousands of evaluations will fill the disk. Set "
                   "WRITE_XTC = False in user/usersettings.py; the analysis does not need the "
                   "kept copies.")
            sys.exit(1)
        main()

    elif STATE == 1:
        mprint("Running fixed sequence list (no GA).")
        Provenance.check_expected()
        if not UserSettings.WRITE_XTC:
            mprint("Note: WRITE_XTC is False, so trajectories will not be kept after each run. "
                   "The observables are still computed from them; set it True only if you want "
                   "to inspect the calibration trajectories afterwards.")
        with open("seq.txt") as handle:
            sequences = [line.strip() for line in handle if line.strip()]
        _run_specific_sequences(sequences)
