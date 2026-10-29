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
import glob

if not sys.version_info[0] == 3:
    sys.exit("Script requires Python 3.")

import core.settings            as CoreSettings
import core.simulation          as Sim
import core.system              as System
import core.genetic_algorithm   as GA
import core.provenance          as Provenance
import core.errors              as Errors

import user.usersettings        as UserSettings
import user.membranes           as Membranes

from core.system                import mprint


# The twenty residues PeptideBuilder has geometries for. The check below exists because
# PeptideBuilder does not do it: Geometry.geometry() returns the glycine geometry for anything
# it does not recognise, including every lowercase letter, and raises nothing. A calibration set
# pasted from a FASTA download in lowercase would therefore run to completion as a set of
# polyglycine helices, produce no separation between the known classes, and be read as the
# observable having failed. That verdict decides whether the project continues, so this is
# checked before anything is built.
VALID_RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")


def _validate_sequence(sequence, source=None, line_number=None):
    """Return the sequence upper-cased, or stop naming the offending position."""
    upper = sequence.strip().upper()

    where = source
    if source is not None and line_number is not None:
        where = "%s, line %d" % (source, line_number)

    if not upper:
        Errors.fail("empty sequence", where=where,
                    why="every sequence must have at least one residue",
                    fix="remove the blank entry")

    bad = sorted({c for c in upper if c not in VALID_RESIDUES})
    if bad:
        positions = ", ".join("%r at position %d" % (c, i + 1)
                              for i, c in enumerate(upper) if c in bad)
        Errors.fail("sequence %s contains characters that are not standard residues: %s"
                    % (upper, ", ".join(repr(c) for c in bad)),
                    where=where if where else "sequence %s" % upper,
                    why=("PeptideBuilder has no geometry for these and silently substitutes "
                         "glycine rather than failing, so the peptide would be built wrong and "
                         "the run would look successful.\nOffending positions: %s" % positions),
                    fix="use the one-letter codes %s" % "".join(sorted(VALID_RESIDUES)))
    return upper


# A failed evaluation is scored rather than fatal, so one bad candidate does not end a campaign
# that has been running for days. The score has to be worse than any real candidate could reach,
# so the genetic algorithm never selects a failure as a parent; it is not -inf because the
# progress file and the report tooling parse these as floats and an infinity propagates into
# every summary statistic.
WORST_FITNESS = -1.0e6

# Consecutive failures before giving up. Scattered failures are candidates; an unbroken run of
# them is the setup.
MAX_CONSECUTIVE_FAILURES = 10
_CONSECUTIVE_FAILURES = 0


# Read relative to this file rather than to the current directory. The campaign is usually
# launched from a scratch or output directory, and resolving the calibration set against the
# repository is both what the operator means and what makes the failure reproducible.
# Default seq.txt, resolved against the repository. AMP_SEQ_FILE overrides it, which is how a
# subset (for example only the peptides added to an existing validation_run) is run without
# editing seq.txt or redoing the peptides already measured.
SEQUENCE_FILE = os.environ.get("AMP_SEQ_FILE") or os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "seq.txt")


def _read_sequence_file(path):
    """
    Read the calibration set: one sequence per line, blank lines and # comments ignored.

    Every error here names the file and the line, because this is operator-supplied input and
    the whole go/no-go verdict rests on it being what the operator thinks it is.
    """
    Errors.require_file(
        path,
        why=("STATE = 1 measures a fixed list of peptides whose behaviour is already known, "
             "one sequence per line. It is not shipped because the set is yours to choose."),
        fix=("cp %s %s\nthen edit it, or point SEQUENCE_FILE in run.py elsewhere"
             % (os.path.basename(path) + ".example", os.path.basename(path))))

    sequences = []
    seen = {}
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            text = line.split("#", 1)[0].strip()
            if not text:
                continue
            sequence = _validate_sequence(text, source=path, line_number=number)
            if sequence in seen:
                # Not fatal: the same peptide twice is a duplicated measurement, which is
                # wasteful rather than wrong, but it silently halves an apparent sample size.
                mprint("Note: %s on line %d of %s repeats line %d; it will be simulated twice."
                       % (sequence, number, path, seen[sequence]))
            else:
                seen[sequence] = number
            sequences.append(sequence)

    if not sequences:
        Errors.fail("%s contains no sequences" % os.path.basename(path),
                    where=os.path.abspath(path),
                    why="every line was blank or a comment",
                    fix="add one sequence per line, using one-letter residue codes")
    return sequences


def preflight(membranes=True):
    """
    Refuse to start on anything that would fail every candidate.

    All of these used to surface mid-run: a missing binary as a non-zero exit code inside the
    first candidate, a missing membrane file as a FileNotFoundError naming a path in a temporary
    directory that no longer exists by the time it is read. Each cost the operator a diagnosis
    out of proportion to the fault, and under MPI a worker dying this way leaves the master
    blocked in recv with nothing printed at all. Checking here costs about a second.

    Returns the provenance record, so the caller does not capture it twice.
    """
    missing = []

    provenance = Provenance.capture()
    if provenance.get("gmx_version") in (None, "", "UNKNOWN"):
        missing.append("GROMACS: could not run %r" % provenance.get("gmx_binary", "gmx"))
    if provenance.get("martinize2_version") in (None, "", "UNKNOWN"):
        missing.append("martinize2: not on PATH or not runnable")

    try:
        import PeptideBuilder  # noqa: F401
    except ImportError as exc:
        missing.append("PeptideBuilder: %s" % exc)

    if missing:
        Errors.fail("required programs are missing or not runnable",
                    where="\n".join(missing),
                    why=("every candidate calls these, so the campaign would fail on the first "
                         "one and keep failing"),
                    fix=("source ~/amp-env.sh  (loads the GROMACS module and the virtualenv)\n"
                         "then re-run tools/, or see 'Verify the installation' in the README"))

    # Settings that cannot produce a run, checked here rather than discovered by the first
    # candidate after it has already paid for its simulations.
    if UserSettings.OBJECTIVE == "binding" and UserSettings.BINDING_proxy == "combined":
        # The combined proxy z-scores each candidate against its own generation, which needs
        # population statistics. Nothing in the pipeline builds them: score_observables is
        # called per candidate with population_stats=None, so this proxy raises on the first
        # evaluation and has never been runnable.
        Errors.fail("BINDING_proxy = 'combined' is not wired up",
                    where="user/usersettings.py",
                    why=("it z-scores each candidate against its own generation, which needs "
                         "statistics over the whole population. Nothing collects those: the "
                         "objective is called one candidate at a time, so the proxy raises on "
                         "the first evaluation."),
                    fix="use BINDING_proxy = 'interaction_energy' or 'depth'")

    # Separate from the block above: a missing membrane is the operator's own build step not
    # having been copied into place, and the fix is a different one.
    if membranes:
        Membranes.require_all_calibrated()
        Membranes.require_all_built()

    Provenance.check_expected(provenance)
    return provenance


def _setup_simulation_directory(sequence, path_parent_dir, dirname=None):
    # Both entry points pass through here, so this is the one place where an unbuildable
    # sequence can be caught regardless of whether it came from seq.txt or from the GA.
    sequence = _validate_sequence(sequence)
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


def _compute_fitness_simulation(sequence, dirname=None, replica=1):
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
            simulation = Sim.Simulation(sequence_dir, sequence=sequence, replica=replica,
                                        retain_all=getattr(UserSettings, "RETAIN_ALL_TRAJECTORIES", False))

            simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_seq))
            simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_ssd))

            # Short names rather than the function names, because these are what a reader scans
            # in steps.jsonl and in the timing table.
            simulation.add_module(Sim.Module(module_generate_peptide, "build_peptide"))
            simulation.add_module(Sim.Module(module_insert_peptide_nmem, "insert"))
            # module_production_2mem is named for its origin but loops over membrane_types and
            # is already generic over the number of compositions.
            simulation.add_module(Sim.Module(module_production_2mem, "production"))
            simulation.add_module(Sim.Module(module_compute_observables, "observables"))
            simulation.add_module(Sim.Module(module_compute_fitness, "fitness"))

            simulation.run()

            raw = simulation.get_output_variable("ddF")
            fitness = simulation.get_output_variable("fitness")

        except (Errors.OperatorError, Errors.OperatorExit):
            # A misconfiguration with a known fix, already reported in full. It is not this
            # candidate's fault and it will not be any other candidate's fault either, so
            # retrying or scoring it would turn a one-line diagnosis into a campaign of
            # worst-possible scores that looks like a completed run. Listed before SystemExit
            # because OperatorExit is a subclass of it.
            mprint("Stopping: the failure above is a setup problem, not a property of %s."
                   % sequence)
            raise

        except System.SimulationError:
            N_retries += 1
            if N_retries > CoreSettings.N_simulation_retry:
                print("Simulation stopped due to an error. Simulation retries exceeded. Exiting.")
                sys.exit(1)
            else:
                print("Simulation stopped due to an error. Retrying simulation. %d retries left."
                      % (CoreSettings.N_simulation_retry - N_retries))
            continue

        except SystemExit:
            # core.system.check_exitcode calls sys.exit on any non-zero exit from a GROMACS or
            # martinize2 command, which would otherwise end a multi-day campaign because one
            # sequence failed. A single unmappable candidate should cost that candidate, not the
            # run.
            #
            # The consecutive counter is what keeps this from hiding a systemic fault: a missing
            # binary or a broken membrane topology fails every candidate, and that aborts after
            # the threshold rather than producing a campaign of worst-possible scores. Successful
            # evaluations reset it, so scattered failures never accumulate into a false abort.
            global _CONSECUTIVE_FAILURES
            _CONSECUTIVE_FAILURES += 1
            mprint("Evaluation of %s failed (a command exited non-zero). Scoring it as the worst "
                   "possible candidate and continuing. %d consecutive failures."
                   % (sequence, _CONSECUTIVE_FAILURES))
            if _CONSECUTIVE_FAILURES >= MAX_CONSECUTIVE_FAILURES:
                mprint("Aborting: %d evaluations in a row failed, which means the fault is in the "
                       "setup rather than in the sequences. Check %s in the most recent candidate "
                       "directory under %s."
                       % (_CONSECUTIVE_FAILURES, CoreSettings.filename_output_stdout_stderr,
                          CoreSettings.output_dir))
                raise
            return [WORST_FITNESS, 0.0]
        break

    _CONSECUTIVE_FAILURES = 0
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

    # Each sequence is run N_REPLICAS times as independent replicas. The replicas differ in two
    # ways: a distinct initial peptide orientation (a yaw about the membrane normal, set per
    # replica in the insertion module) and the random velocity seed (gen_seed=-1 in the mdps). This makes
    # the replicas genuinely independent starting points, not just reseeded copies, which is what
    # the single-replica validation lacked when near-identical peptides gave very different values.
    n_replicas = max(1, int(getattr(UserSettings, "N_REPLICAS", 1)))
    tasks = []
    for sequence in sequences:
        for replica in range(1, n_replicas + 1):
            tasks.append((sequence, replica))

    # Round-robin the (sequence, replica) tasks across MPI ranks. The directory name carries the
    # replica index, so it is unique across ranks without needing the rank in the name.
    buffer = [[] for _ in range(size)]
    for i, task in enumerate(tasks):
        buffer[i % size].append(task)

    queue = buffer[rank]
    for sequence, replica in queue:
        dirname = "".join(sequence) + "_rep" + str(replica)
        # Resume guard. A task that already wrote observables.json is complete, so skip it:
        # _setup_simulation_directory empties the task directory on entry, so without this a
        # re-run after an interruption would wipe and redo every replica instead of finishing
        # only what is missing. A partial directory (no observables.json) is not skipped and is
        # rebuilt cleanly. The recursive match covers observables.json in the task directory or
        # in its retained work/ copy.
        task_dir = os.path.join(CoreSettings.output_dir, dirname)
        if glob.glob(os.path.join(task_dir, "**", "observables.json"), recursive=True):
            mprint("skip %s: already complete (observables.json present)" % dirname)
            continue
        _compute_fitness_simulation(sequence, dirname, replica=replica)
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
    # Refuse to start on anything that would fail every candidate: a missing binary, a membrane
    # still carrying placeholder reference values, a membrane whose files were never copied into
    # place, or a GROMACS that is not the pinned one. Catching these here costs a second;
    # catching them after a campaign costs the campaign.
    provenance = preflight(membranes=True)

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
    # Overridable by the AMP_STATE environment variable, so the run mode is a launch-time choice
    # (AMP_STATE=1 for the validation) rather than an edit to this file. Keeping it out of the
    # file means the committed default stays 0 (campaign) and a git pull never conflicts on a
    # locally toggled STATE.
    STATE = int(os.environ.get("AMP_STATE", "0"))

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
        if getattr(UserSettings, "RETAIN_ALL_TRAJECTORIES", False):
            mprint("Refusing to start: RETAIN_ALL_TRAJECTORIES is True, which copies the entire "
                   "working directory (every intermediate trajectory) for every candidate. That "
                   "is for the validation, not a campaign of thousands of evaluations, which it "
                   "would fill the disk with. Set RETAIN_ALL_TRAJECTORIES = False in "
                   "user/usersettings.py, or decide a campaign retention policy first.")
            sys.exit(1)
        main()

    elif STATE == 1:
        mprint("Running fixed sequence list (no GA).")
        if not UserSettings.WRITE_XTC:
            mprint("Note: WRITE_XTC is False, so trajectories will not be kept after each run. "
                   "The observables are still computed from them; set it True only if you want "
                   "to inspect the calibration trajectories afterwards.")
        sequences = _read_sequence_file(SEQUENCE_FILE)
        mprint("Read %d sequences from %s." % (len(sequences), os.path.abspath(SEQUENCE_FILE)))
        preflight(membranes=True)
        _run_specific_sequences(sequences)
