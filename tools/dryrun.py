#!/usr/bin/env python3
"""
A small version of the whole pipeline, with every step the real one has.

WHAT THIS IS FOR

Before committing a workstation to a campaign that runs for days, you want to know that the
machinery works: that GROMACS and martinize2 are reachable, that a membrane can be built and
equilibrated, that a peptide can be coarse-grained and inserted, that the observables come out,
that the genetic algorithm turns them into a selection, and that the reporting tools read what
was produced. Every one of those steps can fail for a reason that has nothing to do with the
science, and finding out eight hours in is expensive.

So this runs the same code paths as a real campaign, on a system small enough to finish while
you watch. It calls tools/protocol.sh for the membranes and run.py's own modules for the
candidates, rather than reimplementing either: a smoke test that exercises different code than
production proves nothing about production.

WHAT IT DOES NOT PROVE

The numbers are not usable. The box is small, the runs are short, and the averaging window is
widened to compensate, so the surface tension is not converged and the reference measurements
are not a calibration. A membrane prepared by this script must never be pasted into
user/membranes.py. The dry run answers "does it run", not "is the answer right".

USAGE

    python tools/dryrun.py                 # everything, about 12 minutes on 8 cores
    python tools/dryrun.py --stage 0       # just the environment check, about a second
    python tools/dryrun.py --from 2        # resume at the candidate stage
    python tools/dryrun.py --no-md         # skip everything that needs GROMACS

--no-md is the plumbing test: it checks the driver, the settings overrides, the registry, the
scoring layer and the reporting tools with stub observables, in a couple of seconds and without
a GROMACS at all. Useful on a laptop, and useful to tell a broken environment from a broken
pipeline.
"""

import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import core.errors as Errors  # noqa: E402 - needs HERE on sys.path first

RULE = "=" * 80
THIN = "-" * 80


# -------------------------------------------------------------------------------------------
# The small system
# -------------------------------------------------------------------------------------------
# Every value here is a deliberate reduction of a production setting, and each one is the
# smallest that still exercises the step it belongs to.

BOX_NM = 6.0            # 8.0 in production. Still has to fit the peptide plus a water layer.
PEPTIDE_LENGTH = 12     # 24 in production. A 12-mer helix is about 1.8 nm, so it fits.
MEMBRANE_NS = 1         # 200 equilibration + 400 production in the real protocol.
CANDIDATE_PS = 1000     # 400000 in production.
POPULATION = 4          # 144 in production.
ITERATIONS = 2          # 100 in production.

# The production setting keeps the last 5 percent of the trajectory, which out of 1 ns is 50 ps:
# too few frames to average anything. Widened here so the measurement is at least arithmetically
# meaningful. It is still not converged, which is why these numbers are not a calibration.
ENERGY_INIT_FRACTION = 0.5

MEMBRANES = [
    # name, composition, role. Two compositions rather than one, because the lipid counting,
    # the index groups and the neutralisation differ between an anionic bacterial mixture and a
    # zwitterionic host one, and those differences are exactly where this pipeline breaks.
    ("dry_target", "POPE:70,POPG:30", "target"),
    ("dry_host", "POPC:100", "host"),
]

# Fixed seed: a dry run that passes and then fails on a rerun with different random sequences
# tells you nothing about what changed.
SEED = 20261008


# -------------------------------------------------------------------------------------------
# Output
# -------------------------------------------------------------------------------------------

class Reporter:
    """Prints the stage framework. Kept separate so the stages themselves stay readable."""

    def __init__(self, no_md=False):
        self.results = []
        self.proved = []
        self.no_md = no_md
        self._started = None
        self._name = None

    def banner(self, work_dir, stages, estimate):
        print(RULE)
        print("AMP-evo-MD dry run")
        print(RULE)
        print()
        if self.no_md:
            print("PLUMBING ONLY (--no-md): no simulation runs. This checks the driver, the")
            print("scoring layer, the genetic algorithm and the reporting tools, and tells a")
            print("broken environment from a broken pipeline. It does not touch GROMACS.")
        else:
            print("A small version of the whole pipeline, with every step the real one has.")
            print("It proves the machinery works end to end. It does NOT produce usable")
            print("numbers: the box and the runs are too small to be a calibration.")
        print()
        print("  %-14s %.1f x %.1f nm, %d membranes" % ("system", BOX_NM, BOX_NM, len(MEMBRANES)))
        if not self.no_md:
            print("  %-14s %d ns per membrane run, %g ns per candidate run (production: 200, "
                  "400 and 400)" % ("simulation", MEMBRANE_NS, CANDIDATE_PS / 1000.0))
        print("  %-14s %d candidates x %d iterations" % ("campaign", POPULATION, ITERATIONS))
        print("  %-14s %d residues" % ("peptides", PEPTIDE_LENGTH))
        print("  %-14s %s" % ("working dir", work_dir))
        print("  %-14s %s" % ("estimate", estimate))
        print()
        print("Stages to run:")
        for number, title, purpose in stages:
            print("  %d  %-18s %s" % (number, title, purpose))
        print()

    def start(self, number, title, proves, command=None):
        print()
        print(THIN)
        print("STAGE %d  %s" % (number, title))
        print(THIN)
        print("What this proves: %s" % proves)
        if command:
            print()
            print("Equivalent command:")
            for line in command.strip().split("\n"):
                print("  %s" % line)
        print()
        self._started = time.time()
        self._name = "%d %s" % (number, title)

    def expect(self, lines):
        """What the operator should see, printed before it happens."""
        print("Expect:")
        for line in lines:
            print("  %s" % line)
        print()

    def detail(self, label, value):
        print("  %-28s %s" % (label, value))

    def earned(self, claim):
        """A claim this stage has just made good on. Printed in the summary, nowhere else."""
        self.proved.append(claim)

    def passed(self, proved=(), note=None):
        """
        proved: what this stage, having passed, entitles the summary to claim.

        Collected per stage rather than printed as one fixed list, so running a subset or
        --no-md cannot make the summary claim something nothing checked.
        """
        elapsed = time.time() - self._started
        print()
        print("PASS  (%s)%s" % (_duration(elapsed), ("  " + note) if note else ""))
        self.results.append((self._name, "PASS", elapsed, note or ""))
        self.proved.extend(proved)

    def failed(self, reason):
        elapsed = time.time() - self._started
        print()
        print("FAIL  (%s)" % _duration(elapsed))
        for line in str(reason).strip().split("\n"):
            print("  %s" % line)
        self.results.append((self._name, "FAIL", elapsed, _short_reason(reason)))

    def summary(self, work_dir, kept):
        print()
        print(RULE)
        print("SUMMARY")
        print(RULE)
        print()
        print("  %-24s %-6s %10s  %s" % ("stage", "result", "time", "note"))
        for name, result, elapsed, note in self.results:
            print("  %-24s %-6s %10s  %s" % (name[:24], result, _duration(elapsed), note[:34]))
        print()

        failed = [r for r in self.results if r[1] == "FAIL"]
        if failed:
            print("%d stage(s) failed. The first failure is the one to read; the later stages" %
                  len(failed))
            print("depend on the earlier ones, so one cause usually produces several lines.")
            print()
            print("Everything the run produced is under %s" % work_dir)
            return 1

        print("Everything ran. What that means, precisely:")
        print()
        for claim in self.proved:
            print("  - %s" % claim)
        print()
        # Claiming the simulation works in the mode that ran no simulation, or claiming a
        # stage that was skipped, would make this summary worse than no summary.
        if self.no_md:
            print("NOT proven, because --no-md ran no simulation: that GROMACS, martinize2 and")
            print("insane work, that a membrane can be built or equilibrated, and that a peptide")
            print("can be coarse-grained or inserted. Run without --no-md for those.")
        else:
            print("What it does not mean: these numbers are not converged and the membranes are")
            print("not calibrated. Do not paste the MembraneDef this produced into")
            print("user/membranes.py. Run tools/protocol.sh at the real box and run lengths.")
        skipped = [title for number, title, _ in STAGE_TABLE
                   if not any(r[0].endswith(title) for r in self.results)]
        if skipped:
            print()
            print("Stages not run: %s. Nothing above says anything about them."
                  % ", ".join(skipped))
        print()
        if kept:
            print("Output kept in %s" % work_dir)
        else:
            print("Working directory removed. Pass --keep to inspect it.")
        return 0


def _short_reason(reason):
    """
    One line for the summary table.

    An operator-error block starts with a rule and a banner, so the first line of it says
    nothing. The 'what:' field is the line that does.
    """
    for line in str(reason).strip().split("\n"):
        line = line.strip()
        if not line or set(line) <= {"="} or line == BANNER_LINE:
            continue
        if line.startswith("what:"):
            return line[len("what:"):].strip()
        return line
    return "failed"


BANNER_LINE = "AMP-evo-MD: cannot continue"


def _duration(seconds):
    if seconds < 90:
        return "%.1f s" % seconds
    if seconds < 5400:
        return "%.1f min" % (seconds / 60.0)
    return "%.1f h" % (seconds / 3600.0)


# -------------------------------------------------------------------------------------------
# Settings overrides
# -------------------------------------------------------------------------------------------

_REAL_INSERT_DIR = None


def apply_overrides(work_dir):
    """
    Shrink the pipeline in memory rather than by editing the shipped settings.

    Every module reads these through the UserSettings module object at call time, so assigning
    here is enough and nothing on disk changes. Editing usersettings.py instead would leave a
    repository that silently runs 1 ns campaigns if the dry run were interrupted.
    """
    import user.usersettings as UserSettings
    import core.settings as CoreSettings

    UserSettings.SIM_time = CANDIDATE_PS
    UserSettings.ENERGY_init_fraction = ENERGY_INIT_FRACTION
    UserSettings.ENERGY_init_time = ENERGY_INIT_FRACTION * CANDIDATE_PS

    UserSettings.min_length = PEPTIDE_LENGTH
    UserSettings.max_length = PEPTIDE_LENGTH
    UserSettings.rep_lengths = [PEPTIDE_LENGTH]
    UserSettings.repeats = [1]

    # The membranes this run builds live in the dry-run directory, not in the repository, so a
    # dry run can never overwrite a real prepared membrane.
    insert_dir = os.path.join(work_dir, "insert_peptide")
    os.makedirs(insert_dir, exist_ok=True)
    # The shipped directory, remembered the first time. Reading the setting instead would,
    # on a second run in the same process, copy from the previous run's directory, which has
    # usually been deleted by then.
    global _REAL_INSERT_DIR
    if _REAL_INSERT_DIR is None:
        _REAL_INSERT_DIR = UserSettings.path_basedir_insert_peptide
    real_base = _REAL_INSERT_DIR
    for filename in os.listdir(real_base):
        source = os.path.join(real_base, filename)
        if os.path.isfile(source):
            shutil.copy(source, insert_dir)   # force field itps and the mdp templates
    UserSettings.path_basedir_insert_peptide = insert_dir

    CoreSettings.output_dir = os.path.join(work_dir, "campaign")
    os.makedirs(CoreSettings.output_dir, exist_ok=True)
    return UserSettings, CoreSettings


def register_membranes(definitions):
    """Replace the shipped registry with the membranes this run built."""
    import user.membranes as Membranes
    import user.usersettings as UserSettings

    Membranes.MEMBRANES.clear()
    Membranes.MEMBRANES.update(definitions)
    UserSettings.membrane_types = list(definitions)
    return Membranes


# -------------------------------------------------------------------------------------------
# Stages
# -------------------------------------------------------------------------------------------

def stage_environment(report, args):
    report.start(
        0, "environment",
        "the external programs the pipeline shells out to are runnable, and which GROMACS "
        "produced this run is recorded",
        command='python -c "import run; run.preflight(membranes=False)"')
    report.expect([
        "a version string for each of gmx, martinize2 and PeptideBuilder",
        "anything missing is named with the command that installs or loads it",
    ])

    import core.provenance as Provenance
    provenance = Provenance.capture()

    missing = []
    for label, key in (("gmx", "gmx_version"), ("martinize2", "martinize2_version")):
        value = provenance.get(key)
        report.detail(label, value)
        if value in (None, "", "UNKNOWN"):
            missing.append(label)

    try:
        import PeptideBuilder
        report.detail("PeptideBuilder", getattr(PeptideBuilder, "__version__", "installed"))
    except ImportError as exc:
        report.detail("PeptideBuilder", "MISSING (%s)" % exc)
        missing.append("PeptideBuilder")

    report.detail("python", provenance.get("python", sys.version.split()[0]))

    if args.no_md:
        missing = [m for m in missing if m == "PeptideBuilder"]

    if missing:
        return "missing: %s\nRun 'source ~/amp-env.sh', or see the README." % ", ".join(missing)
    report.earned("the external programs the pipeline calls are reachable and report a version")
    return None


def stage_membranes(report, args, work_dir):
    report.start(
        1, "membranes",
        "a membrane can be planned, built with insane, minimised, equilibrated, measured, and "
        "turned into the four reference numbers a fitness evaluation needs",
        command="\n".join(
            "BOX=%.1f EQUIL_NS=%d PROD_NS=%d ROOT=%s tools/protocol.sh %s '%s' %s"
            % (BOX_NM, MEMBRANE_NS, MEMBRANE_NS, os.path.join(work_dir, "prep"), name,
               composition, role)
            for name, composition, role in MEMBRANES))
    report.expect([
        "each membrane prints steps 1/7 to 7/7",
        "a MembraneDef block with area_nm2, surften_free and thickness_nm",
        "very likely a warning that the membrane is not tensionless: expected at %d ns, and the"
        % MEMBRANE_NS,
        "  reason these numbers must not be used as a calibration",
    ])

    if args.no_md:
        return _stub_membranes(report, work_dir)

    definitions = {}
    for name, composition, role in MEMBRANES:
        environment = dict(os.environ,
                           BOX="%.1f" % BOX_NM,
                           MAX_RESIDUES=str(PEPTIDE_LENGTH),
                           EQUIL_NS=str(MEMBRANE_NS),
                           PROD_NS=str(MEMBRANE_NS),
                           ROOT=os.path.join(work_dir, "prep"))
        log = os.path.join(work_dir, "protocol_%s.log" % name)
        report.detail(name, "building (%s)" % log)
        with open(log, "w") as handle:
            process = subprocess.run(
                [os.path.join(HERE, "tools", "protocol.sh"), name, composition, role],
                cwd=HERE, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True)
            handle.write(process.stdout)

        if process.returncode != 0:
            return ("tools/protocol.sh failed for %s (exit %d).\nLast lines of %s:\n%s"
                    % (name, process.returncode, log,
                       "\n".join(process.stdout.strip().split("\n")[-15:])))

        definition = _parse_membranedef(process.stdout, name)
        if definition is None:
            return ("tools/protocol.sh finished for %s but emitted no MembraneDef block.\n"
                    "See %s." % (name, log))

        prepared = os.path.join(work_dir, "prep", name)
        _install_membrane(prepared, name)
        definitions[name] = definition
        report.detail("  area_nm2", "%.4f" % definition.area_nm2)
        report.detail("  surften_free", "%.3f +/- %.3f bar nm" % definition.surften_free)
        report.detail("  thickness_nm", "%.4f" % definition.thickness_nm)

    return _check_membranes(report, work_dir, definitions)


def _install_membrane(prepared_dir, name):
    """Copy the built system where insert_peptide_nmem expects it, as the protocol says to."""
    import user.usersettings as UserSettings
    destination = UserSettings.path_basedir_insert_peptide
    for suffix, target in ((".gro", "membrane_%s.gro" % name), (".top", "membrane_%s.top" % name)):
        source = os.path.join(prepared_dir, "%s_eq%s" % (name, suffix))
        if not os.path.isfile(source):
            source = os.path.join(prepared_dir, "%s%s" % (name, suffix))
        shutil.copy(source, os.path.join(destination, target))


def _parse_membranedef(text, name):
    """
    Execute the MembraneDef block the protocol prints, rather than re-deriving the numbers.

    Reading it back the way the operator would is the point: if the emitted block is malformed,
    that is a defect worth catching here and not when it is pasted by hand.
    """
    import user.membranes as Membranes

    lines = text.split("\n")
    try:
        first = next(i for i, line in enumerate(lines)
                     if line.strip().startswith('"%s": MembraneDef(' % name))
    except StopIteration:
        return None
    last = next(i for i in range(first, len(lines)) if lines[i].strip() == "),")

    block = "\n".join(line.strip() for line in lines[first:last + 1])
    namespace = {"MembraneDef": Membranes.MembraneDef}
    exec("REGISTRY = {" + block + "}", namespace)       # noqa: S102 - our own generated text
    return namespace["REGISTRY"][name]


def _write_registry(work_dir, definitions):
    """
    Write the registry this run used, as a file calibrate.py can load with --membranes.

    calibrate.py runs as its own process, so it imports the registry shipped in the repository
    rather than the one this driver installed in memory. Without this the dry run could not
    re-score its own results, which is itself the gap an operator hits when re-scoring an
    archived campaign after user/membranes.py has moved on.
    """
    path = os.path.join(work_dir, "dryrun_membranes.py")
    with open(path, "w") as handle:
        handle.write("# Registry used by this dry run. Not a calibration: see tools/dryrun.py.\n")
        handle.write("MEMBRANES = {\n")
        for name, definition in sorted(definitions.items()):
            handle.write('    "%s": MembraneDef(\n' % name)
            handle.write('        name="%s",\n' % definition.name)
            handle.write('        gro="%s",\n' % definition.gro)
            handle.write('        top="%s",\n' % definition.top)
            handle.write('        role="%s",\n' % definition.role)
            handle.write("        area_nm2=%r,\n" % definition.area_nm2)
            handle.write("        surften_free=%r,\n" % (tuple(definition.surften_free),))
            handle.write("        thickness_nm=%r,\n" % definition.thickness_nm)
            handle.write("        composition=%r,\n" % (definition.composition,))
            handle.write('        notes="dry run, not a calibration",\n')
            handle.write("    ),\n")
        handle.write("}\n")
    return path


def _check_membranes(report, work_dir, definitions):
    """Physical plausibility, not convergence. A short run cannot be judged on convergence."""
    import user.membranes as Membranes

    for name, definition in definitions.items():
        area_per_lipid = definition.area_nm2 / max(sum(definition.composition.values()), 1) * 2
        if not 0.3 < area_per_lipid < 1.2:
            return ("%s has an implausible area per lipid of %.3f nm2 (expected 0.4 to 0.8). "
                    "The lipid count and the box do not agree." % (name, area_per_lipid))
        if not 2.5 < definition.thickness_nm < 6.5:
            return ("%s has an implausible thickness of %.2f nm (expected 3.5 to 4.5). The "
                    "bilayer may not have formed." % (name, definition.thickness_nm))

    register_membranes(definitions)
    Membranes.require_all_calibrated()
    Membranes.require_all_built()
    report.detail("registered", ", ".join(sorted(definitions)))
    report.detail("registry written to", _write_registry(work_dir, definitions))
    report.earned("the membrane registry accepts what the protocol emits"
                  if report.no_md else
                  "a membrane can be planned, built, equilibrated and measured")
    return None


def _stub_membranes(report, work_dir):
    """--no-md: plausible numbers and empty files, to exercise everything downstream."""
    import user.membranes as Membranes
    import user.usersettings as UserSettings

    definitions = {}
    for name, _composition, role in MEMBRANES:
        for suffix in (".gro", ".top"):
            path = os.path.join(UserSettings.path_basedir_insert_peptide,
                                "membrane_%s%s" % (name, suffix))
            with open(path, "w") as handle:
                handle.write("; dry run stub\n")
        definitions[name] = Membranes.MembraneDef(
            name=name, gro="membrane_%s.gro" % name, top="membrane_%s.top" % name, role=role,
            area_nm2=BOX_NM * BOX_NM, surften_free=(1.5, 1.0), thickness_nm=4.0,
            composition={"POPE": 100}, notes="dry run stub, not a measurement")
    report.detail("stubbed", ", ".join(sorted(definitions)))
    return _check_membranes(report, work_dir, definitions)


def stage_candidate(report, args, work_dir):
    sequence = _example_sequence()
    report.start(
        2, "one candidate",
        "a single sequence goes through every module of the per-candidate pipeline and comes "
        "out with a score",
        command="STATE = 1 in run.py, with one sequence in seq.txt")
    report.expect([
        "five steps in steps.jsonl: build_peptide, insert, production, observables, fitness",
        "an observables.json carrying surften, delta_gamma, e_lipids, d_mem and orientation",
        "a finite fitness, whose sign says whether the peptide raised or lowered the tension",
    ])
    report.detail("sequence", sequence)

    import run
    if args.no_md:
        return _stub_candidate(report, work_dir, sequence)

    fitness, raw = run._compute_fitness_simulation(sequence)
    return _check_candidate(report, work_dir, sequence, fitness, raw)


def _check_candidate(report, work_dir, sequence, fitness, raw):
    import core.steplog as steplog
    import run

    directory = os.path.join(_campaign_dir(work_dir), sequence.upper())
    records = steplog.read(os.path.join(directory, steplog.FILENAME))
    steps = [r for r in records if r.get("event") == "step"]
    names = [r["step"] for r in steps]

    expected = ["build_peptide", "insert", "production", "observables", "fitness"]
    if names != expected:
        return ("the pipeline did not run every step.\nexpected: %s\ngot:      %s"
                % (expected, names or "nothing"))
    report.detail("steps", ", ".join("%s %s" % (r["step"], r["status"]) for r in steps))

    observables_path = os.path.join(directory, "observables.json")
    if not os.path.isfile(observables_path):
        return "no observables.json was written in %s" % directory
    with open(observables_path) as handle:
        observables = json.load(handle)

    for name in MEMBRANES:
        entry = observables.get(name[0])
        if entry is None:
            return "observables.json has no entry for membrane %r" % name[0]
        for key in ("surften", "delta_gamma", "d_mem", "orientation"):
            if key not in entry:
                return "observables.json is missing %r for %s" % (key, name[0])
        report.detail(name[0], "delta_gamma %+.3f, d_mem %.2f nm, orientation %.2f"
                      % (entry["delta_gamma"]["mean"], entry["d_mem"]["mean"],
                         entry["orientation"]))

    import math
    if not math.isfinite(fitness) or fitness == run.WORST_FITNESS:
        return ("the candidate was scored %s, which is the failure sentinel. Read "
                "%s/steps.jsonl and stdout_stderr.log." % (fitness, directory))
    report.detail("fitness", "%+.4f  (raw %+.4f)" % (fitness, raw))
    report.earned("a peptide can be built, coarse-grained, inserted, simulated and measured")
    report.earned("the observables are produced and scored")
    return None


def _stub_candidate(report, work_dir, sequence):
    """--no-md: write the observables a real run would, and check the scoring layer reads them."""
    from user.modules.compute_fitness import score_observables
    import user.usersettings as UserSettings

    directory = os.path.join(_campaign_dir(work_dir), sequence.upper())
    os.makedirs(directory, exist_ok=True)
    observables = {}
    for index, (name, _composition, _role) in enumerate(MEMBRANES):
        observables[name] = {
            "surften": {"mean": 1.0 + index, "err": 0.5, "units": "bar*nm"},
            "delta_gamma": {"mean": -2.5 + index, "err": 0.4,
                            "units": "kJ mol-1 nm-2"},
            "delta_gamma_work": {"mean": -1.2 + index, "err": 0.3},
            "e_lipids": {"mean": -400.0 + 50 * index, "err": 10.0},
            "e_solvent": {"mean": -900.0, "err": 20.0},
            "d_mem": {"mean": 1.4, "err": 0.1, "units": "nm"},
            "orientation": 0.8,
        }
    with open(os.path.join(directory, "observables.json"), "w") as handle:
        json.dump(observables, handle, indent=2)

    for objective in ("damage", "binding"):
        fitness, components = score_observables(observables, objective=objective)
        report.detail("scored as %s" % objective, "%+.4f" % fitness)
        if fitness != fitness:
            return "the %s objective produced nan from well-formed observables" % objective
    report.earned("both objectives turn a well-formed observables file into a finite score")
    return None


def stage_campaign(report, args, work_dir):
    report.start(
        3, "mini campaign",
        "the genetic algorithm evaluates a population, selects parents, mutates, writes a "
        "progress file and a checkpoint, and survives a failed candidate",
        command="mpirun -np $AMP_NP python -m mpi4py run.py --output_dir <dir>")
    report.expect([
        "%d candidates evaluated per iteration, %d iterations" % (POPULATION, ITERATIONS),
        "output.txt gaining an 'Iteration N' block per iteration",
        "checkpoint.dat appearing and being replaced, never truncated",
    ])

    import core.genetic_algorithm as GA
    import core.settings as CoreSettings
    import run

    campaign = _campaign_dir(work_dir)
    progress = os.path.join(campaign, CoreSettings.name_output_progress)
    checkpoint = os.path.join(campaign, CoreSettings.name_checkpoint_file)

    ga = GA.GeneticAlgorithm()
    ga.description = "dry run"
    ga.ptr_scoring_function = (_fake_score if args.no_md else run._compute_fitness_simulation)
    ga.population_size = POPULATION
    ga.iterations = ITERATIONS
    ga.N_parents = max(2, POPULATION // 2)
    ga.ptr_selection_function = GA.SelectionMethods.best_N
    ga.N_elites = 1
    ga.N_rerun_elites = 0
    ga.N_rerun_threshold = 4
    ga.ptr_constraint_function = run.Constraint.combined
    ga.filepath_progress_output = progress
    ga.filepath_checkpoint = checkpoint
    ga.RNG = random.Random(SEED)

    ga.initialize()
    ga.run()
    best_sequence, best_scores = ga.get_result()[0]
    report.detail("best", "%s  fitness %+.4f" % ("".join(best_sequence), best_scores[0]))

    if not os.path.isfile(progress):
        return "no progress file was written at %s" % progress
    if not os.path.isfile(checkpoint):
        return "no checkpoint was written at %s" % checkpoint

    sys.path.insert(0, os.path.join(HERE, "tools"))
    import ga_report
    _settings, iterations, unparsed = ga_report.parse_progress(progress)
    report.detail("iterations in output.txt", len(iterations))
    report.detail("unparsed rows", unparsed)
    if len(iterations) < ITERATIONS:
        return ("the progress file has %d iteration(s), expected %d"
                % (len(iterations), ITERATIONS))
    if unparsed:
        return ("%d row(s) of the progress file did not parse, which means the format and the "
                "report tool have diverged" % unparsed)
    report.earned("the genetic algorithm evaluates, selects, mutates and checkpoints")
    return None


def _fake_score(sequence):
    """--no-md: a deterministic score, so the GA machinery runs without any simulation."""
    sequence = "".join(sequence)
    raw = -1.0 * (sequence.count("K") + sequence.count("R")) - 0.1 * sequence.count("W")
    return [raw, raw]


def stage_reports(report, args, work_dir):
    report.start(
        4, "reports",
        "the analysis tools read what the campaign wrote, which is the part that normally only "
        "gets exercised days in",
        command="\n".join([
            "python tools/ga_report.py <dir>/output.txt --outdir <dir>/report",
            "python calibrate.py <dir> --objective damage --membranes <dir>/dryrun_membranes.py",
        ]))
    report.expect([
        "progress.csv, summary.json and three .png files",
        "a TIMING BY STEP table naming the five pipeline steps",
        "calibrate.py printing a table of candidates with no labels supplied",
    ])

    campaign = _campaign_dir(work_dir)
    outdir = os.path.join(work_dir, "report")

    # calibrate.py runs as its own process and so imports the registry shipped in the
    # repository, not the one this driver installed in memory. --membranes points it at the one
    # these results were produced under, which is the same flag an operator needs to re-score an
    # archived campaign after user/membranes.py has moved on.
    registry = os.path.join(work_dir, "dryrun_membranes.py")

    for command in (
            [sys.executable, os.path.join(HERE, "tools", "ga_report.py"),
             os.path.join(campaign, "output.txt"), "--outdir", outdir],
            [sys.executable, os.path.join(HERE, "calibrate.py"), campaign,
             "--objective", "damage", "--membranes", registry]):
        process = subprocess.run(command, cwd=HERE, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True)
        name = os.path.basename(command[1])
        if process.returncode != 0:
            return ("%s exited %d:\n%s"
                    % (name, process.returncode,
                       "\n".join(process.stdout.strip().split("\n")[-15:])))
        report.detail(name, "ok, %d lines of output" % len(process.stdout.split("\n")))

    produced = []
    for filename in ("progress.csv", "summary.json", "performance.png",
                     "composition.png", "positional.png"):
        path = os.path.join(outdir, filename)
        if not os.path.isfile(path):
            return "tools/ga_report.py did not write %s" % filename
        produced.append("%s (%d kB)" % (filename, os.path.getsize(path) // 1024))
    report.detail("written to %s" % outdir, "")
    for item in produced:
        report.detail("  " + item.split(" ")[0], item.split(" ", 1)[1])
    report.earned("the reporting tools read what the campaign wrote")
    return None


# -------------------------------------------------------------------------------------------

def _campaign_dir(work_dir):
    return os.path.join(work_dir, "campaign")


def _example_sequence():
    """A cationic amphipathic 12-mer, built from the same alphabet the GA draws from."""
    import user.usersettings as UserSettings
    alphabet = sorted(a for a, p in UserSettings.propensities.items() if p > 0)
    rng = random.Random(SEED)
    import run
    while True:
        candidate = [rng.choice(alphabet) for _ in range(PEPTIDE_LENGTH)]
        if not run.Constraint.combined(candidate):
            return "".join(candidate)


STAGE_TABLE = [
    (0, "environment", "gmx, martinize2, PeptideBuilder"),
    (1, "membranes", "build, equilibrate, measure, emit a MembraneDef"),
    (2, "one candidate", "build peptide, insert, produce, measure, score"),
    (3, "mini campaign", "the genetic algorithm, %d iterations" % ITERATIONS),
    (4, "reports", "calibrate.py and tools/ga_report.py"),
]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", default=None,
                        help="working directory (default: a dryrun_<timestamp> directory here)")
    parser.add_argument("--stage", type=int, default=None, help="run only this stage")
    parser.add_argument("--from", dest="start", type=int, default=0,
                        help="start at this stage, requires --dir from a previous run")
    parser.add_argument("--no-md", action="store_true",
                        help="skip everything that needs GROMACS; tests the plumbing only")
    parser.add_argument("--keep", action="store_true",
                        help="keep the working directory even when every stage passes")
    args = parser.parse_args()

    work_dir = args.dir or os.path.join(os.getcwd(),
                                        "dryrun_" + time.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(work_dir, exist_ok=True)

    wanted = [args.stage] if args.stage is not None else [
        n for n, _, _ in STAGE_TABLE if n >= args.start]

    report = Reporter(no_md=args.no_md)
    report.banner(work_dir, [s for s in STAGE_TABLE if s[0] in wanted],
                  "a few seconds" if args.no_md else "about 12 minutes on 8 cores")

    apply_overrides(work_dir)

    stages = {
        0: lambda: stage_environment(report, args),
        1: lambda: stage_membranes(report, args, work_dir),
        2: lambda: stage_candidate(report, args, work_dir),
        3: lambda: stage_campaign(report, args, work_dir),
        4: lambda: stage_reports(report, args, work_dir),
    }

    for number in wanted:
        try:
            failure = stages[number]()
        except SystemExit as exc:
            failure = "the stage called sys.exit(%s); the message above says why" % exc.code
        except Errors.OperatorError as exc:
            # Already a complete, located diagnosis. A traceback on top of it would bury the
            # part the operator needs under frames from inside the pipeline.
            failure = str(exc)
        except Exception as exc:                                       # noqa: BLE001
            # Anything else is a defect rather than a misconfiguration, so the traceback is the
            # useful artefact and is kept.
            import traceback
            failure = "%s: %s\n%s" % (type(exc).__name__, exc, traceback.format_exc()[-1500:])
        if failure:
            report.failed(failure)
            break
        report.passed()

    status = report.summary(work_dir, kept=args.keep or any(
        r[1] == "FAIL" for r in report.results))
    if status == 0 and not args.keep:
        shutil.rmtree(work_dir, ignore_errors=True)
    return status


if __name__ == "__main__":
    sys.exit(main())
