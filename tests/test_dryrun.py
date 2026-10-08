"""
Tests for the dry run itself.

The dry run is what an operator trusts before committing a workstation to days of compute, so a
dry run that passes when something is broken is worse than no dry run. Two properties matter:

  it fails when the pipeline is broken, rather than reporting success
  it claims only what it actually checked

The second is easy to get wrong and invisible when it is: a summary listing "a membrane can be
built" after a run that built no membrane reads exactly like one that did.

These run the --no-md path, which exercises the driver, the settings overrides, the registry,
the scoring layer, the genetic algorithm and the reporting tools without needing GROMACS.
"""

import contextlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tools"))

import dryrun


@contextlib.contextmanager
def run_dryrun(*extra, patch=None):
    """
    Run the dry run as a subprocess, yielding (status, text, directory).

    A subprocess rather than a call into dryrun.main(), for two reasons. It is the command an
    operator actually runs, so this tests that rather than an approximation of it. And the
    driver deliberately mutates UserSettings and the membrane registry in the process it runs
    in, which is right for a driver and ruinous for a shared test process: in-process runs left
    min_length, SIM_time and the registry pointing at the dry run's values for every test that
    happened to run afterwards.

    patch: optional Python source executed inside the subprocess before main(), for the tests
    that need to break something and check the dry run notices.
    """
    directory = tempfile.mkdtemp()
    source = [
        "import sys",
        "sys.path.insert(0, %r)" % HERE,
        "sys.path.insert(0, %r)" % os.path.join(HERE, "tools"),
        "import dryrun",
        "sys.argv = %r" % (["dryrun.py", "--no-md", "--keep", "--dir", directory] + list(extra)),
        patch or "",
        "sys.exit(dryrun.main())",
    ]
    try:
        process = subprocess.run([sys.executable, "-c", "\n".join(source)],
                                 cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 text=True)
        yield process.returncode, process.stdout, directory
    finally:
        shutil.rmtree(directory, ignore_errors=True)


class TestDryRunPasses(unittest.TestCase):

    def test_every_stage_passes_without_gromacs(self):
        with run_dryrun() as (status, text, _directory):
            self.assertEqual(status, 0, "the --no-md dry run failed:\n" + text[-3000:])
            self.assertNotIn("FAIL", text)

    def test_it_runs_all_five_stages(self):
        with run_dryrun() as (_status, text, _directory):
            for number, title, _purpose in dryrun.STAGE_TABLE:
                self.assertIn("STAGE %d  %s" % (number, title), text)

    def test_it_writes_a_registry_its_own_calibrate_can_load(self):
        # calibrate.py runs as a separate process and would otherwise import the shipped
        # registry, which these results were not produced under.
        with run_dryrun() as (_status, _text, directory):
            registry = os.path.join(directory, "dryrun_membranes.py")
            self.assertTrue(os.path.isfile(registry))
            import user.membranes as Membranes
            original = dict(Membranes.MEMBRANES)
            try:
                loaded = Membranes.load_registry(registry)
                self.assertEqual(sorted(loaded), sorted(name for name, _c, _r in dryrun.MEMBRANES))
            finally:
                Membranes.MEMBRANES.clear()
                Membranes.MEMBRANES.update(original)


class TestDryRunIsHonest(unittest.TestCase):

    def test_no_md_does_not_claim_the_simulation_works(self):
        with run_dryrun() as (_status, text, _directory):
            summary = text.split("SUMMARY", 1)[1]
            for claim in ("a membrane can be planned, built",
                          "coarse-grained, inserted, simulated"):
                self.assertNotIn(claim, summary,
                                 "--no-md claimed something it did not check")
            self.assertIn("NOT proven", summary)

    def test_a_single_stage_claims_only_that_stage(self):
        with run_dryrun("--stage", "0") as (_status, text, _directory):
            summary = text.split("SUMMARY", 1)[1]
            self.assertIn("external programs", summary)
            self.assertNotIn("genetic algorithm evaluates", summary)
            self.assertIn("Stages not run:", summary)

    def test_the_banner_says_no_simulation_will_run(self):
        with run_dryrun() as (_status, text, _directory):
            self.assertIn("PLUMBING ONLY", text.split("STAGE 0", 1)[0])


class TestDryRunDetectsBreakage(unittest.TestCase):
    """The property that matters: a broken pipeline must not pass."""

    # A results set whose observables do not carry every registered membrane: what re-scoring
    # an archived campaign against a changed registry looks like.
    DROP_A_MEMBRANE = """
import json, os
_real = dryrun._stub_candidate
def _drop(report, work_dir, sequence):
    _real(report, work_dir, sequence)
    path = os.path.join(dryrun._campaign_dir(work_dir), sequence.upper(), "observables.json")
    with open(path) as handle:
        data = json.load(handle)
    data.pop(dryrun.MEMBRANES[1][0])
    from user.modules.compute_fitness import score_observables
    score_observables(data, objective="damage")
dryrun._stub_candidate = _drop
"""

    BREAK_SCORING = """
from user.modules import compute_fitness
def _broken(*args, **kwargs):
    raise RuntimeError("a distinctive failure string")
compute_fitness.score_observables = _broken
"""

    def test_a_missing_membrane_entry_fails_the_candidate_stage(self):
        with run_dryrun(patch=self.DROP_A_MEMBRANE) as (status, text, _directory):
            self.assertEqual(status, 1, "a broken scoring layer still passed:\n" + text[-2000:])
            self.assertIn("FAIL", text)
            self.assertIn("one candidate", text.split("SUMMARY", 1)[1])

    def test_a_scoring_failure_fails_rather_than_reporting_success(self):
        with run_dryrun(patch=self.BREAK_SCORING) as (status, _text, _directory):
            self.assertEqual(status, 1)

    def test_the_summary_note_names_the_cause(self):
        with run_dryrun(patch=self.BREAK_SCORING) as (_status, text, _directory):
            self.assertIn("distinctive failure string", text)


class TestShortReason(unittest.TestCase):
    """The summary table has one line per stage, so it has to pick the line that says something."""

    def test_an_operator_error_block_yields_its_what_field(self):
        import core.errors as Errors
        block = Errors.format_block("the thing that went wrong", where="/somewhere",
                                    why="because", fix="do this")
        self.assertEqual(dryrun._short_reason(block), "the thing that went wrong")

    def test_a_plain_message_yields_its_first_line(self):
        self.assertEqual(dryrun._short_reason("something broke\nand more detail"),
                         "something broke")


if __name__ == "__main__":
    unittest.main(verbosity=2)
