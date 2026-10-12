"""
Tests for resuming a campaign from a checkpoint.

A restart is used exactly when a campaign has already cost days, so a defect here destroys the
most expensive thing in the project. Three properties:

  the checkpoint on disk is always readable, whenever the process dies
  resuming does not destroy the history of the run it continues
  the resumed generations are distinguishable from the ones they continue

The second and third were broken together, and compounded: run() rewrote the progress file from
scratch on every start, including a restart into the documented output directory, so the
per-generation history was truncated. The checkpoint survived, so the run resumed correctly and
silently lost the only human-readable record of how it got there. Iteration numbering also
restarted at zero, so what did survive collided with itself in tools/ga_report.py, which keys
generations by that number.
"""

import os
import pickle
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "tools"))

import core.genetic_algorithm as GA
import run


def build_ga(directory, iterations=2, population=4, seed=1):
    import random
    ga = GA.GeneticAlgorithm()
    ga.description = "restart test"
    ga.ptr_scoring_function = lambda sequence: [float(len(set(sequence))), 1.0]
    ga.population_size = population
    ga.iterations = iterations
    ga.N_parents = 2
    ga.ptr_selection_function = GA.SelectionMethods.best_N
    ga.N_elites = 1
    ga.N_rerun_elites = 0
    ga.N_rerun_threshold = 4
    ga.ptr_constraint_function = run.Constraint.combined
    ga.filepath_progress_output = os.path.join(directory, "output.txt")
    ga.filepath_checkpoint = os.path.join(directory, "checkpoint.dat")
    ga.RNG = random.Random(seed)
    return ga


class TestRestartKeepsTheHistory(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.progress = os.path.join(self.dir, "output.txt")
        self.checkpoint = os.path.join(self.dir, "checkpoint.dat")

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _iteration_numbers(self):
        import ga_report
        _settings, iterations, _unparsed = ga_report.parse_progress(self.progress)
        return sorted(iterations)

    def test_a_first_run_writes_its_generations(self):
        ga = build_ga(self.dir, iterations=3)
        ga.initialize()
        ga.run()
        self.assertEqual(self._iteration_numbers(), [0, 1, 2])

    def test_resuming_does_not_truncate_the_progress_file(self):
        first = build_ga(self.dir, iterations=3)
        first.initialize()
        first.run()
        before = open(self.progress).read()

        second = build_ga(self.dir, iterations=2, seed=2)
        second.initialize(path_restart_file=self.checkpoint)
        second.run()
        after = open(self.progress).read()

        self.assertIn(before.strip(), after,
                      "the history of the campaign being continued was destroyed")
        self.assertGreater(len(after), len(before))

    def test_resumed_generations_get_new_numbers(self):
        first = build_ga(self.dir, iterations=3)
        first.initialize()
        first.run()
        second = build_ga(self.dir, iterations=2, seed=2)
        second.initialize(path_restart_file=self.checkpoint)
        second.run()

        numbers = self._iteration_numbers()
        self.assertEqual(numbers, [0, 1, 2, 3, 4],
                         "resumed generations collided with the ones they continue")
        self.assertEqual(len(numbers), len(set(numbers)))

    def test_the_report_tool_sees_every_generation_after_a_restart(self):
        # The property that matters downstream: the merge used to be invisible, because two
        # generations keyed 0 simply became one bucket with twice the candidates.
        import ga_report
        first = build_ga(self.dir, iterations=2)
        first.initialize()
        first.run()
        second = build_ga(self.dir, iterations=2, seed=2)
        second.initialize(path_restart_file=self.checkpoint)
        second.run()

        _settings, iterations, unparsed = ga_report.parse_progress(self.progress)
        self.assertEqual(unparsed, 0)
        self.assertEqual(len(iterations), 4)


class TestCheckpointFormat(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_the_checkpoint_records_the_iteration(self):
        ga = build_ga(self.dir, iterations=3)
        ga.initialize()
        ga.run()
        with open(os.path.join(self.dir, "checkpoint.dat"), "rb") as handle:
            loaded = pickle.load(handle)
        self.assertEqual(len(loaded), 4, "the iteration is not in the checkpoint")
        self.assertEqual(loaded[3], 2)

    def test_an_older_three_element_checkpoint_still_loads(self):
        # Checkpoints written before the iteration was recorded must not become unreadable:
        # that would turn a format change into exactly the loss this file exists to prevent.
        ga = build_ga(self.dir)
        ga.initialize()
        ga.run()
        path = os.path.join(self.dir, "checkpoint.dat")
        with open(path, "rb") as handle:
            pool, runs, population, _iteration = pickle.load(handle)
        with open(path, "wb") as handle:
            pickle.dump((pool, runs, population), handle)

        resumed = build_ga(self.dir, iterations=1, seed=3)
        resumed.initialize(path_restart_file=path)
        self.assertEqual(resumed.restart_from_iteration, 0)

    def test_no_temporary_file_is_left_behind(self):
        ga = build_ga(self.dir, iterations=2)
        ga.initialize()
        ga.run()
        leftovers = [name for name in os.listdir(self.dir) if name.endswith(".tmp")]
        self.assertEqual(leftovers, [])

    def test_the_checkpoint_is_replaced_atomically(self):
        """
        A kill during the write must leave the previous checkpoint intact, not a truncated one.

        Simulated by failing the write to the temporary file: the real checkpoint must be
        untouched, which is the whole point of writing beside it and moving it into place.
        """
        ga = build_ga(self.dir, iterations=1)
        ga.initialize()
        ga.run()
        path = os.path.join(self.dir, "checkpoint.dat")
        with open(path, "rb") as handle:
            good = handle.read()

        real_replace = os.replace
        os.replace = lambda *args, **kwargs: (_ for _ in ()).throw(OSError("disk full"))
        try:
            resumed = build_ga(self.dir, iterations=1, seed=4)
            resumed.initialize(path_restart_file=path)
            with self.assertRaises(OSError):
                resumed.run()
        finally:
            os.replace = real_replace

        with open(path, "rb") as handle:
            self.assertEqual(handle.read(), good,
                             "a failed write damaged the previous checkpoint")


if __name__ == "__main__":
    unittest.main(verbosity=2)
