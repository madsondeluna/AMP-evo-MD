"""
Tests for how a failed candidate evaluation is handled.

core.system.check_exitcode calls sys.exit on any non-zero exit from a GROMACS or martinize2
command. Left uncaught, that ends a campaign that may have been running for days because one
sequence could not be built, and the loss is silent: the output directory simply stops growing.

A single failure should cost that candidate. An unbroken run of failures is not the sequences,
it is the setup, and continuing then produces a whole campaign of worst-possible scores that
looks like a completed run. These tests pin both halves.
"""

import os
import sys
import types
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import run


class _FakeSimulation:
    """Stands in for core.simulation.Simulation, failing on the evaluations we choose."""

    fail_on = set()
    counter = 0

    def __init__(self, directory):
        pass

    def add_input_file(self, path):
        pass

    def add_module(self, module):
        pass

    def run(self):
        _FakeSimulation.counter += 1
        if _FakeSimulation.counter in _FakeSimulation.fail_on:
            sys.exit(1)

    def get_output_variable(self, key):
        return 2.5 if key == "fitness" else 1.0


class TestFailedEvaluation(unittest.TestCase):

    def setUp(self):
        self._real_sim = run.Sim
        self._real_setup = run._setup_simulation_directory
        run.Sim = types.SimpleNamespace(Simulation=_FakeSimulation, Module=lambda f: f)
        run._setup_simulation_directory = lambda seq, parent, dirname=None: "/tmp"
        _FakeSimulation.counter = 0
        _FakeSimulation.fail_on = set()
        run._CONSECUTIVE_FAILURES = 0

    def tearDown(self):
        run.Sim = self._real_sim
        run._setup_simulation_directory = self._real_setup
        run._CONSECUTIVE_FAILURES = 0

    def test_single_failure_scores_the_candidate_and_continues(self):
        _FakeSimulation.fail_on = {1}
        first = run._compute_fitness_simulation("KWLRKAFKKLASKALQSLWSALAQ")
        self.assertEqual(first[0], run.WORST_FITNESS)
        second = run._compute_fitness_simulation("KWLRKAFKKLASKALQSLWSALAQ")
        self.assertEqual(second[0], 2.5)

    def test_failure_score_is_worse_than_any_plausible_candidate(self):
        # The genetic algorithm maximises fitness, so a failure must never be selectable as a
        # parent. Real delta_gamma values are of order 1 to 10 kJ mol-1 nm-2.
        self.assertLess(run.WORST_FITNESS, -1000)

    def test_failure_score_is_finite(self):
        # The progress file and tools/ga_report.py parse these as floats; an infinity would
        # propagate into every summary statistic and into the figures.
        import math
        self.assertTrue(math.isfinite(run.WORST_FITNESS))

    def test_success_resets_the_consecutive_counter(self):
        _FakeSimulation.fail_on = {1, 2}
        run._compute_fitness_simulation("AAAA")
        run._compute_fitness_simulation("AAAA")
        self.assertEqual(run._CONSECUTIVE_FAILURES, 2)
        run._compute_fitness_simulation("AAAA")
        self.assertEqual(run._CONSECUTIVE_FAILURES, 0)

    def test_scattered_failures_never_accumulate_to_an_abort(self):
        # Every other evaluation fails, far more than the threshold in total, but never
        # consecutively. The campaign must survive.
        _FakeSimulation.fail_on = set(range(1, 60, 2))
        for _ in range(50):
            run._compute_fitness_simulation("AAAA")
        self.assertLess(run._CONSECUTIVE_FAILURES, run.MAX_CONSECUTIVE_FAILURES)

    def test_systemic_failure_aborts(self):
        # A missing binary or a broken membrane fails everything. Continuing would produce a
        # campaign of worst-possible scores that looks like a completed run.
        _FakeSimulation.fail_on = set(range(1, 500))
        with self.assertRaises(SystemExit):
            for _ in range(run.MAX_CONSECUTIVE_FAILURES + 5):
                run._compute_fitness_simulation("AAAA")
        self.assertGreaterEqual(run._CONSECUTIVE_FAILURES, run.MAX_CONSECUTIVE_FAILURES)


if __name__ == "__main__":
    unittest.main(verbosity=2)
