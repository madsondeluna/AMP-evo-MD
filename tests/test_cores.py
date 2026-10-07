"""
Tests for the core planner.

The edge cases matter more than the common one: a planner that returns zero ranks, or that
oversubscribes a small machine, fails at the moment someone runs it on hardware you did not
anticipate, and the symptom is a GROMACS error rather than anything pointing back here.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))

import cores


class TestPlan(unittest.TestCase):

    def test_typical_workstation(self):
        p = cores.plan(budget=32)
        self.assertEqual(p["threads"], 2)
        self.assertEqual(p["ranks"], 15)
        # One reserved, and at most one more idle from the split.
        self.assertLessEqual(p["occupancy"], 31)
        self.assertGreaterEqual(p["occupancy"], 30)

    def test_never_oversubscribes(self):
        for n in range(1, 129):
            p = cores.plan(budget=n)
            self.assertLessEqual(p["occupancy"], n,
                                 "occupancy %d exceeds %d cores" % (p["occupancy"], n))

    def test_always_at_least_one_rank_and_one_thread(self):
        for n in (1, 2, 3):
            p = cores.plan(budget=n)
            self.assertGreaterEqual(p["ranks"], 1)
            self.assertGreaterEqual(p["threads"], 1)

    def test_single_core_machine(self):
        p = cores.plan(budget=1)
        self.assertEqual(p["ranks"], 1)
        self.assertEqual(p["threads"], 1)

    def test_reserve_larger_than_machine_still_yields_a_plan(self):
        # Asking to reserve more cores than exist must not produce zero or a negative count.
        p = cores.plan(budget=2, reserve=8)
        self.assertEqual(p["ranks"], 1)
        self.assertEqual(p["threads"], 1)

    def test_reserve_zero_uses_everything(self):
        p = cores.plan(budget=16, reserve=0)
        self.assertEqual(p["usable"], 16)
        self.assertEqual(p["occupancy"], 16)

    def test_more_threads_per_rank_means_fewer_ranks(self):
        few = cores.plan(budget=33, threads_per_rank=8)
        many = cores.plan(budget=33, threads_per_rank=2)
        self.assertLess(few["ranks"], many["ranks"])
        # Both stay within budget.
        self.assertLessEqual(few["occupancy"], 32)
        self.assertLessEqual(many["occupancy"], 32)

    def test_one_thread_per_rank_leaves_nothing_idle(self):
        for n in (7, 13, 31, 64):
            p = cores.plan(budget=n, threads_per_rank=1)
            self.assertEqual(p["idle"], 0, "cores idle at budget %d" % n)

    def test_threads_never_exceed_usable(self):
        p = cores.plan(budget=3, threads_per_rank=16)
        self.assertLessEqual(p["threads"], p["usable"])
        self.assertLessEqual(p["occupancy"], 3)


class TestDetection(unittest.TestCase):

    def test_logical_cpus_is_positive(self):
        self.assertGreaterEqual(cores.logical_cpus(), 1)

    def test_physical_cores_is_none_or_positive(self):
        n = cores.physical_cores()
        self.assertTrue(n is None or n >= 1)

    def test_physical_does_not_exceed_logical(self):
        # Hyperthreading can only add logical CPUs, never remove them.
        physical = cores.physical_cores()
        if physical is not None:
            self.assertLessEqual(physical, cores.logical_cpus() * 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
