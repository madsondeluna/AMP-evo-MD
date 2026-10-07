"""
Tests for the two objectives over synthetic observables.

No GROMACS involved: the point of separating measurement from scoring is that the scoring is a
pure function of observables.json, so it can be tested on hand-written input. These tests fix
the sign conventions and the host-term behaviour, which are the two things easiest to get
backwards and hardest to notice once a campaign is running.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import user.usersettings as UserSettings
import user.membranes as Membranes
from user.modules import fitness_damage, fitness_binding
from user.modules.compute_fitness import score_observables


def entry(role, delta_gamma, e_lipids=-250.0, d_mem=1.75, orientation=1.0, area=45.0):
    return {
        "role": role,
        "area_nm2": area,
        "b_nm": 1.75,
        "delta_gamma": {"mean": delta_gamma, "err": 0.3},
        "delta_gamma_work": {"mean": delta_gamma * area, "err": 0.3 * area},
        "e_lipids": {"mean": e_lipids, "err": 5.0},
        "e_solvent": {"mean": -400.0, "err": 5.0},
        "d_mem": {"mean": d_mem, "err": 0.02},
        "orientation": orientation,
    }


def observables(dg_target, dg_host, **kwargs):
    return {
        "bacterial": entry("target", dg_target, **kwargs.get("target", {})),
        "host": entry("host", dg_host, **kwargs.get("host", {})),
    }


class TestDamageObjective(unittest.TestCase):

    def setUp(self):
        UserSettings.DAMAGE_host_tolerance = 0.0
        UserSettings.DAMAGE_host_weight = 1.0
        UserSettings.APPLY_orientation_factor = True

    def test_more_strain_on_target_scores_higher(self):
        weak, _ = fitness_damage.fitness_from_observables(observables(0.5, 0.0))
        strong, _ = fitness_damage.fitness_from_observables(observables(2.0, 0.0))
        self.assertGreater(strong, weak)

    def test_host_strain_is_penalised(self):
        clean, _ = fitness_damage.fitness_from_observables(observables(2.0, 0.0))
        dirty, _ = fitness_damage.fitness_from_observables(observables(2.0, 1.0))
        self.assertGreater(clean, dirty)

    def test_host_term_is_a_hinge_not_a_subtraction(self):
        # Below the tolerance the host term must not move the fitness at all. A symmetric
        # subtraction would reward peptides merely for being inert on the host, including
        # peptides inert on everything.
        UserSettings.DAMAGE_host_tolerance = 0.5
        a, _ = fitness_damage.fitness_from_observables(observables(2.0, 0.0))
        b, _ = fitness_damage.fitness_from_observables(observables(2.0, 0.4))
        self.assertAlmostEqual(a, b)
        c, _ = fitness_damage.fitness_from_observables(observables(2.0, 1.5))
        self.assertLess(c, a)

    def test_orientation_factor_takes_the_minimum(self):
        obs = observables(2.0, 0.0)
        obs["host"]["orientation"] = 0.25
        fitness, components = fitness_damage.fitness_from_observables(obs)
        self.assertAlmostEqual(components["orientation_factor"], 0.25)
        self.assertAlmostEqual(fitness, 0.25 * components["raw"])

    def test_detached_peptide_scores_zero(self):
        obs = observables(5.0, 0.0)
        obs["bacterial"]["orientation"] = 0.0
        fitness, _ = fitness_damage.fitness_from_observables(obs)
        self.assertAlmostEqual(fitness, 0.0)


class TestBindingObjective(unittest.TestCase):

    def setUp(self):
        UserSettings.BINDING_proxy = "interaction_energy"
        UserSettings.BINDING_host_weight = 1.0
        UserSettings.APPLY_orientation_factor = True

    def test_stronger_target_interaction_scores_higher(self):
        # Interaction energies are negative when bound, so more negative must score higher.
        weak, _ = fitness_binding.fitness_from_observables(
            observables(1.0, 0.0, target={"e_lipids": -100.0}))
        strong, _ = fitness_binding.fitness_from_observables(
            observables(1.0, 0.0, target={"e_lipids": -400.0}))
        self.assertGreater(strong, weak)

    def test_host_binding_reduces_selectivity(self):
        selective, _ = fitness_binding.fitness_from_observables(
            observables(1.0, 0.0, target={"e_lipids": -400.0}, host={"e_lipids": -50.0}))
        promiscuous, _ = fitness_binding.fitness_from_observables(
            observables(1.0, 0.0, target={"e_lipids": -400.0}, host={"e_lipids": -400.0}))
        self.assertGreater(selective, promiscuous)

    def test_depth_proxy_prefers_deeper_insertion(self):
        UserSettings.BINDING_proxy = "depth"
        shallow, _ = fitness_binding.fitness_from_observables(
            observables(1.0, 0.0, target={"d_mem": 2.4}, host={"d_mem": 2.4}))
        deep, _ = fitness_binding.fitness_from_observables(
            observables(1.0, 0.0, target={"d_mem": 1.5}, host={"d_mem": 2.4}))
        self.assertGreater(deep, shallow)

    def test_missing_energies_raise_rather_than_default(self):
        # A zero interaction energy is indistinguishable from a peptide that never touched the
        # membrane, so the absence must be loud.
        obs = observables(1.0, 0.0)
        obs["bacterial"]["e_lipids"] = None
        obs["bacterial"]["energy_warning"] = "energygrps not set"
        with self.assertRaises(fitness_binding.BindingProxyUnavailable):
            fitness_binding.fitness_from_observables(obs)

    def test_combined_proxy_refuses_without_population_stats(self):
        UserSettings.BINDING_proxy = "combined"
        with self.assertRaises(fitness_binding.BindingProxyUnavailable):
            fitness_binding.fitness_from_observables(observables(1.0, 0.0))


class TestDispatch(unittest.TestCase):

    def test_same_observables_score_under_both_objectives(self):
        # The whole reason for one repository: one measurement, two scalarisations.
        UserSettings.BINDING_proxy = "interaction_energy"
        obs = observables(2.0, 0.3)
        damage, dcomp = score_observables(obs, objective="damage")
        binding, bcomp = score_observables(obs, objective="binding")
        self.assertEqual(dcomp["objective"], "damage")
        self.assertEqual(bcomp["objective"], "binding")
        self.assertIsInstance(damage, float)
        self.assertIsInstance(binding, float)

    def test_unknown_objective_is_rejected(self):
        with self.assertRaises(ValueError):
            score_observables(observables(1.0, 0.0), objective="curvature")

    def test_rankings_can_disagree(self):
        # Constructed so that the peptide generating more strain is the weaker binder. That
        # the two objectives can disagree is the premise of running both; if they never could,
        # one of them would be redundant.
        UserSettings.BINDING_proxy = "interaction_energy"
        p1 = observables(3.0, 0.0, target={"e_lipids": -150.0}, host={"e_lipids": -50.0})
        p2 = observables(1.0, 0.0, target={"e_lipids": -500.0}, host={"e_lipids": -50.0})
        d1, _ = score_observables(p1, objective="damage")
        d2, _ = score_observables(p2, objective="damage")
        b1, _ = score_observables(p1, objective="binding")
        b2, _ = score_observables(p2, objective="binding")
        self.assertGreater(d1, d2)
        self.assertGreater(b2, b1)


class TestRegistryRoles(unittest.TestCase):

    def test_exactly_one_target(self):
        self.assertEqual(Membranes.the_target().name, "bacterial")

    def test_host_is_optional_but_unique(self):
        self.assertEqual(Membranes.the_host().name, "host")


if __name__ == "__main__":
    unittest.main(verbosity=2)
