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

    def test_detached_peptide_is_rejected_rather_than_scored_zero(self):
        # It used to score exactly 0.0, which on a maximised fitness outranks every attached
        # peptide with a negative raw score. Scaling cannot express "this measurement does not
        # mean what the objective assumes"; rejection can.
        obs = observables(5.0, 0.0)
        obs["bacterial"]["orientation"] = 0.0
        fitness, components = fitness_damage.fitness_from_observables(obs)
        self.assertEqual(fitness, UserSettings.DETACHED_FITNESS)
        self.assertIn("rejected", components)


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




class TestMetadataTolerance(unittest.TestCase):
    """
    Metadata keys in observables.json must not reach the objectives as membranes.

    compute_observables writes provenance alongside the per-membrane entries. Before the
    filtering helper existed, both objectives iterated over every key to find the minimum
    orientation factor, so a metadata entry would raise KeyError or, worse, contribute a
    meaningless value. This is the kind of break that only shows up against real output.
    """

    def setUp(self):
        UserSettings.APPLY_orientation_factor = True
        UserSettings.BINDING_proxy = "interaction_energy"

    def _with_metadata(self):
        obs = observables(2.0, 0.3)
        obs["_provenance"] = {"gmx_version": "2024.3", "modules": ["gromacs/2024.3"]}
        return obs

    def test_damage_ignores_metadata(self):
        plain, _ = fitness_damage.fitness_from_observables(observables(2.0, 0.3))
        tagged, _ = fitness_damage.fitness_from_observables(self._with_metadata())
        self.assertAlmostEqual(plain, tagged)

    def test_binding_ignores_metadata(self):
        plain, _ = fitness_binding.fitness_from_observables(observables(2.0, 0.3))
        tagged, _ = fitness_binding.fitness_from_observables(self._with_metadata())
        self.assertAlmostEqual(plain, tagged)

    def test_dispatch_ignores_metadata(self):
        for objective in ("damage", "binding"):
            plain, _ = score_observables(observables(2.0, 0.3), objective=objective)
            tagged, _ = score_observables(self._with_metadata(), objective=objective)
            self.assertAlmostEqual(plain, tagged, msg="objective=%s" % objective)


class TestOrientationWithNegativeRaw(unittest.TestCase):
    """
    A detached peptide must never outrank an attached one.

    The orientation factor is in [0, 1] and was applied as a plain multiplier. On a maximised
    fitness that inverts its purpose for any negative raw score: a peptide that had left the
    membrane (factor 0) scored exactly 0, beating every attached peptide that merely performed
    badly, so the optimiser was rewarded for producing peptides that detach. Upstream guarded
    against this by applying the factor only when the raw score was positive; the rewrite dropped
    the guard and these tests pin it back.
    """

    def setUp(self):
        UserSettings.APPLY_orientation_factor = True
        UserSettings.BINDING_proxy = "interaction_energy"
        UserSettings.DAMAGE_host_tolerance = 0.0
        UserSettings.DAMAGE_host_weight = 1.0

    def _pair(self, dg):
        attached = observables(dg, 0.0)
        detached = observables(dg, 0.0)
        for name in ("bacterial", "host"):
            detached[name]["orientation"] = 0.0
        return attached, detached

    def test_damage_detached_does_not_beat_attached(self):
        for dg in (-3.0, -1.0, -0.1):
            attached, detached = self._pair(dg)
            fa, _ = fitness_damage.fitness_from_observables(attached)
            fd_, _ = fitness_damage.fitness_from_observables(detached)
            self.assertLessEqual(fd_, fa, "detached outranked attached at delta_gamma=%s" % dg)

    def test_binding_detached_does_not_beat_attached(self):
        for dg in (-3.0, -1.0, -0.1):
            attached, detached = self._pair(dg)
            attached["bacterial"]["e_lipids"] = {"mean": 50.0, "err": 5.0}
            detached["bacterial"]["e_lipids"] = {"mean": 50.0, "err": 5.0}
            fa, _ = fitness_binding.fitness_from_observables(attached)
            fb, _ = fitness_binding.fitness_from_observables(detached)
            self.assertLessEqual(fb, fa)

    def test_penalty_still_applies_to_good_candidates(self):
        # The guard must not disable the orientation penalty where it belongs.
        attached, detached = self._pair(3.0)
        fa, _ = fitness_damage.fitness_from_observables(attached)
        fd_, _ = fitness_damage.fitness_from_observables(detached)
        self.assertGreater(fa, fd_)
        self.assertEqual(fd_, UserSettings.DETACHED_FITNESS)

    def test_a_detached_peptide_loses_to_every_attached_one(self):
        """
        The property the whole guard exists for, stated directly rather than pairwise.

        The earlier tests compared a detached and an attached peptide at the SAME raw score,
        which the old scaling already satisfied. The reachable failure was across raw scores: a
        detached peptide at raw +3 scored 0.0 and beat an attached one at raw -0.5, so in early
        generations, where most random sequences score negative, detaching was worth more than
        most honest attempts.
        """
        detached = observables(3.0, 0.0)
        for name in ("bacterial", "host"):
            detached[name]["orientation"] = 0.0
        worst_detached, _ = fitness_damage.fitness_from_observables(detached)

        for dg in (-5.0, -2.0, -0.5, -0.01):
            attached = observables(dg, 0.0)
            attached_fitness, _ = fitness_damage.fitness_from_observables(attached)
            self.assertLess(worst_detached, attached_fitness,
                            "a detached peptide outranked an attached one at delta_gamma=%s" % dg)

    def test_a_candidate_just_above_the_threshold_is_scaled_not_rejected(self):
        # The rejection must be a floor on the degenerate case, not a cliff that swallows
        # candidates which are merely imperfectly oriented.
        obs = observables(4.0, 0.0)
        obs["bacterial"]["orientation"] = UserSettings.ORIENTATION_reject_below + 0.01
        obs["host"]["orientation"] = 1.0
        fitness, components = fitness_damage.fitness_from_observables(obs)
        self.assertNotIn("rejected", components)
        self.assertGreater(fitness, UserSettings.DETACHED_FITNESS)

    def test_partial_orientation_scales_a_positive_score(self):
        obs = observables(4.0, 0.0)
        obs["bacterial"]["orientation"] = 0.5
        fitness, components = fitness_damage.fitness_from_observables(obs)
        self.assertAlmostEqual(components["orientation_factor"], 0.5)
        self.assertAlmostEqual(fitness, 0.5 * components["raw"])


class TestBindingProxyGuards(unittest.TestCase):

    def setUp(self):
        UserSettings.APPLY_orientation_factor = True

    def test_depth_proxy_reports_missing_depth(self):
        UserSettings.BINDING_proxy = "depth"
        obs = observables(1.0, 0.0)
        del obs["bacterial"]["d_mem"]
        with self.assertRaises(fitness_binding.BindingProxyUnavailable):
            fitness_binding.fitness_from_observables(obs)


if __name__ == "__main__":
    unittest.main(verbosity=2)


class TestOrientationIsAveragedOverFrames(unittest.TestCase):
    """
    The factor must be the mean of the per-frame factors, not the factor at the mean depth.

    The cosine window is strongly non-linear and symmetric about b, so the two differ most in
    exactly the cases the factor exists to catch. A peptide spending half the window buried
    between the leaflets and half dissolved in water has a mean depth of almost precisely b and
    was given 1.0: the maximum possible weight for a configuration that was never at the
    interface at all.

    compute_observables does the averaging over the real trajectory; this pins the arithmetic
    property that makes it necessary, so a later simplification back to factor(mean) fails here.
    """

    def setUp(self):
        from user.modules.compute_observables import orientation_factor
        self.factor = orientation_factor
        self.b = 1.75

    def _both(self, depths):
        import numpy as np
        depths = np.asarray(depths, dtype=float)
        at_mean = self.factor(float(depths.mean()), self.b)
        averaged = float(np.mean([self.factor(float(d), self.b) for d in depths]))
        return at_mean, averaged

    def test_half_buried_half_dissolved_is_not_a_perfect_score(self):
        at_mean, averaged = self._both([0.1] * 50 + [3.4] * 50)
        self.assertAlmostEqual(at_mean, 1.0, places=3,
                               msg="the setup no longer reproduces the mean-depth trap")
        self.assertLess(averaged, 0.05,
                        "a peptide never at the interface must not score near 1.0")

    def test_a_peptide_leaving_during_the_window_is_discounted(self):
        import numpy as np
        at_mean, averaged = self._both(np.linspace(0.25, 3.25, 100))
        self.assertGreater(at_mean, 0.99)
        self.assertLess(averaged, 0.7)

    def test_a_stable_peptide_is_essentially_unaffected(self):
        import numpy as np
        depths = np.random.default_rng(0).normal(self.b, 0.10, 200)
        at_mean, averaged = self._both(depths)
        self.assertGreater(averaged, 0.95)
        self.assertLess(abs(at_mean - averaged), 0.05)


class TestOnlyScoredRolesAffectTheFitness(unittest.TestCase):
    """
    user/membranes.py documents role="probe" as carried along for comparison and ignored by the
    objectives. The orientation minimum took every entry in the file, so adding one probe at a
    low orientation rescaled every fitness in the campaign, which looks exactly like the
    peptides having got worse.
    """

    def setUp(self):
        UserSettings.APPLY_orientation_factor = True
        UserSettings.DAMAGE_host_tolerance = 0.0
        UserSettings.DAMAGE_host_weight = 1.0
        self._registry = dict(Membranes.MEMBRANES)

    def tearDown(self):
        Membranes.MEMBRANES.clear()
        Membranes.MEMBRANES.update(self._registry)

    def _with_probe(self, orientation):
        obs = observables(2.0, 0.0)
        obs["probe1"] = {"role": "probe",
                         "delta_gamma": {"mean": 0.0, "err": 0.1},
                         "orientation": orientation,
                         "d_mem": {"mean": 1.0, "err": 0.1}}
        return obs

    def test_a_probe_membrane_does_not_scale_the_fitness(self):
        Membranes.MEMBRANES["probe1"] = Membranes.MembraneDef(
            name="probe1", gro="g", top="t", role="probe",
            area_nm2=64.0, surften_free=(0.0, 1.0), thickness_nm=4.0)
        without, _ = fitness_damage.fitness_from_observables(observables(2.0, 0.0))
        with_probe, _ = fitness_damage.fitness_from_observables(self._with_probe(0.10))
        self.assertAlmostEqual(without, with_probe)

    def test_a_probe_cannot_reject_a_candidate(self):
        Membranes.MEMBRANES["probe1"] = Membranes.MembraneDef(
            name="probe1", gro="g", top="t", role="probe",
            area_nm2=64.0, surften_free=(0.0, 1.0), thickness_nm=4.0)
        fitness, components = fitness_damage.fitness_from_observables(self._with_probe(0.0))
        self.assertNotIn("rejected", components)
        self.assertGreater(fitness, UserSettings.DETACHED_FITNESS)

    def test_the_target_still_counts(self):
        obs = observables(2.0, 0.0)
        obs["bacterial"]["orientation"] = 0.5
        obs["host"]["orientation"] = 1.0
        _fitness, components = fitness_damage.fitness_from_observables(obs)
        self.assertAlmostEqual(components["orientation_factor"], 0.5)

    def test_an_archived_run_falls_back_to_the_role_in_the_file(self):
        # Re-scoring results whose membranes are no longer in the registry: the role recorded
        # alongside the measurement is what decides, so a probe stays ignored.
        obs = self._with_probe(0.0)
        self.assertNotIn("probe1", Membranes.MEMBRANES)
        fitness, components = fitness_damage.fitness_from_observables(obs)
        self.assertNotIn("rejected", components)
        self.assertGreater(fitness, UserSettings.DETACHED_FITNESS)
