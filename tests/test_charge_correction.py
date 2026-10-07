"""
Regression tests for the rewritten charge_correction.

Each case below is something the upstream implementation got wrong. They are kept as tests
rather than as a note in the commit message because the failure modes are silent: upstream did
not raise on the unreachable-window case, it looped, and the NameError only appeared for the
specific combination of a neutral sequence and a window that excludes zero.
"""

import random
import sys
import os
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.genetic_algorithm import _StringChromosome


# Ten helix-prone residues plus R, matching user/usersettings.py.
PROPENSITIES = {"A": 0.1, "C": 0.0, "D": 0.0, "E": 0.1, "F": 0.1, "G": 0.0, "H": 0.0,
                "I": 0.0, "K": 0.1, "L": 0.1, "M": 0.1, "N": 0.0, "P": 0.0, "Q": 0.1,
                "R": 0.1, "S": 0.1, "T": 0.0, "V": 0.0, "W": 0.1, "Y": 0.1}

# No cationic residue available, used to check the give-up path.
PROPENSITIES_NO_CATION = dict(PROPENSITIES, K=0.0, R=0.0)


def chromosome(sequence):
    return _StringChromosome(list(sequence))


class TestChargeCorrection(unittest.TestCase):

    def setUp(self):
        self.rng = random.Random(20261007)

    def test_neutral_sequence_into_cationic_window(self):
        # Upstream raised NameError here: net_charge == 0 entered the loop and neither of its
        # two unguarded ifs assigned charge_list.
        c = chromosome("ALAQSALAQSALAQSALAQSALAQ")
        self.assertEqual(c.net_charge(), 0)
        ok = c.charge_correction(1, self.rng, PROPENSITIES, [3, 9])
        self.assertTrue(ok)
        self.assertTrue(3 <= c.net_charge() <= 9)

    def test_can_raise_charge(self):
        # Upstream could only move charge toward zero, so a positive lower bound was
        # unsatisfiable regardless of the starting sequence.
        c = chromosome("EEALQSALAQSALAQSALAQSALA")
        self.assertLess(c.net_charge(), 3)
        ok = c.charge_correction(1, self.rng, PROPENSITIES, [4, 8])
        self.assertTrue(ok)
        self.assertTrue(4 <= c.net_charge() <= 8)

    def test_can_lower_charge(self):
        # The one direction upstream did handle; kept so the rewrite does not regress it.
        c = chromosome("KKKKKKKKKKKKALAQSALAQSAL")
        self.assertGreater(c.net_charge(), 9)
        ok = c.charge_correction(1, self.rng, PROPENSITIES, [3, 9])
        self.assertTrue(ok)
        self.assertTrue(3 <= c.net_charge() <= 9)

    def test_already_inside_window_is_untouched(self):
        c = chromosome("KALAQSKLAQSALAQRALAQSALA")
        before = list(c.sequence)
        self.assertTrue(3 <= c.net_charge() <= 9)
        ok = c.charge_correction(1, self.rng, PROPENSITIES, [3, 9])
        self.assertTrue(ok)
        self.assertEqual(before, c.sequence)

    def test_unreachable_window_gives_up_instead_of_looping(self):
        # No cationic residue in the alphabet and a window that requires one. Upstream would
        # spin forever; this must return False in bounded time.
        c = chromosome("ALAQSALAQSALAQSALAQSALAQ")
        ok = c.charge_correction(1, self.rng, PROPENSITIES_NO_CATION, [4, 8])
        self.assertFalse(ok)

    def test_length_is_preserved(self):
        for seq in ("ALAQSALAQSALAQSALAQSALAQ", "EEEEALAQSALAQSALAQSALAQ", "KKKKKKKKALAQSALAQSALAQS"):
            c = chromosome(seq)
            c.charge_correction(1, self.rng, PROPENSITIES, [3, 9])
            self.assertEqual(len(c.sequence), len(seq))

    def test_alphabet_is_respected(self):
        # Zero-weight residues must never be introduced.
        allowed = {aa for aa, w in PROPENSITIES.items() if w > 0}
        c = chromosome("EEEEEEEEALAQSALAQSALAQS")
        c.charge_correction(1, self.rng, PROPENSITIES, [4, 8])
        introduced = set(c.sequence) - set("EEEEEEEEALAQSALAQSALAQS")
        self.assertTrue(introduced <= allowed, "introduced disallowed residues: %s"
                        % (introduced - allowed))

    def test_repeat_factor_scales_the_window(self):
        # charge_range applies to the full sequence (block * repeat), not to the block.
        c = chromosome("KALAQSALAQSA")
        ok = c.charge_correction(2, self.rng, PROPENSITIES, [4, 8])
        self.assertTrue(ok)
        self.assertTrue(4 <= c.net_charge() * 2 <= 8)

    def test_many_random_starts_converge(self):
        alphabet = [aa for aa, w in PROPENSITIES.items() if w > 0]
        for _ in range(300):
            seq = [self.rng.choice(alphabet) for _ in range(24)]
            c = chromosome(seq)
            ok = c.charge_correction(1, self.rng, PROPENSITIES, [3, 9])
            self.assertTrue(ok, "failed to correct %s" % "".join(seq))
            self.assertTrue(3 <= c.net_charge() <= 9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
