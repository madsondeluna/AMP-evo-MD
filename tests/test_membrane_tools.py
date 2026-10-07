"""
Tests for lipid counting and membrane planning.

The counting rule is the one piece of this tooling that is wrong in a way nobody would notice:
an inflated or deflated lipid count feeds the measured area per lipid, which feeds the second
planning pass, which sets the box the production runs are held at. A ten percent error there
produces a pre-stressed bilayer and every peptide is then scored against the wrong baseline,
with nothing in the output to say so.

Two counting mistakes have already been made here. Matching the bead name "PO4" exactly missed
cardiolipin entirely, because its beads are PO41 and PO42. Replacing that with a whitelist of
cardiolipin residue names then double-counted CDL1, which the whitelist omitted. The rule is now
by bead-name shape, and these tests pin every lipid family in the shipped library.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))

import membrane


def lipid(resname, beads, z=5.9):
    return [(resname, bead, z) for bead in beads]


class TestCountLipids(unittest.TestCase):
    """Every lipid counts exactly once, whatever its head group looks like."""

    CASES = [
        # name, beads as they appear in the Martini lipid library, expected count
        ("POPE", ["NH3", "PO4", "GL1", "GL2", "C1A", "C2A"], 1),
        ("POPG", ["GL0", "PO4", "GL1", "GL2", "C1A", "C2A"], 1),
        ("POPC", ["NC3", "PO4", "GL1", "GL2", "C1A", "C2A"], 1),
        ("POPS", ["CNO", "PO4", "GL1", "GL2", "C1A", "C2A"], 1),
        ("DPSM", ["NC3", "PO4", "AM1", "AM2", "C1A"], 1),
        # All three cardiolipin protonation states carry PO41 and PO42.
        ("CDL0", ["GL0", "PO41", "GL11", "GL21", "PO42", "GL12"], 1),
        ("CDL1", ["GL0", "PO41", "GL11", "GL21", "PO42", "GL12"], 1),
        ("CDL2", ["GL0", "PO41", "GL11", "GL21", "PO42", "GL12"], 1),
        # Sterols have no phosphate at all.
        ("CHOL", ["ROH", "R1", "R2", "R3", "R4", "R5"], 1),
    ]

    def test_each_lipid_counts_once(self):
        for name, beads, expected in self.CASES:
            counts, _zs = membrane.count_lipids(lipid(name, beads))
            self.assertEqual(counts.get(name, 0), expected,
                             "%s counted %d times" % (name, counts.get(name, 0)))

    def test_unlisted_multi_phosphate_lipid_counts_once(self):
        # The rule must not depend on knowing the residue name in advance, which is exactly how
        # CDL1 came to be double-counted.
        counts, _zs = membrane.count_lipids(lipid("XXXX", ["GL0", "PO41", "PO42"]))
        self.assertEqual(counts.get("XXXX"), 1)

    def test_mixture_totals(self):
        atoms = (lipid("POPE", ["NH3", "PO4"]) * 126
                 + lipid("POPG", ["GL0", "PO4"]) * 36
                 + lipid("CDL2", ["GL0", "PO41", "PO42"]) * 18)
        counts, zs = membrane.count_lipids(atoms)
        self.assertEqual(sum(counts.values()), 180)
        self.assertEqual(len(zs), 180)

    def test_sterol_rich_mixture(self):
        atoms = (lipid("POPC", ["NC3", "PO4"]) * 114
                 + lipid("DPSM", ["NC3", "PO4"]) * 20
                 + lipid("CHOL", ["ROH", "R1"]) * 56)
        counts, _zs = membrane.count_lipids(atoms)
        self.assertEqual(sum(counts.values()), 190)
        self.assertEqual(counts["CHOL"], 56)

    def test_glycerol_linker_is_not_counted(self):
        # GL0 is a linking bead on POPG and on cardiolipin; counting it would double every one.
        counts, _zs = membrane.count_lipids(lipid("POPG", ["GL0", "PO4"]))
        self.assertEqual(counts["POPG"], 1)


class TestMeasureThickness(unittest.TestCase):

    def _bilayer(self, species, upper_z=5.9, lower_z=2.0, per_leaflet=10):
        atoms = []
        for z in (upper_z, lower_z):
            for _ in range(per_leaflet):
                atoms.extend(lipid(*species, z=z))
        return atoms

    def test_plain_phospholipid(self):
        thickness, upper, lower = membrane.measure_thickness(
            self._bilayer(("POPE", ["NH3", "PO4"])))
        self.assertAlmostEqual(thickness, 3.9)
        self.assertEqual((upper, lower), (10, 10))

    def test_cardiolipin_contributes_one_phosphate(self):
        # Both phosphates sit at the same depth, so counting both would weight the lipid double.
        thickness, upper, lower = membrane.measure_thickness(
            self._bilayer(("CDL1", ["GL0", "PO41", "PO42"])))
        self.assertAlmostEqual(thickness, 3.9)
        self.assertEqual((upper, lower), (10, 10))

    def test_sterols_are_excluded(self):
        # Insertion depth is referenced to the phosphate plane, so a sterol must not shift it.
        atoms = self._bilayer(("POPC", ["NC3", "PO4"]))
        atoms.extend(lipid("CHOL", ["ROH"], z=4.5) * 20)
        thickness, upper, lower = membrane.measure_thickness(atoms)
        self.assertAlmostEqual(thickness, 3.9)
        self.assertEqual((upper, lower), (10, 10))

    def test_no_phosphates_fails_loudly(self):
        with self.assertRaises(SystemExit):
            membrane.measure_thickness(lipid("CHOL", ["ROH", "R1"]) * 10)


class TestMeanApl(unittest.TestCase):

    def test_mean_override_replaces_the_whole_mixture(self):
        comp = membrane.parse_composition("POPE:70,POPG:20,CDL2:10")
        value, measured = membrane.mean_apl(comp, {"MEAN": 0.711})
        self.assertAlmostEqual(value, 0.711)
        self.assertTrue(measured)

    def test_seed_values_are_not_reported_as_measured(self):
        comp = membrane.parse_composition("POPE:70,POPG:20,CDL2:10")
        _value, measured = membrane.mean_apl(comp, {})
        self.assertFalse(measured)

    def test_per_lipid_override_is_measured(self):
        comp = membrane.parse_composition("POPE:70,POPG:30")
        value, measured = membrane.mean_apl(comp, {"POPE": 0.5})
        self.assertTrue(measured)
        self.assertAlmostEqual(value, 0.7 * 0.5 + 0.3 * membrane.APL_SEED_NM2["POPG"])

    def test_unknown_lipid_without_a_seed_fails(self):
        comp = membrane.parse_composition("POPE:50,ZZZZ:50")
        with self.assertRaises(SystemExit):
            membrane.mean_apl(comp, {})


class TestParseComposition(unittest.TestCase):

    def test_fractions_normalise(self):
        comp = membrane.parse_composition("POPE:70,POPG:20,CDL2:10")
        self.assertAlmostEqual(sum(comp.values()), 1.0)
        self.assertAlmostEqual(comp["POPE"], 0.7)

    def test_relative_values_need_not_sum_to_100(self):
        comp = membrane.parse_composition("POPE:7,POPG:3")
        self.assertAlmostEqual(comp["POPE"], 0.7)

    def test_malformed_entry_raises(self):
        with self.assertRaises(ValueError):
            membrane.parse_composition("POPE")


if __name__ == "__main__":
    unittest.main(verbosity=2)
