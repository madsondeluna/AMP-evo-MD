"""
Tests for reading molecule charges out of Martini topologies.

This is load-bearing in a way that is easy to miss. TOP_get_system_charge decides how many
counterions GRO_TOP_neutralize_system adds, and it adds them by converting solvent beads. Get the
charge wrong and every production box runs charged under PME, which shifts the surface tension
and the interaction energies, which are the two observables the whole pipeline is built on. There
is no downstream check that can catch it: a charged box produces perfectly plausible numbers.

The defect these tests pin: a molecule was committed to the charge table only when a section
header that was neither moleculetype nor atoms came along. Every molecule followed directly by
the next [ moleculetype ] was dropped, which is the shape of every single-bead species. Of the 18
moleculetypes in the shipped Martini 3 ion library, 6 were parsed; the 12 lost included NA and
CL, the counterions this code exists to count.
"""

import glob
import io
import os
import re
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import core.gromacs as Gromacs

FORCE_FIELD_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               "user", "base", "insert_peptide")

MOLECULETYPE_RE = re.compile(r"^\s*\[\s*moleculetype\s*\]\s*\n(?:\s*;[^\n]*\n)*\s*(\S+)", re.M)


def charges_in(itp_text, directory=None):
    """Parse an itp given as text, returning {name: charge}."""
    directory = directory or tempfile.mkdtemp()
    itp = os.path.join(directory, "molecules.itp")
    with open(itp, "w") as handle:
        handle.write(itp_text)

    names = MOLECULETYPE_RE.findall(itp_text)
    top = os.path.join(directory, "system.top")
    with open(top, "w") as handle:
        handle.write('#include "molecules.itp"\n\n[ system ]\ntest\n\n[ molecules ]\n')

    result = {}
    for name in names:
        one = os.path.join(directory, "one.top")
        with open(one, "w") as handle:
            handle.write('#include "molecules.itp"\n\n[ system ]\ntest\n\n'
                         '[ molecules ]\n%s 1\n' % name)
        result[name] = Gromacs.TOP_get_system_charge(one, TOP_itp_directory_list=[directory])
    return result


class TestEveryMoleculeIsParsed(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_consecutive_moleculetypes_are_all_kept(self):
        # The exact shape that was dropped: single-bead molecules, one after another, with no
        # intervening section of any other kind.
        itp = (
            "[ moleculetype ]\nNA 1\n\n[ atoms ]\n1 Q5 1 NA NA 1 1.0\n\n"
            "[ moleculetype ]\nCL 1\n\n[ atoms ]\n1 Q5 1 CL CL 1 -1.0\n\n"
            "[ moleculetype ]\nCA 1\n\n[ atoms ]\n1 D 1 CA CA 1 2.0\n"
        )
        self.assertEqual(charges_in(itp, self.dir), {"NA": 1, "CL": -1, "CA": 2})

    def test_a_molecule_followed_by_other_sections_is_kept(self):
        itp = (
            "[ moleculetype ]\nPOPG 1\n\n"
            "[ atoms ]\n1 Q5 1 POPG PO4 1 -1.0\n2 N4a 1 POPG GL1 2 0.0\n\n"
            "[ bonds ]\n1 2 1 0.47 1250\n\n"
            "[ moleculetype ]\nW 1\n\n[ atoms ]\n1 P4 1 W W 1 0.0\n"
        )
        self.assertEqual(charges_in(itp, self.dir), {"POPG": -1, "W": 0})

    def test_the_last_molecule_in_a_file_is_kept(self):
        itp = "[ moleculetype ]\nNA 1\n\n[ atoms ]\n1 Q5 1 NA NA 1 1.0\n"
        self.assertEqual(charges_in(itp, self.dir), {"NA": 1})

    def test_comments_and_blank_lines_do_not_break_the_parse(self):
        itp = (
            "; a leading comment\n\n"
            "[ moleculetype ]\n; name  nrexcl\nNA 1\n\n"
            "[ atoms ]\n; id type resnr residue atom cgnr charge\n"
            "1 Q5 1 NA NA 1 1.0   ; trailing comment\n"
        )
        self.assertEqual(charges_in(itp, self.dir), {"NA": 1})

    def test_charges_accumulate_across_a_multi_bead_molecule(self):
        itp = ("[ moleculetype ]\nCDL2 1\n\n[ atoms ]\n"
               "1 Q5 1 CDL2 PO41 1 -1.0\n2 Q5 1 CDL2 PO42 2 -1.0\n"
               "3 N4a 1 CDL2 GL1 3 0.0\n")
        self.assertEqual(charges_in(itp, self.dir), {"CDL2": -2})


class TestShippedForceFieldFiles(unittest.TestCase):
    """
    The real files, because the synthetic cases above are what the author imagined and these are
    what the pipeline will actually read.
    """

    def _parse_whole_file(self, path):
        with open(path) as handle:
            text = handle.read()
        names = MOLECULETYPE_RE.findall(text)

        directory = tempfile.mkdtemp()
        try:
            shutil.copy(path, directory)
            top = os.path.join(directory, "system.top")
            with open(top, "w") as handle:
                handle.write('#include "%s"\n\n[ system ]\ntest\n\n[ molecules ]\n'
                             % os.path.basename(path))
                for name in names:
                    handle.write("%s 1\n" % name)

            captured = io.StringIO()
            stderr, sys.stderr = sys.stderr, captured
            try:
                Gromacs.TOP_get_system_charge(top, TOP_itp_directory_list=[directory])
            finally:
                sys.stderr = stderr
            unparsed = re.findall(r"(\S+) x 1", captured.getvalue())
        finally:
            shutil.rmtree(directory, ignore_errors=True)
        return names, unparsed

    def test_no_moleculetype_is_dropped_from_any_shipped_itp(self):
        files = sorted(glob.glob(os.path.join(FORCE_FIELD_DIR, "*.itp")))
        self.assertTrue(files, "no force field files found at %s" % FORCE_FIELD_DIR)

        checked = 0
        for path in files:
            names, unparsed = self._parse_whole_file(path)
            if not names:
                continue
            checked += len(names)
            self.assertEqual(unparsed, [],
                             "%s: %d of %d moleculetypes were not parsed: %s"
                             % (os.path.basename(path), len(unparsed), len(names),
                                sorted(unparsed)[:10]))
        self.assertGreater(checked, 300, "far fewer moleculetypes than expected were checked")

    def test_the_counterions_have_their_known_charges(self):
        # Against physics, not against the parser: these are the species the neutralisation
        # path adds, so a sign or magnitude error here is the one that charges the box.
        directory = tempfile.mkdtemp()
        try:
            for name in ("martini_v3.0.0_ions_v1.itp", "martini_v3.0.0_phospholipids_v1.itp",
                         "martini_v3.0.0_solvents_v2.itp"):
                shutil.copy(os.path.join(FORCE_FIELD_DIR, name), directory)
            includes = "".join('#include "%s"\n' % name for name in
                               ("martini_v3.0.0_ions_v1.itp",
                                "martini_v3.0.0_phospholipids_v1.itp",
                                "martini_v3.0.0_solvents_v2.itp"))

            for molecule, expected in (("NA", 1), ("CL", -1), ("CA", 2),
                                       ("W", 0), ("POPE", 0), ("POPG", -1), ("POPC", 0)):
                top = os.path.join(directory, "one.top")
                with open(top, "w") as handle:
                    handle.write(includes + "\n[ system ]\ntest\n\n[ molecules ]\n%s 1\n"
                                 % molecule)
                self.assertEqual(
                    Gromacs.TOP_get_system_charge(top, TOP_itp_directory_list=[directory]),
                    expected, "%s parsed with the wrong charge" % molecule)
        finally:
            shutil.rmtree(directory, ignore_errors=True)

    def test_a_neutral_salted_system_needs_no_counterions(self):
        """
        The end-to-end case that was wrong: a membrane built by insane with its own salt is
        already neutral, and used to compute as charged because the ions were invisible.
        """
        directory = tempfile.mkdtemp()
        try:
            for name in ("martini_v3.0.0_ions_v1.itp", "martini_v3.0.0_phospholipids_v1.itp",
                         "martini_v3.0.0_solvents_v2.itp"):
                shutil.copy(os.path.join(FORCE_FIELD_DIR, name), directory)
            top = os.path.join(directory, "system.top")
            with open(top, "w") as handle:
                handle.write(
                    "".join('#include "%s"\n' % name for name in
                            ("martini_v3.0.0_ions_v1.itp",
                             "martini_v3.0.0_phospholipids_v1.itp",
                             "martini_v3.0.0_solvents_v2.itp"))
                    + "\n[ system ]\ntest\n\n[ molecules ]\n"
                      "POPE 100\nPOPG 40\nW 2000\nNA 46\nCL 6\n")
            # POPG 40 x -1, NA 46 x +1, CL 6 x -1 = 0
            self.assertEqual(
                Gromacs.TOP_get_system_charge(top, TOP_itp_directory_list=[directory]), 0)
        finally:
            shutil.rmtree(directory, ignore_errors=True)


class TestMalformedTopologies(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _parse(self, itp_text):
        """
        Parse an itp directly, without going through the name-discovery helper above.

        The helper skips a file it finds no moleculetype in, which is exactly the malformed
        case here, so using it would make these tests pass without parsing anything.
        """
        itp = os.path.join(self.dir, "molecules.itp")
        with open(itp, "w") as handle:
            handle.write(itp_text)
        top = os.path.join(self.dir, "system.top")
        with open(top, "w") as handle:
            handle.write('#include "molecules.itp"\n\n[ system ]\ntest\n\n[ molecules ]\n')
        return Gromacs.TOP_get_system_charge(top, TOP_itp_directory_list=[self.dir])

    def test_an_atoms_section_before_any_moleculetype_is_named(self):
        with self.assertRaises(ValueError) as caught:
            self._parse("[ atoms ]\n1 Q5 1 NA NA 1 1.0\n")
        self.assertIn("before any", str(caught.exception))

    def test_a_short_atom_line_names_the_file_and_line(self):
        with self.assertRaises(ValueError) as caught:
            self._parse("[ moleculetype ]\nNA 1\n\n[ atoms ]\n1 Q5 1 NA\n")
        message = str(caught.exception)
        self.assertIn("line 5", message)
        self.assertIn("7 columns", message)


if __name__ == "__main__":
    unittest.main(verbosity=2)
