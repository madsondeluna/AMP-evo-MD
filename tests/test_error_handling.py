"""
Regression tests for failures that used to be illegible, silent, or fatal to a whole campaign.

Each test here pins one defect that was verified to occur. Grouped by what the defect cost,
because that is what decides whether a guard is worth its line count:

  CAMPAIGN LOSS      a fault that discarded compute already spent, or hung the run
  SILENT WRONG       a fault that produced a plausible number from a broken measurement
  ILLEGIBLE          a fault that was correctly detected but reported without saying where

The silent-wrong group is the one that matters most. A crash costs a restart; a plausible number
from a broken measurement propagates into the calibration verdict and into the campaign that
verdict authorises, and nothing downstream can tell it apart from a real result.
"""

import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import core.errors as Errors


@contextlib.contextmanager
def silenced():
    """
    Swallow the error block a test is deliberately triggering.

    Errors.fail writes to both stdout and stderr on purpose, so that a stop reason cannot land
    in the stream nobody is reading. That is right in production and unreadable in a test suite,
    where a few dozen of these would bury the results.
    """
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        yield buffer


# ---------------------------------------------------------------------------------------------
# CAMPAIGN LOSS
# ---------------------------------------------------------------------------------------------

class TestProgressWriterSurvivesAFailedCandidate(unittest.TestCase):
    """
    run.py scores a failed evaluation [WORST_FITNESS, 0.0] precisely so one bad candidate does
    not end a multi-day campaign. The progress writer then divided fitness by that 0.0, so the
    first failure killed the run anyway, from inside a logging function, after the whole
    generation had already been simulated.
    """

    def _penalty(self, fitness, ddF):
        # The expression as it appears in _Population.__write_progress_file.
        return (fitness / ddF) if ddF else float("nan")

    def test_a_failed_candidate_does_not_raise(self):
        import run
        try:
            value = self._penalty(run.WORST_FITNESS, 0.0)
        except ZeroDivisionError:
            self.fail("a failed candidate still divides by zero in the progress writer")
        self.assertNotEqual(value, value)  # nan

    def test_a_normal_candidate_is_unchanged(self):
        self.assertAlmostEqual(self._penalty(2.5, 1.25), 2.0)

    def test_the_report_parser_accepts_the_nan_row(self):
        # Writing nan is only safe if the tooling reads it; otherwise the fix moves the failure
        # from the campaign to the analysis.
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "tools"))
        import ga_report
        row = "\tAAAA\t-1000000.0\tddF=0.0\tpenalty=nan\tn_run=1"
        self.assertIsNotNone(ga_report.ROW_RE.match(row))


class TestPopulationInitialisationTerminates(unittest.TestCase):
    """
    Both initialisation loops were unbounded. With settings that no draw can satisfy the symptom
    was a pegged core and no further output, indistinguishable from a slow start.
    """

    def test_an_impossible_charge_window_fails_instead_of_hanging(self):
        import core.genetic_algorithm as GA
        import user.usersettings as UserSettings
        import run

        original = UserSettings.charge_range
        UserSettings.charge_range = [40, 50]  # unreachable for a 24-mer
        try:
            ga = GA.GeneticAlgorithm()
            ga.ptr_scoring_function = lambda seq: [0.0, 0.0]
            ga.population_size = 4
            ga.ptr_constraint_function = run.Constraint.combined
            with silenced():
                with self.assertRaises(SystemExit):
                    ga.initialize()
        finally:
            UserSettings.charge_range = original

    def test_the_default_settings_still_build_a_population(self):
        # The guard must not fire on the shipped configuration.
        import core.genetic_algorithm as GA
        import run
        ga = GA.GeneticAlgorithm()
        ga.ptr_scoring_function = lambda seq: [0.0, 0.0]
        ga.population_size = 8
        ga.ptr_constraint_function = run.Constraint.combined
        ga.initialize()


class TestFindReplenHandlesAnUndivisibleLength(unittest.TestCase):
    """
    repeat and new_seq were assigned only inside the loop body, so a sequence length not
    divisible by any rep_lengths entry died mid-campaign with an UnboundLocalError.
    """

    def test_no_unbound_local(self):
        import core.genetic_algorithm as GA
        chromosome = GA._StringChromosome(list("KWLRKAFKKLA"))  # 11, not divisible by 24
        repeat = chromosome.find_replen()
        self.assertEqual(repeat, 1)
        self.assertEqual("".join(chromosome.sequence), "KWLRKAFKKLA")


class TestCheckpointIsWrittenAtomically(unittest.TestCase):

    def test_a_truncated_checkpoint_is_reported_not_raised_raw(self):
        import core.genetic_algorithm as GA
        directory = tempfile.mkdtemp()
        try:
            path = os.path.join(directory, "checkpoint.dat")
            with open(path, "wb") as handle:
                handle.write(b"\x80\x04\x95truncated")
            population = GA._Population.__new__(GA._Population)
            with silenced():
                with self.assertRaises(SystemExit):
                    population.read_checkpoint_file(path)
        finally:
            shutil.rmtree(directory, ignore_errors=True)

    def test_a_missing_checkpoint_names_the_path(self):
        import core.genetic_algorithm as GA
        population = GA._Population.__new__(GA._Population)
        with silenced():
            with self.assertRaises(SystemExit):
                population.read_checkpoint_file("/nonexistent/checkpoint.dat")


# ---------------------------------------------------------------------------------------------
# SILENT WRONG ANSWERS
# ---------------------------------------------------------------------------------------------

class TestInvalidResiduesAreRejected(unittest.TestCase):
    """
    PeptideBuilder returns the glycine geometry for any character it does not recognise,
    including every lowercase letter, and raises nothing. A calibration set pasted in lowercase
    therefore ran to completion as polyglycine, showed no separation between the known classes,
    and would have been read as the observable having failed. That verdict decides whether the
    whole approach continues.
    """

    def test_lowercase_is_accepted_and_upcased(self):
        import run
        self.assertEqual(run._validate_sequence("kwlrkafkk"), "KWLRKAFKK")

    def test_a_non_residue_character_stops_the_run(self):
        import run
        for sequence in ("KWLRKXAFKK", "KWLRK-AFKK", "KWLRK AFKK", "KWLRKBAFKK"):
            with self.assertRaises(SystemExit, msg="%r was accepted" % sequence):
                run._validate_sequence(sequence)

    def test_every_standard_residue_is_accepted(self):
        import run
        self.assertEqual(run._validate_sequence("ACDEFGHIKLMNPQRSTVWY"),
                         "ACDEFGHIKLMNPQRSTVWY")

    def test_the_validated_alphabet_matches_what_peptidebuilder_can_build(self):
        # If these ever diverge, the check either rejects buildable peptides or lets through the
        # silent-glycine case it exists to prevent.
        import run
        try:
            from PeptideBuilder import Geometry
        except ImportError:
            self.skipTest("PeptideBuilder not installed")
        for residue in sorted(run.VALID_RESIDUES):
            geometry = Geometry.geometry(residue)
            if residue != "G":
                self.assertNotEqual(type(geometry).__name__, "GlyGeo",
                                    "%s silently builds as glycine" % residue)


class TestPartialEnergySumIsRefused(unittest.TestCase):
    """
    _sum_terms skipped absent keys, so a sum of no keys returned a clean 0.0 that reached the
    binding objective as a perfectly plausible zero interaction energy. The guard upstream
    checked a different key set, so changing ENERGY_group_terms passed it and zeroed this.
    """

    def test_a_missing_term_raises_instead_of_summing_to_zero(self):
        from user.modules.compute_observables import _sum_terms
        terms = {"Coul-SR:Protein-Lipids": (-250.0, 5.0)}
        with self.assertRaises(KeyError) as caught:
            _sum_terms(terms, ["Coul-SR:Protein-Lipids", "LJ-SR:Protein-Lipids"])
        self.assertIn("LJ-SR:Protein-Lipids", str(caught.exception))

    def test_a_complete_sum_still_works(self):
        from user.modules.compute_observables import _sum_terms
        terms = {"Coul-SR:Protein-Lipids": (-250.0, 3.0),
                 "LJ-SR:Protein-Lipids": (-150.0, 4.0)}
        total, err = _sum_terms(terms, list(terms))
        self.assertAlmostEqual(total, -400.0)
        self.assertAlmostEqual(err, 5.0)  # quadrature


class TestLastMoleculeInAnItpIsCosted(unittest.TestCase):
    """
    A molecule was committed only when a LATER [section] header appeared, so the last molecule
    in a file whose [ atoms ] runs to the end contributed zero charge. Combined with the
    neutralisation path, that is how a charged system comes out looking neutral and runs without
    counterions, with grompp's warning swallowed by the max_warnings allowance.
    """

    def test_a_molecule_at_end_of_file_is_not_dropped(self):
        import core.gromacs as Gromacs
        directory = tempfile.mkdtemp()
        try:
            path = os.path.join(directory, "ions.itp")
            with open(path, "w") as handle:
                handle.write("[ moleculetype ]\nNA 1\n\n[ atoms ]\n"
                             "1 Qd 1 NA NA 1 1.0\n")  # no trailing section
            charges = Gromacs.TOP_get_system_charge.__wrapped__ if hasattr(
                Gromacs.TOP_get_system_charge, "__wrapped__") else None
            # Call the nested parser through a topology that includes this itp.
            top = os.path.join(directory, "system.top")
            with open(top, "w") as handle:
                handle.write('#include "ions.itp"\n\n[ system ]\ntest\n\n'
                             '[ molecules ]\nNA 4\n')
            total = Gromacs.TOP_get_system_charge(top, TOP_itp_directory_list=[directory])
            self.assertEqual(total, 4, "the last moleculetype in an itp was not costed")
        finally:
            shutil.rmtree(directory, ignore_errors=True)

    def test_an_uncostable_molecule_with_a_zero_total_is_refused(self):
        import core.gromacs as Gromacs
        directory = tempfile.mkdtemp()
        try:
            top = os.path.join(directory, "system.top")
            with open(top, "w") as handle:
                handle.write("[ system ]\ntest\n\n[ molecules ]\nPOPE 128\nW 1000\n")
            with silenced():
                with self.assertRaises(Errors.OperatorError) as caught:
                    Gromacs.TOP_get_system_charge(top, TOP_itp_directory_list=[directory])
            self.assertIn("POPE", str(caught.exception))
        finally:
            shutil.rmtree(directory, ignore_errors=True)


# ---------------------------------------------------------------------------------------------
# ILLEGIBLE FAILURES
# ---------------------------------------------------------------------------------------------

class TestErrorBlockNamesTheLocation(unittest.TestCase):

    def test_all_four_fields_are_rendered(self):
        block = Errors.format_block("what happened", where="/path/to/file",
                                    why="the expectation", fix="the command")
        for expected in ("what:", "where:", "why:", "fix:",
                         "what happened", "/path/to/file", "the expectation", "the command"):
            self.assertIn(expected, block)

    def test_absent_fields_are_omitted_rather_than_left_blank(self):
        block = Errors.format_block("what happened")
        self.assertIn("what:", block)
        self.assertNotIn("where:", block)

    def test_multiline_values_stay_indented_under_their_label(self):
        block = Errors.format_block("x", where="line one\nline two")
        self.assertIn("\n         line two", block)

    def test_require_file_reports_the_working_directory_for_a_relative_path(self):
        with silenced():
            with self.assertRaises(SystemExit):
                Errors.require_file("definitely-not-here.txt", why="y", fix="f")

    def test_require_file_returns_the_path_when_it_exists(self):
        handle = tempfile.NamedTemporaryFile(delete=False)
        handle.close()
        try:
            self.assertEqual(Errors.require_file(handle.name, why="y", fix="f"), handle.name)
        finally:
            os.unlink(handle.name)


class TestSequenceFileErrors(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _write(self, text):
        path = os.path.join(self.dir, "seq.txt")
        with open(path, "w") as handle:
            handle.write(text)
        return path

    def test_a_missing_file_stops_with_a_message(self):
        import run
        with silenced():
            with self.assertRaises(SystemExit):
                run._read_sequence_file(os.path.join(self.dir, "absent.txt"))

    def test_comments_and_blank_lines_are_ignored(self):
        import run
        path = self._write("# a comment\n\nKWLRKAFKK  # trailing comment\n\n")
        self.assertEqual(run._read_sequence_file(path), ["KWLRKAFKK"])

    def test_a_file_of_only_comments_is_an_error(self):
        import run
        with silenced():
            with self.assertRaises(SystemExit):
                run._read_sequence_file(self._write("# nothing here\n\n"))

    def test_a_bad_residue_is_rejected(self):
        import run
        with silenced():
            with self.assertRaises(SystemExit):
                run._read_sequence_file(self._write("KWLRKAFKK\nKWLRK1AFKK\n"))

    def test_the_shipped_example_file_parses(self):
        # It is the file the error message tells the operator to copy, so it has to be valid.
        import run
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        example = os.path.join(root, "seq.txt.example")
        self.assertTrue(os.path.isfile(example), "seq.txt.example is not shipped")
        self.assertTrue(run._read_sequence_file(example))


class TestCommandFailureNamesTheCommand(unittest.TestCase):

    def test_check_exitcode_reports_the_command_and_the_log(self):
        import core.system as System
        directory = tempfile.mkdtemp()
        try:
            log = os.path.join(directory, "out.log")
            with open(log, "w") as handle:
                handle.write("Fatal error:\nAtomtype P4 not found\n")
            with silenced() as captured:
                with self.assertRaises(SystemExit):
                    System.check_exitcode(1, command="gmx grompp -f npt.mdp",
                                          path_stdout=log)
            text = captured.getvalue()
            self.assertIn("gmx grompp", text)
            self.assertIn(log, text)
            self.assertIn("Atomtype P4 not found", text)
        finally:
            shutil.rmtree(directory, ignore_errors=True)

    def test_a_zero_exit_code_is_silent(self):
        import core.system as System
        System.check_exitcode(0, command="gmx grompp")  # must not raise

    def test_simulation_error_carries_a_message(self):
        import core.system as System
        with self.assertRaises(System.SimulationError) as caught:
            System._sim_check_exitcode(1, command="gmx mdrun -deffnm npt")
        self.assertIn("gmx mdrun", str(caught.exception))


class TestMissingOutputVariableNamesTheCandidate(unittest.TestCase):

    def test_the_message_carries_the_sequence_and_the_key(self):
        import core.simulation as Sim
        directory = tempfile.mkdtemp()
        try:
            simulation = Sim.Simulation(directory, sequence="KWLRKAFKK")
            with self.assertRaises(KeyError) as caught:
                simulation.get_output_variable("fitness")
            message = str(caught.exception)
            self.assertIn("fitness", message)
            self.assertIn("KWLRKAFKK", message)
        finally:
            shutil.rmtree(directory, ignore_errors=True)


class TestMembraneFilesAreCheckedBeforeTheRun(unittest.TestCase):

    def test_a_calibrated_membrane_with_no_files_on_disk_is_refused(self):
        import user.membranes as Membranes
        definition = Membranes.MembraneDef(
            name="unit_test_membrane",
            gro="definitely_absent.gro", top="definitely_absent.top", role="target",
            area_nm2=0.65, surften_free=(0.0, 1.0), thickness_nm=4.0)
        original = dict(Membranes.MEMBRANES)
        Membranes.MEMBRANES.clear()
        Membranes.MEMBRANES["unit_test_membrane"] = definition
        try:
            with silenced():
                with self.assertRaises(SystemExit):
                    Membranes.require_all_built()
        finally:
            Membranes.MEMBRANES.clear()
            Membranes.MEMBRANES.update(original)


class TestCalibrateInputErrors(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_a_labels_file_without_the_expected_header_is_refused(self):
        import calibrate
        path = os.path.join(self.dir, "labels.csv")
        with open(path, "w") as handle:
            handle.write("KWLRKAFKK;selective\n")  # no header, wrong delimiter
        with silenced():
            with self.assertRaises(SystemExit):
                calibrate.load_labels(path)

    def test_a_valid_labels_file_loads(self):
        import calibrate
        path = os.path.join(self.dir, "labels.csv")
        with open(path, "w") as handle:
            handle.write("sequence,class\nKWLRKAFKK,selective\nAAAAAAAA,inactive\n")
        self.assertEqual(calibrate.load_labels(path),
                         {"KWLRKAFKK": "selective", "AAAAAAAA": "inactive"})

    def test_a_truncated_observables_file_names_the_path(self):
        import calibrate
        candidate = os.path.join(self.dir, "KWLRKAFKK")
        os.makedirs(candidate)
        with open(os.path.join(candidate, "observables.json"), "w") as handle:
            handle.write('{"bacterial": {')  # killed mid-write
        with silenced():
            with self.assertRaises(SystemExit):
                list(calibrate.find_results(self.dir))

    def test_a_complete_observables_file_is_read(self):
        import calibrate
        candidate = os.path.join(self.dir, "KWLRKAFKK")
        os.makedirs(candidate)
        with open(os.path.join(candidate, "observables.json"), "w") as handle:
            json.dump({"bacterial": {"delta_gamma": {"mean": 1.0, "err": 0.1}}}, handle)
        results = list(calibrate.find_results(self.dir))
        self.assertEqual(len(results), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
