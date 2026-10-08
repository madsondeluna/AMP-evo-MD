"""
Tests for per-step logging.

What this exists to guarantee, in order of importance:

A failure must be attributable to a step. The pipeline was previously silent about which of its
five stages it was in, so a candidate that hung or died gave no indication of where, and the only
artefact was one concatenated stdout file with no separators.

Logging must never break the run it observes. A full disk or a read-only directory should cost
the log and nothing else, which is why every write is guarded and the tests check that a
failure to log does not propagate.

The format must survive a kill. JSON Lines is chosen so that every completed record is valid on
its own; a line truncated mid-write is skipped and everything before it still reads.
"""

import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import core.simulation as Sim
import core.steplog as steplog


def _noop(simulation):
    pass


def _raise_system_exit(simulation):
    # How core.system.check_exitcode reports a non-zero GROMACS command.
    sys.exit(3)


def _raise_runtime(simulation):
    raise RuntimeError("something went wrong")


class TestStepLogFormat(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_each_step_produces_one_record(self):
        log = steplog.StepLog(self.dir, "KWLRK")
        for name in ("a", "b", "c"):
            with log.step(name):
                pass
        records = [r for r in steplog.read(os.path.join(self.dir, steplog.FILENAME))
                   if r.get("event") == "step"]
        self.assertEqual([r["step"] for r in records], ["a", "b", "c"])
        self.assertTrue(all(r["status"] == "ok" for r in records))

    def test_records_carry_the_sequence(self):
        log = steplog.StepLog(self.dir, "KWLRKAFKK")
        with log.step("a"):
            pass
        records = steplog.read(os.path.join(self.dir, steplog.FILENAME))
        self.assertEqual(records[0]["sequence"], "KWLRKAFKK")

    def test_failure_is_recorded_and_re_raised(self):
        log = steplog.StepLog(self.dir, "AAAA")
        with self.assertRaises(RuntimeError):
            with log.step("insert"):
                raise RuntimeError("boom")
        record = steplog.read(os.path.join(self.dir, steplog.FILENAME))[0]
        self.assertEqual(record["status"], "failed")
        self.assertEqual(record["error_type"], "RuntimeError")
        self.assertIn("boom", record["error"])

    def test_systemexit_is_recorded(self):
        # The case that matters: a non-zero GROMACS command.
        log = steplog.StepLog(self.dir, "AAAA")
        with self.assertRaises(SystemExit):
            with log.step("production"):
                sys.exit(1)
        record = steplog.read(os.path.join(self.dir, steplog.FILENAME))[0]
        self.assertEqual(record["status"], "failed")
        self.assertEqual(record["error_type"], "SystemExit")

    def test_a_truncated_line_does_not_lose_the_rest(self):
        path = os.path.join(self.dir, steplog.FILENAME)
        with open(path, "w") as handle:
            handle.write(json.dumps({"event": "step", "step": "a", "status": "ok"}) + "\n")
            handle.write('{"event": "step", "step": "b"')  # killed mid-write
        records = steplog.read(path)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["step"], "a")

    def test_unwritable_directory_does_not_raise(self):
        # Logging must never be the reason a campaign stops.
        log = steplog.StepLog("/proc/nonexistent/cannot/write", "AAAA")
        with log.step("a"):
            pass
        log.event("finish")


class TestSummarise(unittest.TestCase):

    def test_totals_and_share(self):
        records = [
            {"event": "step", "step": "production", "status": "ok", "seconds": 100.0},
            {"event": "step", "step": "production", "status": "ok", "seconds": 300.0},
            {"event": "step", "step": "insert", "status": "ok", "seconds": 50.0},
            {"event": "step", "step": "insert", "status": "failed", "seconds": 50.0},
        ]
        rows = steplog.summarise(records)
        by_name = {r["step"]: r for r in rows}
        self.assertEqual(by_name["production"]["n"], 2)
        self.assertEqual(by_name["production"]["total_seconds"], 400.0)
        self.assertEqual(by_name["insert"]["failed"], 1)
        self.assertAlmostEqual(sum(r["share"] for r in rows), 1.0)

    def test_sorted_by_total_time(self):
        # The slowest step must come first: it is the only one worth optimising.
        records = [
            {"event": "step", "step": "cheap", "status": "ok", "seconds": 1.0},
            {"event": "step", "step": "expensive", "status": "ok", "seconds": 1000.0},
        ]
        self.assertEqual(steplog.summarise(records)[0]["step"], "expensive")

    def test_non_step_events_are_ignored(self):
        records = [
            {"event": "start", "modules": ["a"]},
            {"event": "step", "step": "a", "status": "ok", "seconds": 1.0},
            {"event": "finish", "status": "ok"},
        ]
        rows = steplog.summarise(records)
        self.assertEqual(len(rows), 1)

    def test_empty_input(self):
        self.assertEqual(steplog.summarise([]), [])

    def test_duration_formatting_picks_a_sensible_unit(self):
        self.assertIn("s", steplog.format_duration(12))
        self.assertIn("min", steplog.format_duration(600))
        self.assertIn("h", steplog.format_duration(7200))
        self.assertIn("d", steplog.format_duration(400000))


class TestSimulationIntegration(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _records(self):
        return steplog.read(os.path.join(self.dir, steplog.FILENAME))

    def test_successful_run_logs_every_step(self):
        sim = Sim.Simulation(self.dir, sequence="KWLRK")
        for name in ("build_peptide", "insert", "production"):
            sim.add_module(Sim.Module(_noop, name))
        sim.run()
        steps = [r["step"] for r in self._records() if r.get("event") == "step"]
        self.assertEqual(steps, ["build_peptide", "insert", "production"])
        self.assertTrue(any(r.get("event") == "finish" for r in self._records()))

    def test_module_name_defaults_to_the_function_name(self):
        sim = Sim.Simulation(self.dir)
        sim.add_module(Sim.Module(_noop))
        sim.run()
        steps = [r["step"] for r in self._records() if r.get("event") == "step"]
        self.assertEqual(steps, ["_noop"])

    def test_failure_names_the_step_and_stops_there(self):
        sim = Sim.Simulation(self.dir, sequence="AAAA")
        sim.add_module(Sim.Module(_noop, "build_peptide"))
        sim.add_module(Sim.Module(_raise_system_exit, "insert"))
        sim.add_module(Sim.Module(_noop, "production"))
        with self.assertRaises(SystemExit):
            sim.run()

        records = self._records()
        steps = {r["step"]: r for r in records if r.get("event") == "step"}
        self.assertEqual(steps["build_peptide"]["status"], "ok")
        self.assertEqual(steps["insert"]["status"], "failed")
        # The step after the failure must not have run.
        self.assertNotIn("production", steps)

        aborts = [r for r in records if r.get("event") == "abort"]
        self.assertEqual(len(aborts), 1)
        self.assertEqual(aborts[0]["step"], "insert")

    def test_intermediates_are_preserved_on_failure(self):
        # TemporaryDirectory deletes on exit, which removes the evidence at exactly the moment
        # it is needed.
        def writes_then_fails(simulation):
            with open(os.path.join(simulation._temp_dir, "intermediate.txt"), "w") as handle:
                handle.write("state at the moment of failure")
            raise RuntimeError("boom")

        sim = Sim.Simulation(self.dir, sequence="AAAA")
        sim.add_module(Sim.Module(writes_then_fails, "insert"))
        with self.assertRaises(RuntimeError):
            sim.run()

        preserved = os.path.join(self.dir, "failed_insert", "intermediate.txt")
        self.assertTrue(os.path.exists(preserved), "intermediates were not preserved")
        with open(preserved) as handle:
            self.assertIn("moment of failure", handle.read())

    def test_preservation_can_be_switched_off(self):
        sim = Sim.Simulation(self.dir, sequence="AAAA", keep_on_failure=False)
        sim.add_module(Sim.Module(_raise_runtime, "insert"))
        with self.assertRaises(RuntimeError):
            sim.run()
        self.assertFalse(os.path.exists(os.path.join(self.dir, "failed_insert")))

    def test_collect_walks_a_campaign_directory(self):
        campaign = tempfile.mkdtemp()
        try:
            for name in ("CAND1", "CAND2"):
                d = os.path.join(campaign, name)
                os.makedirs(d)
                sim = Sim.Simulation(d, sequence=name)
                sim.add_module(Sim.Module(_noop, "production"))
                sim.run()
            records = steplog.collect(campaign)
            steps = [r for r in records if r.get("event") == "step"]
            self.assertEqual(len(steps), 2)
            self.assertTrue(all("_dir" in r for r in records))
        finally:
            shutil.rmtree(campaign, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
