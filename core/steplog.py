"""
Per-step logging for a candidate evaluation.

The problem this solves: a campaign runs for days and the pipeline was silent about which of
its five stages it was in. A candidate that took an hour gave no way to tell whether the time
went to martinize2, to grompp, or to one of the two 400 ns production runs, and a candidate that
failed left only a single concatenated stdout file with no marker for where one command ended
and the next began.

Format is JSON Lines, one object per line, appended as each step finishes. That choice is
deliberate: the file is readable after a crash mid-campaign because every completed line is
already flushed and valid on its own, which a single JSON array would not be. It also appends
without reading, so a long campaign never rewrites what it has already written.

Each candidate directory gets its own steps.jsonl, so MPI ranks never contend for the same file.

A record looks like:

    {"step": "production", "index": 3, "status": "ok", "started": "...", "seconds": 4821.3,
     "sequence": "KWLRK...", "files": ["production_bacterial.edr", ...]}

and on failure carries "status": "failed" with the exception type and message.
"""

import datetime
import json
import os
import time
import traceback


FILENAME = "steps.jsonl"


def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


class StepLog:
    """
    Append-only step record for one candidate.

    Never raises. Logging that can break the run it is observing is worse than no logging, so
    every write is guarded: a full disk or a read-only directory costs the log, not the campaign.
    """

    def __init__(self, directory, sequence=None):
        self.directory = directory
        self.sequence = sequence
        self.path = os.path.join(directory, FILENAME) if directory else None
        self._failed = False

    def _write(self, record):
        if self.path is None or self._failed:
            return
        try:
            os.makedirs(self.directory, exist_ok=True)
            with open(self.path, "a") as handle:
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                # Flush and fsync so a record survives a scheduler kill. The cost is a few
                # milliseconds per step against steps that run for minutes to hours.
                handle.flush()
                os.fsync(handle.fileno())
        except OSError:
            # Stop trying after the first failure rather than attempting a write per step on a
            # directory that is not writable.
            self._failed = True

    def event(self, kind, **fields):
        record = {"event": kind, "time": _now()}
        if self.sequence:
            record["sequence"] = self.sequence
        record.update(fields)
        self._write(record)

    def step(self, name, index=None):
        """Context manager timing one step and recording how it ended."""
        return _Step(self, name, index)


class _Step:
    def __init__(self, log, name, index):
        self.log = log
        self.name = name
        self.index = index
        self.started = None
        self.t0 = None

    def __enter__(self):
        self.started = _now()
        self.t0 = time.time()
        return self

    def __exit__(self, exc_type, exc, tb):
        record = {
            "event": "step",
            "step": self.name,
            "started": self.started,
            "seconds": round(time.time() - self.t0, 2),
            "time": _now(),
        }
        if self.index is not None:
            record["index"] = self.index
        if self.log.sequence:
            record["sequence"] = self.log.sequence

        if exc_type is None:
            record["status"] = "ok"
        else:
            record["status"] = "failed"
            record["error_type"] = exc_type.__name__
            record["error"] = str(exc)[:2000]
            # The last few frames are what identify the failing call; the whole stack would
            # bloat a file that one line per step is meant to keep scannable.
            record["traceback"] = "".join(traceback.format_exception(exc_type, exc, tb))[-2000:]

        self.log._write(record)
        # Never suppress: the caller's error handling decides what a failure means.
        return False


def read(path):
    """
    Parse a steps.jsonl, skipping anything unparseable.

    A record truncated by a kill mid-write is expected and is not an error; the point of the
    format is that everything before it is still usable.
    """
    records = []
    try:
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    records.append(json.loads(line))
                except ValueError:
                    continue
    except OSError:
        return []
    return records


def collect(root, filename=FILENAME):
    """Every step record under a campaign directory, with the candidate directory attached."""
    out = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if filename not in filenames:
            continue
        for record in read(os.path.join(dirpath, filename)):
            record["_dir"] = dirpath
            out.append(record)
    return out


def summarise(records):
    """
    Per-step totals: count, failures, and the time distribution.

    Reports the median and the 90th percentile rather than only the mean, because what matters
    operationally is the tail: one candidate in fifty taking ten times as long is what decides
    whether a campaign finishes this week, and a mean hides it.
    """
    import collections
    import statistics

    by_step = collections.OrderedDict()
    for r in records:
        if r.get("event") != "step":
            continue
        name = r.get("step", "?")
        by_step.setdefault(name, {"n": 0, "failed": 0, "seconds": []})
        entry = by_step[name]
        entry["n"] += 1
        if r.get("status") == "failed":
            entry["failed"] += 1
        seconds = r.get("seconds")
        if isinstance(seconds, (int, float)):
            entry["seconds"].append(float(seconds))

    rows = []
    for name, entry in by_step.items():
        times = sorted(entry["seconds"])
        row = {
            "step": name,
            "n": entry["n"],
            "failed": entry["failed"],
            "total_seconds": sum(times),
            "mean_seconds": statistics.fmean(times) if times else 0.0,
            "median_seconds": statistics.median(times) if times else 0.0,
            "p90_seconds": times[min(len(times) - 1, int(0.9 * len(times)))] if times else 0.0,
            "max_seconds": max(times) if times else 0.0,
        }
        rows.append(row)

    grand_total = sum(r["total_seconds"] for r in rows) or 1.0
    for row in rows:
        row["share"] = row["total_seconds"] / grand_total
    rows.sort(key=lambda r: -r["total_seconds"])
    return rows


def format_duration(seconds):
    """Human-readable duration, chosen by magnitude rather than a fixed unit."""
    seconds = float(seconds)
    if seconds < 90:
        return "%.1f s" % seconds
    if seconds < 5400:
        return "%.1f min" % (seconds / 60.0)
    if seconds < 172800:
        return "%.1f h" % (seconds / 3600.0)
    return "%.1f d" % (seconds / 86400.0)
