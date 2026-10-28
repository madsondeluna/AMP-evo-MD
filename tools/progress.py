"""
Progress of a STATE = 1 calibration or validation run: which (sequence, replica) tasks have
finished and which are still pending, as a per-sequence grid with a total.

A task is counted:
  done  it has written observables.json (the file calibrate.py reads)
  run   its directory exists but no observables.json yet (in progress, or failed mid-write)
  wait  the directory has not been created yet (still queued)

This is pure filesystem inspection: it starts no simulation and does not touch a running job, so
it is safe to call at any time, including under watch(1).

Usage:
    python tools/progress.py [results_dir] [--seq seq.txt] [--nrep 3]

Defaults: results_dir = validation_run, seq = seq.txt, nrep = 3. Paths resolve from the current
working directory, so run it from the repository root (or pass absolute paths).
"""

import argparse
import glob
import os


def read_sequences(seq_path):
    # One sequence per line, one-letter codes; blank lines and anything after a # are ignored.
    seqs = []
    with open(seq_path) as handle:
        for line in handle:
            s = line.split("#", 1)[0].strip()
            if s:
                seqs.append(s.upper())
    return seqs


def task_status(results_dir, seq, replica):
    # Task directory is named <SEQUENCE>_rep<r> (see run.py _run_specific_sequences).
    task_dir = os.path.join(results_dir, "%s_rep%d" % (seq, replica))
    if not os.path.isdir(task_dir):
        return "wait"
    # observables.json may sit directly in the task dir or under its retained work/ copy;
    # the recursive glob matches both.
    if glob.glob(os.path.join(task_dir, "**", "observables.json"), recursive=True):
        return "done"
    return "run"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results_dir", nargs="?", default="validation_run")
    parser.add_argument("--seq", default="seq.txt")
    parser.add_argument("--nrep", type=int, default=3)
    args = parser.parse_args()

    if not os.path.exists(args.seq):
        parser.error("sequence file not found: %s" % os.path.abspath(args.seq))

    seqs = read_sequences(args.seq)
    reps = range(1, args.nrep + 1)

    # Header: one column per replica.
    cols = "".join("%6s" % ("rep%d" % r) for r in reps)
    header = "%-28s%s" % ("sequence", cols)
    print(header)
    print("-" * len(header))

    counts = {"done": 0, "run": 0, "wait": 0}
    for seq in seqs:
        row = [task_status(args.results_dir, seq, r) for r in reps]
        for st in row:
            counts[st] += 1
        print("%-28s%s" % (seq[:28], "".join("%6s" % st for st in row)))

    total = len(seqs) * args.nrep
    print("-" * len(header))
    print("done %d / %d    run %d    wait %d"
          % (counts["done"], total, counts["run"], counts["wait"]))


if __name__ == "__main__":
    main()
