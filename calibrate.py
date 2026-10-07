"""
Score an existing set of observables.json files under either objective, and report whether
the observables separate the calibration peptides the way the literature does.

This is deliberately a separate script from run.py and it runs no simulations. Two uses:

  1. After the calibration set has been measured (run.py with STATE = 1), check whether the
     measured observables rank the known peptides correctly. That is the go/no-go for the
     whole approach: if delta_gamma cannot tell a lytic peptide from a selective one from an
     inactive one, no amount of optimiser sophistication downstream will fix it.

  2. After any campaign, re-score it under the other objective. Because both objectives read
     the same observables, comparing the binding ranking with the damage ranking is a
     post-processing step, not a second campaign. Whether the two rankings agree is itself a
     result: agreement means the two observables are redundant in this system, disagreement
     is the interesting case and is the thing worth writing about.

Usage:
    python calibrate.py <results_dir> [--objective damage|binding] [--labels labels.csv]

<results_dir> is searched recursively for observables.json files. The sequence is taken from
the input.seq file next to each one.

labels.csv, if given, has a header and two columns: sequence,class
where class is one of inactive, selective, lytic. Supplying it turns the output from a bare
ranking into a check against what is already known.
"""

import argparse
import csv
import json
import os
import sys

import user.usersettings as UserSettings
from user.modules.compute_fitness import score_observables


CLASS_ORDER = {"inactive": 0, "selective": 1, "lytic": 2}


def find_results(root):
    """Yield (sequence, observables) for every observables.json under root."""
    for dirpath, _dirnames, filenames in os.walk(root):
        if "observables.json" not in filenames:
            continue
        with open(os.path.join(dirpath, "observables.json")) as handle:
            observables = json.load(handle)
        sequence = read_sequence(dirpath)
        yield sequence, observables, dirpath


def read_sequence(dirpath):
    """Read the FASTA written by _setup_simulation_directory, falling back to the dir name."""
    seq_path = os.path.join(dirpath, UserSettings.filename_seq)
    if os.path.exists(seq_path):
        lines = [line.strip() for line in open(seq_path) if line.strip()]
        body = [line for line in lines if not line.startswith(">")]
        if body:
            return "".join(body)
    return os.path.basename(dirpath)


def load_labels(path):
    labels = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            seq = row["sequence"].strip().upper()
            cls = row["class"].strip().lower()
            if cls not in CLASS_ORDER:
                sys.exit("Unknown class %r for %s; expected one of %s"
                         % (cls, seq, ", ".join(CLASS_ORDER)))
            labels[seq] = cls
    return labels


def concordance(rows):
    """
    Fraction of labelled pairs ordered correctly by fitness.

    A plain rank-agreement statistic over the pairs whose true classes differ. Pairs within a
    class are skipped, since the labels carry no ordering inside a class. This is a coarse
    check and is meant to be: at this stage the question is whether the signal exists at all,
    not how strong it is.
    """
    labelled = [r for r in rows if r["class"] is not None]
    concordant = 0
    total = 0
    for i in range(len(labelled)):
        for j in range(i + 1, len(labelled)):
            a, b = labelled[i], labelled[j]
            ra, rb = CLASS_ORDER[a["class"]], CLASS_ORDER[b["class"]]
            if ra == rb:
                continue
            total += 1
            if (a["fitness"] - b["fitness"]) * (ra - rb) > 0:
                concordant += 1
    if total == 0:
        return None, 0
    return concordant / float(total), total


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results_dir")
    parser.add_argument("--objective", choices=("damage", "binding"), default=None,
                        help="defaults to UserSettings.OBJECTIVE")
    parser.add_argument("--labels", default=None,
                        help="CSV with columns sequence,class (inactive|selective|lytic)")
    args = parser.parse_args()

    labels = load_labels(args.labels) if args.labels else {}

    rows = []
    failures = []
    for sequence, observables, dirpath in find_results(args.results_dir):
        try:
            fitness, components = score_observables(observables, objective=args.objective)
        except Exception as exc:
            # A proxy that is unavailable for this dataset is a configuration result, not a
            # crash: report it per candidate and carry on.
            failures.append((sequence, dirpath, str(exc)))
            continue
        rows.append({
            "sequence": sequence,
            "fitness": fitness,
            "components": components,
            "class": labels.get(sequence.upper()),
        })

    if failures:
        print("Could not score %d candidate(s):" % len(failures))
        for sequence, dirpath, reason in failures[:5]:
            print("  %-26s %s" % (sequence[:26], reason))
        if len(failures) > 5:
            print("  ... and %d more" % (len(failures) - 5))
        print()

    if not rows:
        sys.exit("No scorable observables.json found under %s" % args.results_dir)

    rows.sort(key=lambda r: r["fitness"], reverse=True)

    objective = args.objective or UserSettings.OBJECTIVE
    print("objective = %s" % objective)
    if objective == "binding":
        print("proxy     = %s" % UserSettings.BINDING_proxy)
    print("n         = %d" % len(rows))
    print()

    header = "%-28s %10s %10s %7s  %s" % ("sequence", "fitness", "raw", "orient", "class")
    print(header)
    print("-" * len(header))
    for row in rows:
        print("%-28s %10.3f %10.3f %7.3f  %s" % (
            row["sequence"][:28],
            row["fitness"],
            row["components"]["raw"],
            row["components"]["orientation_factor"],
            row["class"] or "",
        ))

    if labels:
        score, pairs = concordance(rows)
        print()
        if score is None:
            print("No cross-class pairs to compare; add peptides from more than one class.")
        else:
            print("Pairwise concordance with the labels: %.2f over %d pairs" % (score, pairs))
            print()
            # The interpretation is the point of this script, so it is stated rather than left
            # to the reader of a bare number.
            if score >= 0.8:
                print("The observable orders the calibration set. Proceed to set the window "
                      "boundaries from where the classes change over.")
            elif score >= 0.6:
                print("Weak separation. Before changing the objective, check the obvious "
                      "confounders: averaging window length, whether any candidate has a low "
                      "orientation factor (peptide left the surface), and whether the host "
                      "and target membranes were equilibrated to the same standard.")
            else:
                print("No usable separation. The optimiser cannot recover a signal that is "
                      "not in the observable; revisit the observable, the membranes or the "
                      "imposed-helix assumption before running a GA campaign.")


if __name__ == "__main__":
    main()
