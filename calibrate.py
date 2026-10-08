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

import core.errors as Errors
import user.usersettings as UserSettings
from user.modules.compute_fitness import score_observables


# The expected fitness order of the classes depends on the objective, and getting this wrong
# produces a confident, wrong go/no-go verdict on a working configuration.
#
# damage  maximises strain on the target membrane, so a lytic peptide should outrank a selective
#         one, which should outrank an inactive one.
# binding maximises affinity for the target MINUS affinity for the host. A non-selective lytic
#         peptide binds the host strongly and is penalised for it, so it belongs between the
#         inactive and the selective peptides, not above both.
CLASS_ORDER_BY_OBJECTIVE = {
    "damage":  {"inactive": 0, "selective": 1, "lytic": 2},
    "binding": {"inactive": 0, "lytic": 1, "selective": 2},
}
CLASS_NAMES = ("inactive", "selective", "lytic")


def find_results(root):
    """Yield (sequence, observables, dirpath) for every observables.json under root."""
    for dirpath, _dirnames, filenames in os.walk(root):
        if "observables.json" not in filenames:
            continue
        path = os.path.join(dirpath, "observables.json")
        try:
            with open(path) as handle:
                observables = json.load(handle)
        except ValueError as exc:
            # A candidate killed mid-write leaves a partial file. Walking hundreds of
            # directories, the bare JSONDecodeError names a line and column of an unnamed file.
            Errors.fail("observables.json is not valid JSON (%s)" % exc,
                        where=os.path.abspath(path),
                        why=("it was most likely truncated by a kill while it was being "
                             "written, so this candidate has no usable measurement"),
                        fix=("delete this candidate directory and re-run that sequence, or "
                             "exclude it and re-run calibrate.py"))
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
    Errors.require_file(
        path,
        why="the labels say which experimental class each calibration peptide belongs to",
        fix="write a CSV with a header line: sequence,class")

    labels = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)

        # Checked once, before the rows, because the failure mode otherwise is a KeyError on
        # the first row that names the missing column but not the file, and the three usual
        # causes (no header line, different column names, semicolon delimiter) all look the same
        # from there.
        fields = [name.strip().lower() for name in (reader.fieldnames or [])]
        missing = [name for name in ("sequence", "class") if name not in fields]
        if missing:
            Errors.fail("the labels file is missing the column(s): %s" % ", ".join(missing),
                        where="%s\nheader read as: %s"
                              % (os.path.abspath(path),
                                 ", ".join(reader.fieldnames or []) or "(no header line)"),
                        why=("it is read as a CSV with a header naming at least 'sequence' and "
                             "'class'. A file with no header, or one separated by semicolons "
                             "or tabs, reads as a single unnamed column."),
                        fix="first line: sequence,class")

        for number, row in enumerate(reader, start=2):
            seq = (row.get("sequence") or "").strip().upper()
            cls = (row.get("class") or "").strip().lower()
            if not seq:
                continue
            if cls not in CLASS_NAMES:
                Errors.fail("unknown class %r for %s" % (cls, seq),
                            where="%s, line %d" % (os.path.abspath(path), number),
                            why="the class decides which side of a comparison a peptide is on",
                            fix="use one of: %s" % ", ".join(CLASS_NAMES))
            labels[seq] = cls

    if not labels:
        Errors.fail("no labels were read from %s" % os.path.basename(path),
                    where=os.path.abspath(path),
                    why="the file has a usable header but no data rows",
                    fix="add one row per labelled peptide: SEQUENCE,class")
    return labels


def concordance(rows, objective):
    """
    Fraction of labelled pairs ordered correctly by fitness, for this objective.

    A plain rank-agreement statistic over the pairs whose true classes differ. Pairs within a
    class are skipped, since the labels carry no ordering inside a class. This is a coarse check
    and is meant to be: at this stage the question is whether the signal exists at all, not how
    strong it is. The expected class order comes from the objective, because the two objectives
    do not rank a non-selective lytic peptide the same way.
    """
    order = CLASS_ORDER_BY_OBJECTIVE[objective]
    labelled = [r for r in rows if r["class"] is not None]
    concordant = 0
    total = 0
    for i in range(len(labelled)):
        for j in range(i + 1, len(labelled)):
            a, b = labelled[i], labelled[j]
            ra, rb = order[a["class"]], order[b["class"]]
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
            "gmx_version": observables.get("_provenance", {}).get("gmx_version"),
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

    # A set of results produced by more than one GROMACS build is not internally comparable,
    # and the shift it introduces looks exactly like a real effect of the sequences.
    builds = {r["gmx_version"] for r in rows if r["gmx_version"]}
    if len(builds) > 1:
        print("WARNING: these results were produced by more than one GROMACS build:")
        for b in sorted(builds):
            n = sum(1 for r in rows if r["gmx_version"] == b)
            print("  %-28s %d candidate(s)" % (b, n))
        print("  Surface tension is build-sensitive, so the values are not directly")
        print("  comparable across them. Re-run the minority set under the pinned build.")
        print()

    # How many of the supplied labels actually found a result, and the reverse. A label that
    # matches nothing is silently dropped, and the concordance below is then computed over
    # whatever subset happened to line up, printed with the same confidence as a full one. The
    # usual cause is that the calibration run names its directories SEQUENCE + rank + "_" +
    # count, so a candidate whose input.seq is missing reports a "sequence" like KWLRK0_3 that
    # can never match a label.
    if labels:
        matched = {r["sequence"].upper() for r in rows if r["class"]}
        unmatched_labels = sorted(set(labels) - matched)
        unlabelled_rows = sorted(r["sequence"] for r in rows if not r["class"])
        print("labels    = %d supplied, %d matched a result" % (len(labels), len(matched)))
        if unmatched_labels:
            print("  %d label(s) matched nothing: %s%s"
                  % (len(unmatched_labels), ", ".join(unmatched_labels[:5]),
                     "" if len(unmatched_labels) <= 5 else ", ..."))
        if unlabelled_rows:
            print("  %d result(s) have no label: %s%s"
                  % (len(unlabelled_rows), ", ".join(s[:26] for s in unlabelled_rows[:5]),
                     "" if len(unlabelled_rows) <= 5 else ", ..."))
        if unmatched_labels or unlabelled_rows:
            print("  Everything below is computed on the %d matched peptides only." % len(matched))

    objective = args.objective or UserSettings.OBJECTIVE
    print("objective = %s" % objective)
    if builds:
        print("gmx       = %s" % ", ".join(sorted(builds)))
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
        score, pairs = concordance(rows, objective)
        print()
        order = CLASS_ORDER_BY_OBJECTIVE[objective]
        print("expected order for '%s': %s"
              % (objective, " < ".join(sorted(order, key=order.get))))
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
