"""
Generate a synthetic GA progress file in the exact format core/genetic_algorithm.py writes.

Used to exercise ga_report.py end to end without running simulations. The simulated campaign
converges: fitness rises and then plateaus while diversity falls, which is the behaviour the
report is meant to detect and describe, so a run over this file is also a check that the
written interpretation fires correctly.
"""

import argparse
import random
from datetime import datetime


ALPHABET = list("AEFKLMQRSWY")   # the 11 residues enabled in user/usersettings.py


def make_population(rng, size, length, seed_motif=None, conservation=0.0):
    pop = []
    for _ in range(size):
        seq = []
        for pos in range(length):
            if seed_motif and rng.random() < conservation:
                seq.append(seed_motif[pos % len(seed_motif)])
            else:
                seq.append(rng.choice(ALPHABET))
        pop.append("".join(seq))
    return pop


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="fake_progress.txt")
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--population", type=int, default=144)
    parser.add_argument("--length", type=int, default=24)
    parser.add_argument("--parents", type=int, default=36)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    motif = "WLWKEFFKLAEFWLKEFFLWKAEF"[:args.length]

    lines = [
        "###########",
        "## EVOMD ##",
        "###########",
        "",
        "Date:",
        datetime.now().strftime("%Y/%m/%d %H:%M:%S"),
        "",
        "Description:",
        "\tAMP design run, objective=damage",
        "",
        "Settings:",
        "\tGenes:\t\t\t\t\t%s" % ALPHABET,
        "\tIterations:\t\t\t\t%d" % args.iterations,
        "\tPopulation Size:\t\t\t%d" % args.population,
        "\tNumber of parents:\t\t\t%d" % args.parents,
        "\tNumber of elites:\t\t\t1",
        "\tNumber of rerun elites:\t\t\t1",
        "\tRerun elites threshold:\t\t\t4",
        "\tScoring function:\t\t\t_compute_fitness_simulation",
        "\tSelection method:\t\t\tbest_N",
        "\tDiversity function:\t\t\t__default_diversity_function",
        "\tConstraint function:\t\t\tcombined",
        "\tP(crossover):\t\t\t\t1",
        "\tP_mutate function:\t\t\tdefault_P_mutate_function",
        "",
    ]

    run_counts = {}
    for it in range(args.iterations):
        progress = it / max(1, args.iterations - 1)
        # Saturating improvement, so the campaign plateaus rather than rising forever.
        conservation = 0.92 * (1 - pow(2.718281828, -3.2 * progress))
        pop = make_population(rng, args.population, args.length, motif, conservation)

        lines.append("Iteration %d:" % it)
        lines.append("finished:")
        lines.append(datetime.now().strftime("%Y/%m/%d %H:%M:%S"))

        scored = []
        for seq in pop:
            # Fitness correlates with how much of the motif the sequence carries, plus noise of
            # a magnitude comparable to what a 20 ns surface-tension average actually produces.
            match = sum(1 for a, b in zip(seq, motif) if a == b) / len(motif)
            raw = 1.2 + 7.5 * match + rng.gauss(0, 0.65)
            orientation = min(1.0, max(0.0, 1.0 - abs(rng.gauss(0, 0.06))))
            fitness = raw * orientation
            run_counts[seq] = run_counts.get(seq, 0) + 1
            scored.append((seq, fitness, raw, orientation, run_counts[seq]))

        scored.sort(key=lambda x: -x[1])
        for seq, fitness, raw, orientation, n_run in scored:
            lines.append("\t%s\t%s\tddF=%s\tpenalty=%s\tn_run=%d"
                         % (seq, fitness, raw, orientation, n_run))
        lines.append("")

    with open(args.out, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("wrote %s: %d iterations, %d evaluations"
          % (args.out, args.iterations, args.iterations * args.population))


if __name__ == "__main__":
    main()
