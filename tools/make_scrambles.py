"""
Generate composition-matched, non-amphipathic scramble controls for the calibration set.

Each inactive control is a permutation of an active parent that minimises the Eisenberg
hydrophobic moment of an ideal alpha-helix (100 deg per residue). The permutation preserves the
amino-acid composition exactly and destroys the amphipathic pattern, which is what makes it a
negative control: same residues, no face. Terminal proline is rejected because martinize2 cannot
assign it.

The Eisenberg scale and the moment are identical to tools/ampstyle.py, so the descriptor matches
the one the analysis figures use. The search is seeded per parent (base string "20231023:" plus
the parent sequence), so each scramble is reproducible on its own and independent of the others.

Usage:
    python tools/make_scrambles.py
It prints the parent, the parent moment, the scramble, and the scramble moment, and checks that
the composition is preserved. The scrambles it prints for magainin1, melittin and V681 are the
ones recorded in seq.txt and labels.csv as scram-magainin1, scram-melittin and scram-V681.
"""

import math
import random

# Eisenberg consensus hydrophobicity scale, identical to tools/ampstyle.py _EISENBERG.
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74,
    "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19,
    "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}


def hydrophobic_moment(sequence, angle_deg=100.0):
    """Eisenberg hydrophobic moment per residue, as in tools/ampstyle.py."""
    ang = math.radians(angle_deg)
    sx = sum(EISENBERG.get(a, 0.0) * math.cos(i * ang) for i, a in enumerate(sequence))
    sy = sum(EISENBERG.get(a, 0.0) * math.sin(i * ang) for i, a in enumerate(sequence))
    return math.hypot(sx, sy) / len(sequence)


def scramble(parent, trials=50000):
    """Permutation of parent with the lowest hydrophobic moment, no terminal proline."""
    rng = random.Random("20231023:" + parent)
    chars = list(parent)
    best, best_mu = None, float("inf")
    for _ in range(trials):
        rng.shuffle(chars)
        if chars[0] == "P" or chars[-1] == "P":
            continue
        candidate = "".join(chars)
        mu = hydrophobic_moment(candidate)
        if mu < best_mu:
            best_mu, best = mu, candidate
    return best, best_mu


# Parents whose scrambles this script generates. The first block balanced the core set to 6 active
# / 6 inactive (the earlier scram-magainin2, scram-PGLa and scram-pleurocidin are kept as recorded
# in seq.txt). The second block is the expansion panel (seq_extra.txt): four additional actives,
# all proline-free and alpha-helical, with one composition-matched scramble each, taking the
# calibration set to 10 active / 10 inactive. Verified sequences:
#   aurein 1.2   GLFDIIKKIAESF            selective   Litoria aurea      (C-amidated)
#   citropin 1.1 GLFDVIKKVASVIGGL         selective   Litoria citropa
#   pexiganan    GIGKFLKKAKKFGKAFVKILKK   selective   MSI-78, magainin analog (C-amidated)
#   mastoparan   INLKALAALAKKIL           lytic       Vespula lewisii
PARENTS = {
    "scram-magainin1":   "GIGKFLHSAGKFGKAFVGEIMKS",
    "scram-melittin":    "GIGAVLKVLTTGLPALISWIKRKRQQ",
    "scram-V681":        "KWKSFLKTFKSAVKTVLHTALKAISS",
    "scram-aurein1.2":   "GLFDIIKKIAESF",
    "scram-citropin1.1": "GLFDVIKKVASVIGGL",
    "scram-pexiganan":   "GIGKFLKKAKKFGKAFVKILKK",
    "scram-mastoparan":  "INLKALAALAKKIL",
}


def main():
    print("%-16s %-28s %8s   %-28s %8s  comp_ok"
          % ("name", "parent", "mu_H", "scramble", "mu_H"))
    for name, parent in PARENTS.items():
        s, mu = scramble(parent)
        comp_ok = sorted(s) == sorted(parent)
        print("%-16s %-28s %8.4f   %-28s %8.4f  %s"
              % (name, parent, hydrophobic_moment(parent), s, mu, comp_ok))


if __name__ == "__main__":
    main()
