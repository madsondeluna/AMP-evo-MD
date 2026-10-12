#!/usr/bin/env python3
"""
Membrane preparation and calibration protocol.

Every membrane composition needs four reference measurements taken on the peptide-free bilayer
before any peptide can be evaluated: box area, free-membrane surface tension, bilayer thickness,
and the insertion depth derived from it. They are properties of a composition and are not
transferable between compositions, which is why user/membranes.py refuses to run without them.

This script is the protocol. It has two subcommands, run in sequence, with a GROMACS
equilibration between them:

    plan       compute lipid counts and water for a target box, emit the build commands
    calibrate  read the equilibrated system, measure the four constants, emit a MembraneDef

WHY THE PREPARATION IS ITERATIVE
--------------------------------
Two requirements pull against each other. Each membrane must be tensionless at its own
equilibrium area, otherwise the production run measures a pre-stressed bilayer and the
free-membrane surface tension is not zero. And the two membranes should occupy similar areas,
otherwise the per-peptide perturbation of gamma is diluted differently in each and the
comparison between them acquires a correction factor nobody wants to defend.

Both can be satisfied only by treating the lipid count, not the box, as the free parameter:
fix the target box, then choose how many lipids make that composition's equilibrium area land
there. Since the equilibrium area per lipid is only known after equilibrating, the first build
uses an estimate and the count is corrected on a second pass. Two passes are normally enough.

    plan  ->  build  ->  equilibrate (semi-isotropic, zero tension)  ->  measure area per lipid
          ->  plan again with the measured value  ->  build  ->  equilibrate
          ->  production at constant area  ->  calibrate

THE CONSTANT-AREA ENSEMBLE IS LOAD-BEARING
------------------------------------------
Production runs couple pressure only along z, with zero compressibility in x and y. The area is
therefore fixed and gamma responds to the peptide. Under a semi-isotropic barostat with a
reference tension the area would respond instead and gamma would sit pinned near its reference
value, and the observable this whole pipeline is built on would vanish. The semi-isotropic stage
exists only to find the equilibrium area; it is not the production ensemble.
"""

import argparse
import collections
import math
import os
import re
import sys


# Seed values only. These are starting guesses for the first build and are replaced by the
# measured value on the second pass. The only one anchored to a real Martini 3 measurement is
# POPC, taken from the Evo-MD reference membrane shipped in calibration/reference_popc
# (42.39 nm2 over 64 lipids per leaflet). The others are rough literature-scale estimates and
# should not be quoted anywhere: they exist to get the first box roughly right.
APL_SEED_NM2 = {
    "POPC": 0.662,   # measured, Martini 3, from the Evo-MD reference membrane
    "POPE": 0.60,    # estimate: smaller headgroup than PC
    "POPG": 0.65,    # estimate
    "POPS": 0.62,    # estimate
    "DPSM": 0.55,    # estimate: sphingomyelin packs tighter than PC
    "CHOL": 0.38,    # estimate: sterols occupy little area and condense their neighbours,
                     # so a linear sum over a cholesterol-rich mixture underestimates the
                     # condensing effect and the first pass will be further off than usual
    "CDL2": 1.24,    # estimate: cardiolipin carries four tails, roughly two lipids of area
}

# Martini water bead volume, used to turn a water count into a box height.
W_BEAD_VOLUME_NM3 = 0.12
# Waters per lipid, as beads. The Evo-MD reference uses 1800 W for 128 lipids.
DEFAULT_W_PER_LIPID = 14.0
# Helical rise per residue, for the box-size check.
HELIX_RISE_NM = 0.15


def parse_composition(text):
    """Parse 'POPE:70,POPG:20,CDL2:10' into an ordered {resname: fraction} summing to 1."""
    comp = collections.OrderedDict()
    for part in text.split(","):
        if ":" not in part:
            raise ValueError("composition entry %r should look like POPE:70" % part)
        name, value = part.split(":", 1)
        comp[name.strip().upper()] = float(value)
    total = sum(comp.values())
    if total <= 0:
        raise ValueError("composition fractions must sum to something positive")
    for name in comp:
        comp[name] /= total
    return comp


def mean_apl(comp, overrides):
    """
    Mole-fraction weighted mean area per lipid for a composition, in nm2.

    The special key MEAN overrides the whole mixture with one measured value, which is what the
    second planning pass uses: after equilibrating, the measured quantity is the mean area per
    lipid of the mixture as a whole, not a value per lipid type. Returns (value, measured) so the
    caller can label the output honestly rather than inferring it from whether any override was
    supplied.
    """
    if "MEAN" in overrides:
        return overrides["MEAN"], True

    total = 0.0
    missing = []
    measured = False
    for name, fraction in comp.items():
        if name in overrides:
            apl = overrides[name]
            measured = True
        else:
            apl = APL_SEED_NM2.get(name)
        if apl is None:
            missing.append(name)
            continue
        total += fraction * apl
    if missing:
        raise SystemExit(
            "No area-per-lipid value for: %s. Add it with --apl NAME=VALUE, pass one measured "
            "mixture value with --apl MEAN=VALUE, or add a seed to APL_SEED_NM2 in this file if "
            "it is a lipid you will use repeatedly."
            % ", ".join(missing)
        )
    return total, measured


def available_moleculetypes(force_field_dir=None):
    """Every moleculetype name defined by the force field files shipped for the active ff."""
    import re

    if force_field_dir is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        force_field_dir = os.path.join(here, "user", "base", "insert_peptide")

    try:
        import user.usersettings as UserSettings
        prefix = "martini_v3" if UserSettings.ff == "martini3001" else "martini_v2"
    except Exception:
        prefix = "martini_v3"

    pattern = re.compile(r"^\s*\[\s*moleculetype\s*\]\s*\n(?:\s*;[^\n]*\n)*\s*(\S+)", re.M)
    names = set()
    if not os.path.isdir(force_field_dir):
        return names, prefix
    for filename in sorted(os.listdir(force_field_dir)):
        if filename.startswith(prefix) and filename.endswith(".itp"):
            with open(os.path.join(force_field_dir, filename)) as handle:
                names |= set(pattern.findall(handle.read()))
    return names, prefix


def check_lipids_are_available(composition):
    """
    Stop if a requested lipid has no topology in the active force field.

    Without this the first symptom is a grompp error about an unknown moleculetype, after
    insane has already built a structure, and it names the lipid but not the fact that the
    force field shipped here simply does not define it. Both compositions in user/membranes.py
    hit this: the Martini 3 files here define neither cardiolipin (CDL0/1/2), nor sphingomyelin
    (DPSM), nor cholesterol (CHOL), so neither the bacterial nor the host membrane can be built
    as written without adding those parameter files.
    """
    available, prefix = available_moleculetypes()
    if not available:
        return

    missing = [lipid for lipid in composition if lipid not in available]
    if not missing:
        return

    sys.exit(
        "These lipids have no topology in the %s force field files shipped here: %s\n"
        "\n"
        "insane can place them, but grompp will refuse the system because no moleculetype\n"
        "defines them, after the build has already run.\n"
        "\n"
        "Either add the parameter file that defines them to user/base/insert_peptide/ and\n"
        "user/base/production/ (cardiolipin, sphingomyelin and sterols are distributed\n"
        "separately from the main Martini 3 lipid set), or choose a composition from the\n"
        "lipids that are defined here.\n"
        % (prefix, ", ".join(missing)))


def cmd_plan(args):
    comp = parse_composition(args.composition)
    check_lipids_are_available(comp)
    overrides = {}
    for item in args.apl or []:
        if "=" not in item:
            raise SystemExit("--apl expects NAME=VALUE (got %r)" % item)
        name, value = item.split("=", 1)
        name = name.strip().upper()
        if name != "MEAN" and name not in comp:
            raise SystemExit(
                "--apl %s refers to a lipid that is not in this composition (%s). Use a lipid "
                "name from the composition, or MEAN for one measured value for the whole "
                "mixture." % (name, ", ".join(comp))
            )
        try:
            overrides[name] = float(value)
        except ValueError:
            raise SystemExit("--apl %s: %r is not a number" % (name, value))

    area = args.box ** 2
    apl, measured = mean_apl(comp, overrides)
    per_leaflet = int(round(area / apl))
    counts = {}
    assigned = 0
    names = list(comp)
    # Assign by largest remainder so the totals match per_leaflet exactly rather than drifting
    # by a lipid or two after rounding, which would leave the two leaflets asymmetric.
    exact = {n: comp[n] * per_leaflet for n in names}
    for n in names:
        counts[n] = int(math.floor(exact[n]))
        assigned += counts[n]
    for n in sorted(names, key=lambda n: exact[n] - math.floor(exact[n]), reverse=True):
        if assigned >= per_leaflet:
            break
        counts[n] += 1
        assigned += 1

    total_lipids = 2 * per_leaflet
    n_water = int(round(total_lipids * args.w_per_lipid))
    water_height = (n_water * W_BEAD_VOLUME_NM3) / area
    membrane_height = args.membrane_thickness
    box_z = membrane_height + water_height

    helix = args.max_residues * HELIX_RISE_NM
    gap = args.box - helix

    print("# Membrane plan: %s" % args.name)
    print("#")
    print("# target box          %.2f x %.2f nm   (area %.2f nm2)" % (args.box, args.box, area))
    print("# mean area per lipid %.3f nm2 %s" % (apl, "(measured)" if measured else "(seed estimate)"))
    print("# lipids              %d per leaflet, %d total" % (per_leaflet, total_lipids))
    for n in names:
        print("#   %-6s %6.1f %%   %4d per leaflet" % (n, 100 * comp[n], counts[n]))
    print("# water               %d W beads (%.1f per lipid, %.2f nm of water)"
          % (n_water, args.w_per_lipid, water_height))
    print("# box z               %.2f nm (membrane %.2f + water %.2f)"
          % (box_z, membrane_height, water_height))
    print("#")
    print("# peptide fit         %d residues -> %.2f nm helix, %.2f nm to the periodic image"
          % (args.max_residues, helix, gap))
    if gap < 2.5:
        print("# WARNING: less than 2.5 nm between a fully extended helix and its own image.")
        print("#          A charged peptide will interact with itself across the boundary and")
        print("#          the measured surface tension will not be that of an isolated peptide.")
    print()

    lipid_arg = " ".join("-l %s:%d" % (n, counts[n]) for n in names if counts[n] > 0)
    print("# 1. build")
    print("insane -o %s_raw.gro -p %s.top -pbc rectangular \\" % (args.name, args.name))
    print("    -x %.4f -y %.4f -z %.4f \\" % (args.box, args.box, box_z))
    # -l sets the lower leaflet and -u the upper. The two strings are built from the same counts,
    # so the leaflets are symmetric by construction.
    upper_arg = lipid_arg.replace("-l ", "-u ")
    print("    %s %s -sol W -salt %s -center" % (lipid_arg, upper_arg, args.salt))
    print()
    print("# Note on -salt: this adds counterions plus the requested concentration. The")
    print("# anionic lipids (POPG, CDL2) need neutralising counterions regardless; the salt")
    print("# concentration on top of that is what makes the electrostatics physiological.")
    print("# The upstream Evo-MD systems carried no salt at all, which is the least")
    print("# favourable case for a charged peptide against a charged membrane.")
    print()
    print("# 2. equilibrate at zero tension to find the equilibrium area")
    print("#    Use a semi-isotropic barostat with ref_p = 1 1 bar (NOT the production npt.mdp,")
    print("#    which fixes the area). Run until the box area is flat, typically 100-200 ns.")
    print()
    print("# 3. measure the equilibrium area per lipid and re-plan")
    print("python tools/membrane.py plan --name %s --composition '%s' --box %.2f \\"
          % (args.name, args.composition, args.box))
    print("    --apl <NAME>=<measured value per lipid type, or pass one mean with --apl MEAN=x>")
    print("#    In practice: read the equilibrated area A from the .gro, divide by the lipid")
    print("#    count per leaflet to get the mean APL, and re-run with that. The lipid count")
    print("#    changes; the box does not.")
    print()
    print("# 4. production reference run at constant area, using the production npt.mdp")
    print("#    This is the peptide-free run whose .edr gives the free-membrane surface tension.")
    print()
    print("# 5. calibrate")
    print("python tools/membrane.py calibrate --name %s --gro %s_eq.gro --energy %s_free.ENERGY"
          % (args.name, args.name, args.name))


def read_gro(path):
    """Return (box_x, box_y, box_z, atoms) where atoms is a list of (resname, atomname, z)."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    if len(lines) < 3:
        raise SystemExit("%s does not look like a .gro file" % path)
    n_atoms = int(lines[1].strip())
    atoms = []
    for line in lines[2:2 + n_atoms]:
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            z = float(line[36:44])
        except ValueError:
            continue
        atoms.append((resname, atomname, z))
    box = lines[2 + n_atoms].split()
    return float(box[0]), float(box[1]), float(box[2]), atoms



# Reference bead for counting a lipid. One bead per lipid is what every count here rests on, and
# a plain name match breaks it in two directions: cardiolipin carries PO41 and PO42, so matching
# "PO4" exactly misses it while matching the prefix counts it twice, and sterols carry no
# phosphate at all.
#
# The rule is by bead-name shape rather than by a list of residue names. A whitelist is the wrong
# shape for this: the shipped lipid library defines CDL0, CDL1 and CDL2, and a list that named
# only two of them would count the third twice, silently and in the direction that inflates the
# lipid count. Any future renamed or added multi-phosphate lipid is handled the same way.
STEROL_BEAD = "ROH"
# PO4 followed by a digit: the numbered phosphates of a multi-phosphate lipid.
NUMBERED_PHOSPHATE = re.compile(r"^(PO4)(\d+)$")


def reference_bead(atom, head_bead="PO4"):
    """
    Whether this bead is the one that counts its lipid, by name alone.

    The unnumbered head bead (PO4) for an ordinary phospholipid, the sterol hydroxyl (ROH) for a
    sterol, and only the FIRST numbered phosphate (PO41, not PO42) for a multi-phosphate lipid.
    """
    if atom == head_bead or atom == STEROL_BEAD:
        return True
    m = NUMBERED_PHOSPHATE.match(atom)
    return bool(m) and m.group(1) == head_bead and m.group(2) == "1"


def count_lipids(atoms, head_bead="PO4"):
    """
    Lipids per residue name, and the z position of each one's reference bead.

    Returns (counts, zs). Counting is by reference bead so that a lipid with two phosphate groups
    is counted once and a sterol with no phosphate is counted at all.
    """
    counts = collections.Counter()
    zs = []
    for res, atom, z in atoms:
        if reference_bead(atom, head_bead):
            counts[res] += 1
            zs.append(z)
    return counts, zs


def measure_thickness(atoms, head_bead="PO4"):
    """
    Bilayer thickness: the distance between the mean z of the phosphate beads in each leaflet.

    Measured on phosphates only, which is the right reference for insertion depth, so sterols are
    deliberately excluded. Cardiolipin contributes its first phosphate only, since the two sit at
    the same depth and counting both would weight it double.
    """
    zs = []
    for _res, atom, z in atoms:
        # Phosphates only: sterols are deliberately excluded, since insertion depth is referenced
        # to the phosphate plane. A multi-phosphate lipid contributes its first phosphate only,
        # because both sit at the same depth and counting them twice would weight it double.
        if atom == head_bead:
            zs.append(z)
            continue
        m = NUMBERED_PHOSPHATE.match(atom)
        if m and m.group(1) == head_bead and m.group(2) == "1":
            zs.append(z)
    if len(zs) < 2:
        raise SystemExit(
            "Found %d %s beads. Thickness cannot be measured without head beads in both "
            "leaflets; pass --head-bead if this composition uses a different head group."
            % (len(zs), head_bead)
        )
    midplane = sum(zs) / len(zs)
    upper = [z for z in zs if z >= midplane]
    lower = [z for z in zs if z < midplane]
    if not upper or not lower:
        raise SystemExit("All %s beads fell on one side of the midplane; is the membrane "
                         "split across the periodic boundary? Centre it with trjconv first."
                         % head_bead)
    return (sum(upper) / len(upper)) - (sum(lower) / len(lower)), len(upper), len(lower)


def parse_gmx_energy(path, term="#Surf*SurfTen"):
    """
    Pull (mean, err) for a term out of a gmx energy stdout dump.

    The statistics table gmx prints after the interactive selection has the shape
        <term name>   <average>  <err.est.>  <rmsd>  <tot-drift>
    and the term name contains spaces, so match by prefix rather than splitting blindly.
    """
    try:
        text = open(path).read()
    except OSError as exc:
        raise SystemExit("Could not read %s: %s" % (path, exc))
    # Only surface-tension row labels. Pres-XX used to be in this list as a fallback, which
    # meant that an energy dump without a surface-tension row would silently yield a
    # pressure-tensor component, in bar, written straight into the MembraneDef as if it were the
    # surface tension. Failing is the only safe behaviour here.
    for line in text.splitlines():
        stripped = line.strip()
        for candidate in (term, "Surf*SurfTen", "#Surf*SurfTen"):
            if stripped.startswith(candidate):
                rest = stripped[len(candidate):].split()
                if len(rest) >= 2:
                    try:
                        return float(rest[0]), float(rest[1])
                    except ValueError:
                        continue
    raise SystemExit(
        "No surface tension row found in %s. Produce it with:\n"
        "    echo '#Surf*SurfTen' | gmx energy -f <free membrane>.edr -b <equilibrated time> "
        "-o surften.xvg > free.ENERGY 2>&1" % path
    )


def cmd_calibrate(args):
    box_x, box_y, box_z, atoms = read_gro(args.gro)
    area = box_x * box_y

    # Count lipids by residue, not by bead name: cardiolipin has two phosphates and sterols have
    # none, so a plain bead match would miss one and double the other.
    head_counts, _zs = count_lipids(atoms, args.head_bead)
    n_lipids = sum(head_counts.values())

    thickness, n_upper, n_lower = measure_thickness(atoms, args.head_bead)
    gamma_mean, gamma_err = parse_gmx_energy(args.energy)

    from importlib.machinery import SourceFileLoader
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    b_factor_src = os.path.join(here, "user", "membranes.py")
    b_from_thickness = 0.5 * 0.92
    try:
        src = open(b_factor_src).read()
        m = re.search(r"B_FROM_THICKNESS\s*=\s*([0-9.]+)\s*\*\s*([0-9.]+)", src)
        if m:
            b_from_thickness = float(m.group(1)) * float(m.group(2))
    except OSError:
        pass
    b_nm = b_from_thickness * thickness

    print("# Calibration: %s" % args.name)
    print("#")
    print("# box            %.4f x %.4f x %.4f nm" % (box_x, box_y, box_z))
    print("# area           %.4f nm2" % area)
    print("# lipids         %d total (%d / %d phosphate beads per leaflet)"
          % (n_lipids, n_upper, n_lower))
    if n_upper + n_lower != n_lipids:
        print("#                (phosphate count differs from the lipid count because the")
        print("#                 composition contains sterols, which have none)")
    if n_lipids:
        print("# area per lipid %.4f nm2" % (area / (n_lipids / 2.0)))
    for res, count in sorted(head_counts.items()):
        print("#   %-6s %4d" % (res, count))
    print("# thickness      %.4f nm (head bead to head bead)" % thickness)
    print("# #Surf*SurfTen  %.4f +/- %.4f bar nm" % (gamma_mean, gamma_err))
    if abs(gamma_mean) > 3 * max(gamma_err, 1e-9) and abs(gamma_mean) > 10:
        print("#")
        print("# WARNING: the free membrane is not tensionless. The Evo-MD POPC reference sits at")
        print("#          1.93 +/- 0.99 bar nm, i.e. zero within error. A large value here means the")
        print("#          constant-area box was not set to this composition's equilibrium area,")
        print("#          so every peptide will be evaluated against a pre-stressed bilayer.")
        print("#          Re-run the semi-isotropic equilibration and re-plan the lipid count.")
    if n_upper != n_lower:
        print("#")
        print("# NOTE: leaflets are asymmetric (%d vs %d). Intentional for an outer-membrane"
              % (n_upper, n_lower))
        print("#       model, a problem otherwise: an area mismatch between leaflets shows up as")
        print("#       spontaneous curvature and contaminates the surface tension.")
    print()
    print("# Paste into MEMBRANES in user/membranes.py:")
    print()
    print('    "%s": MembraneDef(' % args.name)
    print('        name="%s",' % args.name)
    print('        gro="membrane_%s.gro",' % args.name)
    print('        top="membrane_%s.top",' % args.name)
    print('        role="%s",' % args.role)
    print("        area_nm2=%.4f," % area)
    print("        surften_free=(%.5f, %.2f)," % (gamma_mean, gamma_err))
    print("        thickness_nm=%.4f,   # b_nm derives to %.4f" % (thickness, b_nm))
    print("        composition={%s}," % ", ".join('"%s": %d' % (r, c) for r, c in sorted(head_counts.items())))
    print('        notes="",')
    print("    ),")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("plan", help="compute lipid counts for a target box and emit build commands")
    p.add_argument("--name", required=True, help="membrane name, e.g. bacterial")
    p.add_argument("--composition", required=True,
                   help="e.g. 'POPE:70,POPG:20,CDL2:10' (values are relative, need not sum to 100)")
    p.add_argument("--box", type=float, default=8.0, help="target lateral box in nm (default 8.0)")
    p.add_argument("--max-residues", type=int, default=30,
                   help="longest peptide the box must accommodate (default 30)")
    p.add_argument("--w-per-lipid", type=float, default=DEFAULT_W_PER_LIPID,
                   help="water beads per lipid (default %.0f, as in the Evo-MD reference)"
                        % DEFAULT_W_PER_LIPID)
    p.add_argument("--membrane-thickness", type=float, default=4.2,
                   help="approximate membrane height in nm, for sizing the box in z (default 4.2)")
    p.add_argument("--salt", default="0.15", help="salt concentration in M (default 0.15)")
    p.add_argument("--apl", action="append", metavar="NAME=VALUE",
                   help="measured area per lipid in nm2, overriding the seed estimate; repeatable")
    p.set_defaults(func=cmd_plan)

    c = sub.add_parser("calibrate", help="measure the reference constants from an equilibrated system")
    c.add_argument("--name", required=True)
    c.add_argument("--gro", required=True, help="equilibrated peptide-free membrane, centred")
    c.add_argument("--energy", required=True,
                   help="gmx energy stdout dump containing the #Surf*SurfTen statistics row")
    c.add_argument("--role", default="target", choices=("target", "host", "probe"))
    c.add_argument("--head-bead", default="PO4", help="head group bead name (default PO4)")
    c.set_defaults(func=cmd_calibrate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
