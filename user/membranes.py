"""
Membrane registry.

In the original Evo-MD the membrane was implicit: a single POPC composition held at two
areas, with its reference quantities written as four bare constants in usersettings.py
(Area_free_low_tension, Area_free_high_tension, SurfTen_free_low_tension,
SurfTen_free_high_tension) plus the insertion depths b.

Here the axis of comparison is composition rather than tension, so every membrane needs its
own copy of those quantities. They are reference measurements on the PEPTIDE-FREE membrane
and they are not transferable between compositions: a value measured on POPC is meaningless
for a POPE/POPG/cardiolipin mixture.

A MembraneDef with any reference field left as None is considered uncalibrated and raises on
use. This is deliberate. A silent placeholder here would propagate into every fitness value
in the run and the numbers would look perfectly reasonable while meaning nothing.

Calibration procedure for each composition, run once, before any peptide:
  1. Build and equilibrate the bare membrane at the target area (constant-area ensemble,
     the same npt.mdp the production runs use).
  2. area_nm2        <- x * y of the equilibrated box (last line of the .gro)
  3. surften_free    <- gmx energy on "#Surf*SurfTen" over the same averaging window the
                        production runs use (UserSettings.ENERGY_init_time onward), as
                        (mean, stderr) in bar*nm, which is the unit of the
                        #Surf*SurfTen column gmx energy reports, not bar
  4. thickness_nm    <- bilayer thickness (phosphate-to-phosphate, or the equivalent
                        headgroup bead distance for lipids without phosphate)
  5. b_nm            <- insertion reference depth; see the note on B_FROM_THICKNESS below
"""

import os

# usersettings imports this module to build its membrane_types list, so importing it back at
# module level would be circular. The two places that need it do a local import instead.


# The original code sets b = 0.5 * 0.92 * thickness (usersettings.py comment:
# "thickness values: 3.601, 3.977 (b_value = 0.5*0.92*thickness)"), i.e. 92% of the
# monolayer thickness. The 0.92 factor does not appear in the main text of van Hilten et al.
# 2023, which quotes b as the monolayer thickness itself (1.90 and 1.79 nm). The factor is
# kept here to stay numerically faithful to the code that produced the published results,
# but it is exposed rather than buried so it can be checked against section S1 and changed
# in one place.
B_FROM_THICKNESS = 0.5 * 0.92


class UncalibratedMembrane(RuntimeError):
    """Raised when a membrane is used before its reference quantities have been measured."""


class MembraneDef:
    """
    One membrane composition plus the reference measurements taken on it without peptide.

    role:
        "target" -> the membrane the peptide is meant to act on (bacterial)
        "host"   -> the membrane it must spare (eukaryotic / erythrocyte mimic)
        "probe"  -> anything else carried along for comparison, ignored by the objectives

    composition:
        LIPID COUNTS in the built system, {resname: number of molecules}, as
        tools/membrane.py calibrate emits them. Not mole fractions.

        This field had two meanings in the repository at once: the placeholder entries below
        were written as fractions summing to 1.0, while the protocol emits counts, and the one
        piece of code that consumes it derives an area per lipid from the sum. With fractions
        that gives 2 * area_nm2, which is off by two orders of magnitude and silently so. The
        counts are the useful form, because area_nm2 divided by them is a quantity that can be
        checked against what a bilayer of that composition should have.

        The intended ratio belongs in notes, where it cannot be mistaken for a measurement.
    """

    def __init__(self, name, gro, top, role,
                 area_nm2=None, surften_free=None, thickness_nm=None, b_nm=None,
                 composition=None, notes=""):
        self.name = name
        self.gro = gro
        self.top = top
        self.role = role
        self.area_nm2 = area_nm2
        # (mean, stderr) in bar*nm, the unit of the #Surf*SurfTen column. Not bar: writing
        # "bar" here invited dividing by the box height to "get" bar, which would be wrong.
        self.surften_free = surften_free
        self.thickness_nm = thickness_nm
        self._b_nm = b_nm
        self.composition = composition or {}
        self.notes = notes

        if role not in ("target", "host", "probe"):
            raise ValueError("role must be one of target, host, probe (got %r)" % role)

    @property
    def b_nm(self):
        # Explicit b wins; otherwise derive it from the measured thickness.
        if self._b_nm is not None:
            return self._b_nm
        if self.thickness_nm is not None:
            return B_FROM_THICKNESS * self.thickness_nm
        return None

    @property
    def is_calibrated(self):
        return (self.area_nm2 is not None
                and self.surften_free is not None
                and self.b_nm is not None)

    def require_calibrated(self):
        if self.is_calibrated:
            return
        missing = [field for field, value in (
            ("area_nm2", self.area_nm2),
            ("surften_free", self.surften_free),
            ("b_nm (or thickness_nm)", self.b_nm),
        ) if value is None]
        raise UncalibratedMembrane(
            "Membrane %r is missing reference measurements: %s. "
            "These must be measured on the peptide-free membrane before any fitness "
            "evaluation; see the calibration procedure in user/membranes.py."
            % (self.name, ", ".join(missing))
        )

    def surften_free_converted(self):
        """Free-membrane surface tension in kJ mol-1 nm-2, with its error."""
        import user.usersettings as UserSettings
        self.require_calibrated()
        mean, err = self.surften_free
        f = UserSettings.conversion_factor
        return mean * f, err * f

    def gro_path(self):
        import user.usersettings as UserSettings
        return os.path.join(UserSettings.path_basedir_insert_peptide, self.gro)

    def top_path(self):
        import user.usersettings as UserSettings
        return os.path.join(UserSettings.path_basedir_insert_peptide, self.top)

    def __repr__(self):
        state = "calibrated" if self.is_calibrated else "UNCALIBRATED"
        return "<MembraneDef %s role=%s %s>" % (self.name, self.role, state)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
# Edit the compositions to match what you actually build. The reference numbers stay None
# until measured; the run will refuse to compute a fitness before then.
#
# The compositions below are a starting point, not a recommendation that survived review:
# they are the common literature mimics and the lipid choices should be argued explicitly in
# the manuscript. Note in particular that neither includes LPS, so "bacterial" here means
# the inner membrane of a Gram-negative or the cytoplasmic membrane of a Gram-positive, and
# any claim about Gram-negative potency has to account for the outer membrane separately.

MEMBRANES = {

    # Calibrated references for this study. The reference surface tension is measured on the
    # peptide-free membrane in the constant-area production ensemble (tools/protocol.sh PASS=2),
    # holding the box at each composition's equilibrium area, so the free membrane is tensionless
    # within the method's noise and the peptide runs are not scored against a pre-stressed bilayer.
    "bacterial": MembraneDef(
        name="bacterial",
        gro="membrane_bacterial.gro",
        top="membrane_bacterial.top",
        role="target",
        area_nm2=62.2264,
        surften_free=(-2.61675, 11.00),
        thickness_nm=3.9312,   # b_nm derives to 1.8084
        composition={"POPE": 140, "POPG": 58},
        notes="Gram-negative inner membrane mimic. Anionic, no sterol. POPE:POPG 70:30.",
    ),

    "host": MembraneDef(
        name="host",
        gro="membrane_host.gro",
        top="membrane_host.top",
        role="host",
        area_nm2=66.2772,
        surften_free=(-9.54955, 4.10),
        thickness_nm=3.8710,   # b_nm derives to 1.7807
        composition={"POPC": 200},
        notes=("Neutral zwitterionic host (pure POPC). Not the sterol-containing erythrocyte "
               "mimic; POPC:DPSM:CHOL upgrade planned."),
    ),
}


def membrane_names():
    """Names in a fixed order, so per-membrane filenames and indices stay reproducible."""
    return sorted(MEMBRANES)


def get(name):
    try:
        return MEMBRANES[name]
    except KeyError:
        raise KeyError("Unknown membrane %r. Registered: %s"
                       % (name, ", ".join(membrane_names())))


def by_role(role):
    """All membranes with a given role, in registry order."""
    return [MEMBRANES[n] for n in membrane_names() if MEMBRANES[n].role == role]


def entry_for(observables, membrane):
    """
    Look up a membrane's entry in an observables dict, explaining a mismatch.

    The objectives take the membrane names from this registry and the measurements from a file
    written by an earlier run, and the two can disagree: re-scoring an archived campaign after
    the registry has changed is the obvious way, and it used to fail as a bare KeyError naming
    a membrane the operator may never have heard of.
    """
    import core.errors as Errors

    try:
        return observables[membrane.name]
    except KeyError:
        present = [k for k in observables if not k.startswith("_")]
        raise Errors.OperatorError(
            "these results were not produced with the membranes currently registered",
            where=("registry (user/membranes.py) expects: %s\n"
                   "the observables file contains:        %s"
                   % (", ".join(sorted(membrane_names())), ", ".join(sorted(present)) or "nothing")),
            why=("the objective reads the %s membrane by name, and this result set does not "
                 "have one called %r" % (membrane.role, membrane.name)),
            fix=("re-score with the registry these results were produced under: "
                 "calibrate.py --membranes <file defining MEMBRANES>, or restore the entries "
                 "in user/membranes.py"))


def load_registry(path):
    """
    Replace the registry with the MEMBRANES dict defined in a Python file.

    For re-scoring a result set whose membranes are no longer the ones in the repository. The
    file is executed, so it is code and is treated as such: it is the operator's own archived
    copy of a registry, not untrusted input.
    """
    import core.errors as Errors

    Errors.require_file(
        path,
        why="--membranes loads the membrane registry a result set was produced under",
        fix=("point it at a file defining MEMBRANES = {...} with MembraneDef entries, such as "
             "an archived copy of user/membranes.py"))

    namespace = {"MembraneDef": MembraneDef, "__name__": "archived_membranes"}
    with open(path) as handle:
        source = handle.read()
    try:
        exec(compile(source, path, "exec"), namespace)                 # noqa: S102
    except Exception as exc:                                           # noqa: BLE001
        Errors.fail("could not load the membrane registry (%s: %s)" % (type(exc).__name__, exc),
                    where=os.path.abspath(path),
                    why="the file is executed as Python and must define MEMBRANES",
                    fix="check it defines MEMBRANES = {name: MembraneDef(...), ...}")

    registry = namespace.get("MEMBRANES")
    if not isinstance(registry, dict) or not registry:
        Errors.fail("the file defines no usable MEMBRANES dict",
                    where=os.path.abspath(path),
                    why="a registry is a non-empty dict of name to MembraneDef",
                    fix="define MEMBRANES = {name: MembraneDef(...), ...}")

    MEMBRANES.clear()
    MEMBRANES.update(registry)

    # usersettings.membrane_types is a snapshot of membrane_names() taken at import, and it is
    # what the measurement modules loop over while the objectives read MEMBRANES. Leaving it
    # stale means a run after a registry swap would measure one set of membranes and score
    # another. Harmless while only calibrate.py swaps registries, since it does not measure,
    # and a trap the first time anything else does.
    try:
        import user.usersettings as UserSettings
        UserSettings.membrane_types = membrane_names()
    except ImportError:
        pass

    return registry


def the_target():
    """The single target membrane. The objectives assume there is exactly one."""
    targets = by_role("target")
    if len(targets) != 1:
        raise ValueError("Expected exactly one membrane with role='target', found %d: %s"
                         % (len(targets), [m.name for m in targets]))
    return targets[0]


def the_host():
    """The single host membrane, or None if the run does not model one."""
    hosts = by_role("host")
    if len(hosts) > 1:
        raise ValueError("Expected at most one membrane with role='host', found %d: %s"
                         % (len(hosts), [m.name for m in hosts]))
    return hosts[0] if hosts else None


def require_all_calibrated():
    """Fail before the first simulation rather than after it."""
    for name in membrane_names():
        MEMBRANES[name].require_calibrated()


def require_all_built():
    """
    Check that the structure and topology each MembraneDef names are actually on disk.

    require_all_calibrated above checks the four reference numbers, which the operator types in
    by hand from the protocol output. It does not check that the equilibrated system those
    numbers describe was ever copied into place, and the two are separate steps that are easy to
    do one of. Without this, the first thing that reads the file is insert_peptide_nmem, long
    after the peptide has been built, and it reports a path inside a temporary directory that no
    longer exists by the time anyone looks:

        FileNotFoundError: '/tmp/tmpk3x9r2a1/membrane_bacterial.gro'

    which names neither the membrane nor the directory the file belongs in.
    """
    import core.errors as Errors

    missing = []
    for name in membrane_names():
        membrane = MEMBRANES[name]
        for kind, path in (("structure", membrane.gro_path()), ("topology", membrane.top_path())):
            if not os.path.isfile(path):
                missing.append("%s (%s of membrane %r)" % (path, kind, name))

    if missing:
        Errors.fail("membrane files named in user/membranes.py are not on disk",
                    where="\n".join(missing),
                    why=("the reference measurements for these membranes are filled in, so the "
                         "equilibrated system they describe was built, but it was not copied "
                         "into place"),
                    fix=("copy the equilibrated .gro and .top there; this is the last step "
                         "tools/protocol.sh prints when it finishes"))

    # area_nm2 is typed in by hand from the protocol output, and it is the denominator of
    # delta_gamma_work. Nothing compared it against the box of the structure it describes, so a
    # mistyped digit, or a membrane rebuilt at a different area after the number was recorded,
    # gave a wrong answer with no symptom: at constant area the reference surface tension was
    # also measured at the other area, so delta_gamma itself stops meaning anything.
    import core.gromacs as Gromacs

    mismatched = []
    for name in membrane_names():
        membrane = MEMBRANES[name]
        try:
            box = Gromacs.GRO_read_box(membrane.gro_path())
            # GRO_read_box returns strings. Converted once here, because formatting them as
            # floats further down is what the error path itself used to fail on.
            box_x, box_y = float(box[0]), float(box[1])
            measured = box_x * box_y
        except Exception:
            # Reading the box is a check, not a requirement; a file this helper cannot parse
            # should not stop a run that would otherwise work. The conversion is inside the
            # guard too, because an unparseable box line fails there rather than in the read.
            continue
        if abs(measured - membrane.area_nm2) > 0.01 * max(measured, 1e-9):
            mismatched.append(
                "%s: area_nm2 = %.4f, but %s has a box of %.4f x %.4f = %.4f nm2 (%.1f%% apart)"
                % (name, membrane.area_nm2, os.path.basename(membrane.gro_path()),
                   box_x, box_y, measured,
                   100.0 * abs(measured - membrane.area_nm2) / max(measured, 1e-9)))

    if mismatched:
        Errors.fail("a membrane's recorded area does not match the structure it names",
                    where="\n".join(mismatched),
                    why=("area_nm2 divides the work term, and the reference surface tension was "
                         "measured on a box of one particular area. A structure at a different "
                         "area makes delta_gamma a difference between two unrelated states."),
                    fix=("re-run tools/protocol.sh for that membrane and paste the MembraneDef "
                         "it emits, or copy the structure the current numbers were measured on"))
