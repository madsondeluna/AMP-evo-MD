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
                        (mean, stderr) in bar
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
    """

    def __init__(self, name, gro, top, role,
                 area_nm2=None, surften_free=None, thickness_nm=None, b_nm=None,
                 composition=None, notes=""):
        self.name = name
        self.gro = gro
        self.top = top
        self.role = role
        self.area_nm2 = area_nm2
        self.surften_free = surften_free      # (mean, stderr) in bar
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

    "bacterial": MembraneDef(
        name="bacterial",
        gro="membrane_bacterial.gro",
        top="membrane_bacterial.top",
        role="target",
        composition={"POPE": 0.70, "POPG": 0.20, "CDL2": 0.10},
        notes="Gram-negative inner membrane mimic. Anionic, no sterol.",
    ),

    "host": MembraneDef(
        name="host",
        gro="membrane_host.gro",
        top="membrane_host.top",
        role="host",
        composition={"POPC": 0.60, "DPSM": 0.10, "CHOL": 0.30},
        notes=("Erythrocyte outer-leaflet mimic. The sterol fraction is the main "
               "determinant of resistance to lysis and is what makes this a hemolysis "
               "proxy rather than a generic neutral membrane."),
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
