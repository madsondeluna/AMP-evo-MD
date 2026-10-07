"""
Record which tools actually produced a run, and refuse to continue if they changed.

The motivation is specific. On a machine with environment modules, `module load gromacs` without
a version loads whatever the site has set as default, and administrators change that default
without announcing it. A campaign that runs for days across several sessions can therefore be
produced by more than one GROMACS build, and nothing in the output would record it: the surface
tension values would simply shift, and the shift would be indistinguishable from a real effect
of the sequences.

Intending to use one version is not the same as having used one. This module makes it checkable:

  capture()        what is on PATH right now, as a dict
  format_block()   that dict as lines for the progress file header
  check_expected() compare against a pinned version and fail loudly on a mismatch

Set GMX_VERSION_EXPECT to pin. With it unset nothing is enforced, only recorded, because
enforcing a version nobody declared would block a first run for no reason.
"""

import os
import re
import subprocess
import sys


def _run(command, timeout=30):
    """Return stdout for a command, or None if it cannot be run."""
    try:
        result = subprocess.run(command, shell=True, capture_output=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        # A non-zero exit means the command did not do what was asked, whatever it printed. The
        # earlier version accepted any stdout, so a wrapper printing usage text and exiting 1
        # was parsed as a banner, found no version line, and recorded "UNKNOWN" for the whole
        # campaign with nothing to indicate the probe had failed.
        return None
    # Several martinize2 releases write the version to stderr rather than stdout.
    out = (result.stdout or b"").decode("utf-8", errors="replace")
    err = (result.stderr or b"").decode("utf-8", errors="replace")
    return out + ("\n" + err if err.strip() else "")


def gmx_version(binary=None):
    """
    Parse the version out of `gmx --version`.

    Returns (version, full_banner). The banner carries the build host, compiler, SIMD
    instruction set and GPU support, all of which change the numbers this pipeline measures and
    none of which are recoverable afterwards from the .edr files.
    """
    binary = binary or os.environ.get("GMX", "gmx")
    out = _run("%s --version" % binary)
    if out is None:
        return None, None

    version = None
    for line in out.splitlines():
        m = re.match(r"\s*GROMACS version:?\s+(.+?)\s*$", line)
        if m:
            version = m.group(1).strip()
            break

    # Everything with a "key: value" shape, rather than a whitelist, so that the CUDA driver and
    # runtime versions and the FFT library are recorded: all of them change what this pipeline
    # measures, and an earlier whitelist dropped exactly those.
    #
    # The skip list is not cosmetic. The banner exists so that a build change can be spotted by
    # comparing progress-file headers across sessions, and anything that varies between runs of
    # the SAME build defeats that. GROMACS prints a random quotation on every invocation, and the
    # paths below differ by where the run happened rather than by what was run.
    skip = ("Executable", "Data prefix", "Working dir", "Command line",
            "GROMACS is written by", "GROMACS reminds you")
    banner = []
    for line in out.splitlines():
        stripped = line.strip()
        if not stripped or ":" not in stripped:
            continue
        if any(stripped.startswith(k) for k in skip):
            continue
        # The ":-)" in the GROMACS header art contains a colon but is not a key/value line.
        if stripped.startswith(":-)") or ":-)" in stripped:
            continue
        banner.append(stripped)
    return version, banner


def martinize2_version(binary=None):
    binary = binary or os.environ.get("MARTINIZE2", "martinize2")
    out = _run("%s --version" % binary)
    if out is None:
        return None
    return out.strip().splitlines()[0] if out.strip() else None


def capture():
    """Everything worth recording about the environment a run was produced in."""
    gmx_binary = os.environ.get("GMX", "gmx")
    version, banner = gmx_version(gmx_binary)

    info = {
        "gmx_binary": gmx_binary,
        "gmx_version": version or "UNKNOWN",
        "gmx_banner": banner or [],
        "martinize2_version": martinize2_version() or "UNKNOWN",
        "python": sys.version.split()[0],
    }

    # The loaded module list, when the machine uses environment modules. This is what actually
    # determines which GROMACS was picked up, and it is not reconstructible later.
    loaded = os.environ.get("LOADEDMODULES")
    if loaded:
        info["modules"] = loaded.split(":")

    for package in ("numpy", "vermouth", "mpi4py"):
        try:
            module = __import__(package)
            info[package] = getattr(module, "__version__", "unknown")
        except ImportError:
            info[package] = "absent"

    return info


def check_expected(info=None, expected=None):
    """
    Compare the running GROMACS against a pinned version and stop on a mismatch.

    Pinned through GMX_VERSION_EXPECT. Matching is a prefix test, so "2024" accepts
    "2024.3-Ubuntu" and "2024.3" does not accept "2024.4". Unset means record only.
    """
    info = info or capture()
    expected = expected or os.environ.get("GMX_VERSION_EXPECT")
    if not expected:
        return True

    found = info.get("gmx_version", "UNKNOWN")
    if found.startswith(expected):
        return True

    sys.stderr.write(
        "\nGROMACS version mismatch.\n"
        "  expected (GMX_VERSION_EXPECT): %s\n"
        "  found:                         %s\n"
        "  binary:                        %s\n\n"
        "Surface tension is sensitive to the build, so mixing versions within one campaign makes\n"
        "the results incomparable in a way that leaves no trace in the output. Load the pinned\n"
        "module, or update GMX_VERSION_EXPECT deliberately if the change is intended.\n\n"
        % (expected, found, info.get("gmx_binary"))
    )
    sys.exit(1)


def format_block(info=None):
    """The provenance as lines, for the genetic algorithm's progress file header."""
    info = info or capture()
    lines = [
        ("GROMACS binary:\t\t\t", info.get("gmx_binary")),
        ("GROMACS version:\t\t", info.get("gmx_version")),
        ("martinize2:\t\t\t", info.get("martinize2_version")),
        ("Python:\t\t\t\t", info.get("python")),
        ("numpy:\t\t\t\t", info.get("numpy")),
        ("vermouth:\t\t\t", info.get("vermouth")),
        ("mpi4py:\t\t\t\t", info.get("mpi4py")),
    ]
    for item in info.get("gmx_banner", []):
        lines.append(("  " + item, ""))
    if "modules" in info:
        lines.append(("Modules loaded:\t\t\t", ", ".join(info["modules"])))
    return lines
