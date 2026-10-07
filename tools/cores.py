#!/usr/bin/env python3
"""
Work out how many MPI ranks and threads per rank to use on this machine.

    eval "$(python tools/cores.py --export)"     # sets AMP_NP and OMP_NUM_THREADS
    python tools/cores.py                        # explain the choice

WHY NOT SIMPLY nproc MINUS ONE
------------------------------
Two reasons, and both change the answer.

`nproc` counts logical CPUs, so on a machine with hyperthreading it returns twice the number of
physical cores. GROMACS treats the extra logical CPUs as real ones and will happily place two
threads on one physical core, where they contend for the same execution units. For the small
coarse-grained systems here that is usually a net loss, so this script counts physical cores by
default.

More importantly, there is no single number to set. Total occupancy is ranks times threads per
rank, and the split between them matters more than the product. This pipeline evaluates many
independent candidates, so what matters is throughput, not how fast any one simulation finishes.
A roughly 5000-particle Martini system stops scaling after a few threads because domain
decomposition and communication overhead grow faster than the work does. Eight ranks of two
threads therefore finishes eight candidates in the time two ranks of eight threads finishes two,
even though both occupy sixteen cores.

So: many ranks, few threads each. This script defaults to two threads per rank and spends the
rest on ranks.

RESERVING A CORE
----------------
One core is left free by default. On a shared machine that keeps it responsive for whoever else
is logged in, and on your own machine it keeps the session usable while a campaign runs. It
costs a few percent of throughput and buys the ability to type.
"""

import argparse
import os
import subprocess
import sys


def logical_cpus():
    """Logical CPUs this process may actually use, respecting affinity and cgroup limits."""
    try:
        # Honours taskset and cgroup CPU pinning, unlike os.cpu_count().
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def cgroup_quota():
    """CPU quota from a cgroup, in whole CPUs, or None when unlimited or absent."""
    # cgroup v2
    try:
        with open("/sys/fs/cgroup/cpu.max") as handle:
            quota, period = handle.read().split()
        if quota != "max":
            return max(1, int(float(quota) / float(period)))
    except (OSError, ValueError):
        pass
    # cgroup v1
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as q, \
             open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as p:
            quota, period = int(q.read().strip()), int(p.read().strip())
        if quota > 0:
            return max(1, quota // period)
    except (OSError, ValueError):
        pass
    return None


def physical_cores():
    """
    Physical cores, from lscpu. Returns None when it cannot be determined.

    sockets times cores-per-socket, which excludes hyperthreads. Parsed from the machine-readable
    lscpu output so the result does not depend on the locale of the display format.
    """
    try:
        out = subprocess.run(["lscpu", "-p=CORE,SOCKET"], capture_output=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    pairs = set()
    for line in out.stdout.decode("utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) >= 2:
            pairs.add((parts[0], parts[1]))
    return len(pairs) or None


def plan(reserve=1, threads_per_rank=2, use_logical=False, budget=None):
    """Return a dict describing the recommended layout and how it was reached."""
    logical = logical_cpus()
    physical = physical_cores()
    quota = cgroup_quota()

    if budget is not None:
        available = budget
        basis = "requested with --cores"
    elif use_logical or physical is None:
        available = logical
        basis = "logical CPUs" if not use_logical else "logical CPUs (--logical)"
    else:
        # Affinity may restrict us below the machine's physical core count.
        available = min(physical, logical)
        basis = "physical cores"

    if quota is not None:
        available = min(available, quota)
        basis += ", capped by cgroup quota"

    usable = max(1, available - reserve)
    threads = max(1, min(threads_per_rank, usable))
    ranks = max(1, usable // threads)

    return {
        "logical": logical,
        "physical": physical,
        "quota": quota,
        "basis": basis,
        "available": available,
        "reserved": reserve,
        "usable": usable,
        "ranks": ranks,
        "threads": threads,
        "occupancy": ranks * threads,
        "idle": usable - ranks * threads,
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reserve", type=int, default=1,
                        help="cores to leave free for the rest of the machine (default 1)")
    parser.add_argument("--threads-per-rank", type=int, default=2,
                        help="OMP threads per MPI rank (default 2; see the note on scaling)")
    parser.add_argument("--logical", action="store_true",
                        help="count logical CPUs including hyperthreads, instead of physical cores")
    parser.add_argument("--cores", type=int, default=None,
                        help="override the detected core count entirely")
    parser.add_argument("--export", action="store_true",
                        help="print shell assignments for eval, nothing else")
    args = parser.parse_args()

    p = plan(reserve=args.reserve, threads_per_rank=args.threads_per_rank,
             use_logical=args.logical, budget=args.cores)

    if args.export:
        print("export AMP_NP=%d" % p["ranks"])
        print("export OMP_NUM_THREADS=%d" % p["threads"])
        return

    print("Detected")
    print("  logical CPUs        %d" % p["logical"])
    print("  physical cores      %s" % (p["physical"] if p["physical"] else "could not determine"))
    if p["quota"]:
        print("  cgroup quota        %d CPUs" % p["quota"])
    print("  counting            %s" % p["basis"])
    print()
    print("Plan")
    print("  available           %d" % p["available"])
    print("  reserved            %d" % p["reserved"])
    print("  usable              %d" % p["usable"])
    print("  MPI ranks           %d" % p["ranks"])
    print("  threads per rank    %d" % p["threads"])
    print("  total occupancy     %d%s" % (p["occupancy"],
          ("  (%d core(s) idle from the split)" % p["idle"]) if p["idle"] else ""))
    print()
    print("Run with")
    print("  export OMP_NUM_THREADS=%d" % p["threads"])
    print("  mpirun -np %d python run.py --output_dir <dir>" % p["ranks"])
    print()
    if p["ranks"] == 1:
        print("  Only one rank fits. Without MPI, drop mpirun and raise OMP_NUM_THREADS to")
        print("  %d, since a single GROMACS process then has the machine to itself." % p["usable"])
    if p["idle"]:
        print("  %d core(s) are idle because the usable count does not divide by the threads per"
              % p["idle"])
        print("  rank. Either accept that, or pass --threads-per-rank 1 to use them.")


if __name__ == "__main__":
    main()
