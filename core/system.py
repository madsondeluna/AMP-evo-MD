"""Miscellaneous IO and other system functions."""

import os
import sys
import shutil
import subprocess
import signal

# Guarded, matching the pattern already used in core/genetic_algorithm.py and the
# "if MPI is not None" check immediately below. Upstream this was an unconditional import, so
# the package could not even be imported without a working MPI stack: no single-process runs,
# no local testing, no scoring of finished results on a workstation.
MPI = None
try:
    from mpi4py import MPI
except ImportError:
    pass


def mprint(string, end="\n"):
    if MPI is not None:
        if MPI.COMM_WORLD.Get_rank() == 0:
            print(string,end=end)
    else:
        print(string,end=end)

class _Process:
    """Keeps track of spawned processes. Handles cleanup of processes on program termination."""
    _active_processes = set()

    @classmethod
    def _new_process(cls, process):
        cls._active_processes.add(process)

        return process

    @classmethod
    def run(cls, command, stdout, stderr, timeout=None, input=None, cwd=None):
        returncode = -1

        process = cls._new_process(subprocess.Popen(command.split(), cwd=cwd, stdin=subprocess.PIPE, stdout=stdout, stderr=stderr))
        try:
            process.communicate(input=input, timeout=timeout)
        except subprocess.TimeoutExpired:
            print("Timeout expired for process id:\t%d. Killing process." % (process.pid))
            process.kill()
            returncode = 2
        else:
            returncode = process.returncode
        cls._active_processes.discard(process)

        return returncode

    # @classmethod
    # def run_in_tmpdir(cls, command, tmpdir, stdout, stderr, timeout=None, input=None):
    #     returncode = -1

    #     process = cls._new_process(subprocess.Popen(command.split(), cwd=str(tmpdir), stdin=subprocess.PIPE, stdout=stdout, stderr=stderr))
    #     try:
    #         process.communicate(input=input, timeout=timeout)
    #     except subprocess.TimeoutExpired:
    #         print("Timeout expired for process id:\t%d. Killing process." % (process.pid))
    #         process.kill()
    #         returncode = 2
    #     else:
    #         returncode = process.returncode
    #     cls._active_processes.discard(process)

    #     return returncode

    @staticmethod
    def _getsig(signum):
        return signal.Signals(signum).name

    @classmethod
    def terminate(cls, signum, frame):
        mprint("\nSignal caught:\t%s\nTerminating simulations...\n" % (cls._getsig(signum)))
        for p in cls._active_processes:
            p.kill()
            print("killed simulation with process id:\t%d" % (p.pid))

        sys.exit(1)
    
    @classmethod
    def register(cls):
        signal.signal(signal.SIGINT, cls.terminate)
        signal.signal(signal.SIGTERM, cls.terminate)

_Process.register()

######################
## COMMON FUNCTIONS ##
######################

def replace_in_file(filepath, tag, replacement):
    """
    Replaces strings 'tag' with string 'replacement' in file 'filepath'.
    """
    data = []
    with open(filepath, 'rt') as file:
        for row in file:
            if tag in row:
                row = row.replace(tag,replacement)
            data.append(row)

    os.remove(filepath)
    
    with open(filepath, 'wt') as file:
        for row in data:
            file.write(row)

def read_text_file(filepath, strip=False):
    data = []
    if strip:
        with open(filepath, 'rt') as file:
            for line in file:
                data.append(line.strip())
    else:
        with open(filepath, 'rt') as file:
            data = file.readlines()
    return data

def write_text_file(filepath, data, add_newline=False):
    if add_newline:
        with open(filepath, 'wt') as file:
            for row in data:
                file.write(row + "\n")
    else:
        with open(filepath, 'wt') as file:
            for row in data:
                file.write(row)

def get_file_dir(file):
    return os.path.dirname(os.path.realpath(file))

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def rm_dir_content(path_to_dir):
    for filename in os.listdir(path_to_dir):
        file_path = os.path.join(path_to_dir, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

def copy_dir_files(src_dir, dest_dir):
    list_files = os.listdir(src_dir)
    for filename in list_files:
        filepath = os.path.join(src_dir, filename)
        if os.path.isfile(filepath):
            shutil.copy(filepath, dest_dir)

def _tail(path, n_lines=20):
    """Last few lines of a log, for quoting into an error message."""
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path, "rt", errors="replace") as handle:
            lines = handle.readlines()
    except OSError:
        return None
    return "".join(lines[-n_lines:]).rstrip() or None


def check_exitcode(exitcode, command=None, path_stdout=None, path_stderr=None):
    """
    Stop on a non-zero exit, naming the command and where its output went.

    This is the single most travelled error path in the pipeline: every grompp, editconf,
    energy, distance, genion and trjconv goes through it. It used to be a bare sys.exit, which
    discarded the command string that run_command had in hand. What reached the operator was
    run.py reporting that some command in some step had exited non-zero for a given candidate,
    leaving them to find the right stdout_stderr.log among a few hundred candidate directories
    and then to work out which of the fifteen or so commands in that file was the one that died.
    """
    if exitcode == 0:
        return

    if command is not None:
        detail = ["", "=" * 78, "A command exited with status %s." % exitcode,
                  "  command: %s" % command]
        for label, path in (("stdout", path_stdout), ("stderr", path_stderr)):
            if path:
                detail.append("  %s:  %s" % (label, os.path.abspath(path)))

        # GROMACS puts the reason on the last few lines; quoting them here saves the operator
        # from opening the log at all in the common cases (a missing index group, an undefined
        # atom type, a file that was never written).
        tail = _tail(path_stderr) or _tail(path_stdout)
        if tail:
            detail.append("  last lines of that output:")
            detail.extend("    " + line for line in tail.split("\n"))
        detail.append("=" * 78)
        sys.stderr.write("\n".join(detail) + "\n")
        sys.stderr.flush()

    sys.exit(exitcode)

################
## SUBPROCESS ##
################

def _try_open(filepath, mode):
    if filepath is not None:
        mgr = open(filepath, mode)
    else:
        mgr = open(os.devnull, "w")
    return mgr

def run_command(command, tTimeout=None, input=None, path_stdout=None, path_stderr=None, cwd=None):
    hdlOut = _try_open(path_stdout, 'at')
    hdlErr = _try_open(path_stderr, 'at')

    with hdlOut as out, hdlErr as err:
        returncode = _Process.run(command, stdout=out, stderr=err, timeout=tTimeout, input=input, cwd=cwd)

    check_exitcode(returncode, command=command, path_stdout=path_stdout, path_stderr=path_stderr)

# def run_command_tmpdir(command, tmpdir, tTimeout=None, input=None, path_stdout=None, path_stderr=None):
#     hdlOut = _try_open(path_stdout, 'at')
#     hdlErr = _try_open(path_stderr, 'at')

#     with hdlOut as out, hdlErr as err:
#         returncode = _Process.run_in_tmpdir(command, tmpdir, stdout=out, stderr=err, timeout=tTimeout, input=input)

#     check_exitcode(returncode)

class SimulationError(Exception):
    """ Should be raised when a GROMACS MDRUN call exits with an error."""
    pass

def _sim_check_exitcode(exitcode, command=None, path_stdout=None, path_stderr=None):
    if exitcode == 0:
        return

    # This exception used to be raised with no message at all, which is what reached the step log
    # as an empty "error" field, and what the retry loop in run.py printed as a bare "Simulation
    # stopped due to an error." An mdrun failure is almost always legible from the last lines of
    # its own log, so they are carried in the exception rather than left on disk.
    parts = ["mdrun exited with status %s" % exitcode]
    if command:
        parts.append("command: %s" % command)
    for label, path in (("stdout", path_stdout), ("stderr", path_stderr)):
        if path:
            parts.append("%s: %s" % (label, os.path.abspath(path)))
    tail = _tail(path_stderr, 12) or _tail(path_stdout, 12)
    if tail:
        parts.append("last lines of that output:\n" + tail)
    raise SimulationError("\n".join(parts))

def run_simulation_command(command, tTimeout=None, path_stdout=None, path_stderr=None):
    hdlOut = _try_open(path_stdout, 'at')
    hdlErr = _try_open(path_stderr, 'at')

    with hdlOut as out, hdlErr as err:
        returncode = _Process.run(command, stdout=out, stderr=err, timeout=tTimeout)

    _sim_check_exitcode(returncode, command=command,
                        path_stdout=path_stdout, path_stderr=path_stderr)

