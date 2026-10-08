"""
One format for every failure an operator can act on.

The pipeline runs unattended for days on a machine the operator is not watching. When it stops,
the only thing that matters is how long it takes to find out where it stopped and why. A bare
traceback answers neither question: it names a line in a file the operator did not write, forty
frames below anything they recognise, and says nothing about which candidate, which membrane or
which step was involved.

So every failure that an operator can do something about is reported through fail() below, which
prints four things and nothing else:

    what   the thing that went wrong, in one line
    where  the file, directory, sequence, membrane or step it happened in
    why    what the pipeline expected, so the operator can tell whether the input or the
           expectation is wrong
    fix    the command or edit that resolves it

The block is delimited by rules so it is findable in a log of GROMACS output, which is the only
place it will ever be read. Tracebacks are left alone for genuine bugs: a KeyError from a typo
in this repository is not an operator error, and dressing it up as one would hide it.
"""

import os
import sys

RULE = "=" * 78
BANNER = "AMP-evo-MD: cannot continue"


class OperatorError(Exception):
    """
    A failure with a known cause and a known fix.

    Raised rather than exited in library code, so a caller holding more context (the candidate
    sequence, the membrane name) can add it before the message reaches the operator. Carries the
    formatted block so whoever prints it does not have to reassemble it.
    """

    def __init__(self, what, where=None, why=None, fix=None):
        self.what = what
        self.where = where
        self.why = why
        self.fix = fix
        super().__init__(format_block(what, where, why, fix))


def format_block(what, where=None, why=None, fix=None):
    """Render the four fields, skipping any that were not supplied."""
    lines = [RULE, BANNER, ""]
    for label, value in (("what", what), ("where", where), ("why", why), ("fix", fix)):
        if value is None:
            continue
        # Multi-line values (a list of missing files, a command) stay aligned under their label
        # rather than wrapping back to column zero, where they would read as separate messages.
        parts = str(value).rstrip().split("\n")
        lines.append("  %-6s %s" % (label + ":", parts[0]))
        for part in parts[1:]:
            lines.append("         " + part)
    lines.append(RULE)
    return "\n".join(lines)


def fail(what, where=None, why=None, fix=None, exit_code=1):
    """
    Report an operator error and stop.

    Writes to stderr and to stdout. Both, deliberately: the operator reads stdout_stderr.log,
    which captures both streams, but an interactive run under a pager or a redirect may show
    only one, and a stop reason that lands in the stream nobody is reading is the failure this
    module exists to prevent.
    """
    block = format_block(what, where, why, fix)
    sys.stderr.write(block + "\n")
    sys.stderr.flush()
    try:
        sys.stdout.write(block + "\n")
        sys.stdout.flush()
    except (OSError, ValueError):
        # A closed or broken stdout must not replace the real message with a BrokenPipeError.
        pass
    sys.exit(exit_code)


def require_file(path, why, fix, what=None):
    """
    Stop unless path is an existing file, naming where it was looked for.

    A relative path is reported with the absolute path it resolved to and the working directory
    it resolved against, because the most common form of this failure is running from the wrong
    directory, and "No such file: seq.txt" does not reveal that.
    """
    if os.path.isfile(path):
        return path

    where = os.path.abspath(path)
    if not os.path.isabs(path):
        where += "\n(resolved from %r relative to the current directory %s)" % (path, os.getcwd())
    elif os.path.isdir(path):
        where += "\n(this path exists but is a directory, not a file)"

    fail(what or ("%s not found" % os.path.basename(path)), where=where, why=why, fix=fix)
