from tempfile import TemporaryDirectory
import shutil
import os
import sys

import core.steplog as StepLog


class Module():
    def __init__(self, module_pointer, name=None):
        self._ptr = module_pointer  # ptr to module function
        # Falls back to the function name, so existing callers that pass only the pointer still
        # produce a readable step name in the log.
        self.name = name or getattr(module_pointer, "__name__", "module")

    def run(self, simulation):
        return self._ptr(simulation)

class Simulation():
    class SharedFiles():  # Contains filenames in the temporary directory that need to be shared between modules/output folder.
        class File():
            def __init__(self, filename, bWrite):
                self.filename = filename
                self.bWrite = bWrite

        def __init__(self, src_dir):
            self.dir = src_dir
            self._files = dict()

        def add_file(self, str_ID, filename, bWrite):
            if str_ID not in self._files:
                self._files[str_ID] = self.File(filename, bWrite)
            else:
                raise ValueError("str_ID already present in the SharedFiles object. Check modules for overlapping str_ID's.")

        def get_filename(self, str_ID):
            return self._files[str_ID].filename

        def write_to_output_dir(self, dst_dir):
            for file_obj in self._files.values():
                if file_obj.bWrite:
                    file_obj.bWrite = False
                    shutil.copyfile(os.path.join(self.dir, file_obj.filename), os.path.join(dst_dir, file_obj.filename))

    def __init__(self, output_dir, sequence=None, keep_on_failure=True,
                 replica=1, retain_all=False):
        self._input_data_paths  = []            # List of paths to input files, will be copied to temp dir upon start of run()
        self._output_dir        = output_dir    # To which certain output (specified by the called module) will be copied from the temp_dir
        self._modules           = []            # Contains module objects in the order of execution
        self._output_variables  = dict()        # Contains variable output from modules

        self._share             = None          # SharedFiles object accessed by the modules
        self._temp_dir          = None          # In which all simulation input/output will be stored, created/destroyed outside the class

        self._sequence          = sequence      # Recorded in the step log, so a log line identifies its candidate
        # Replica index (1-based). The modules read it to vary the initial peptide orientation
        # per replica, so independent replicas start from genuinely different states rather than
        # only different velocity seeds. Default 1 leaves the single-run and GA paths unchanged.
        self.replica            = replica
        # retain_all copies the whole working directory next to the output on SUCCESS, not only
        # on failure, so every intermediate trajectory (minimisation, relax, equilibration,
        # production) is kept for inspection. Disk-heavy: for validation runs, not the campaign.
        self._retain_all        = retain_all
        # On failure, copy the temporary directory next to the output rather than letting
        # TemporaryDirectory delete it. Without this the intermediates of a failed candidate are
        # gone by the time anyone looks, and only the concatenated stdout survives, which is the
        # worst moment to have the least information.
        self._keep_on_failure   = keep_on_failure
        self._log               = StepLog.StepLog(output_dir, sequence)

    def add_input_file(self, filepath):
        if os.path.isfile(filepath):
            self._input_data_paths.append(filepath)
        else:
            print("ERROR: filepath " + str(filepath) + " does not point to a file.\nExiting.\n")
            sys.exit(1)

    def add_module(self, module_pointer):
        self._modules.append(module_pointer)

    def _mark_log(self, step_name, index):
        """
        Write a banner into the concatenated command log before each step starts.

        Every GROMACS and martinize2 call appends its stdout and stderr to one file. Without a
        separator that file is an undifferentiated wall of output, and finding which command
        produced a given error means guessing. These banners are what make it greppable.
        """
        import core.settings as CoreSettings
        if not CoreSettings.write_stdout_stderr_to_file or not self._output_dir:
            return
        path = os.path.join(self._output_dir, CoreSettings.filename_output_stdout_stderr)
        try:
            os.makedirs(self._output_dir, exist_ok=True)
            with open(path, "a") as handle:
                handle.write("\n%s\n=== STEP %d: %s ===\n%s\n"
                             % ("=" * 78, index, step_name, "=" * 78))
        except OSError:
            # Same rule as the step log: observing the run must never break it.
            pass

    def _preserve_temp_dir(self, step_name):
        """Copy the working directory next to the output so a failure can be diagnosed."""
        if not self._keep_on_failure or not self._temp_dir:
            return None
        destination = os.path.join(self._output_dir, "failed_" + str(step_name))
        try:
            if os.path.exists(destination):
                shutil.rmtree(destination)
            shutil.copytree(self._temp_dir, destination)
            return destination
        except OSError:
            # Preserving the evidence must never become the reason the run stops.
            return None

    def _retain_temp_dir(self):
        """Copy the whole working directory into the output on success (retain_all)."""
        if not self._retain_all or not self._temp_dir:
            return None
        destination = os.path.join(self._output_dir, "work")
        try:
            if os.path.exists(destination):
                shutil.rmtree(destination)
            shutil.copytree(self._temp_dir, destination)
            return destination
        except OSError:
            # Keeping the intermediates must never become the reason the run stops.
            return None

    def run(self):
        with TemporaryDirectory() as self._temp_dir:
            for path in self._input_data_paths:
                filename = os.path.basename(path)
                shutil.copyfile(path, os.path.join(self._temp_dir, filename))

            self._share = self.SharedFiles(self._temp_dir)
            self._log.event("start", modules=[m.name for m in self._modules],
                            output_dir=self._output_dir)

            for index, module in enumerate(self._modules):
                self._mark_log(module.name, index)
                try:
                    with self._log.step(module.name, index=index):
                        module.run(self)
                        self._share.write_to_output_dir(self._output_dir)
                except BaseException as exc:
                    # BaseException on purpose: core.system.check_exitcode raises SystemExit on a
                    # non-zero command, and that is exactly the case where the intermediates and
                    # a log line are most needed. The exception is re-raised untouched.
                    kept = self._preserve_temp_dir(module.name)
                    self._log.event("abort", step=module.name, index=index,
                                    error_type=type(exc).__name__, error=str(exc)[:2000],
                                    preserved=kept)
                    raise

            retained = self._retain_temp_dir()
            self._log.event("finish", status="ok", retained=retained)

    def get_output_variable(self, key):
        try:
            return self._output_variables[key]
        except KeyError:
            # Reached when a module returned without publishing what the caller expects, which
            # happens on a partial edit to the module chain. A bare KeyError: 'fitness' names
            # neither the candidate nor the module that was supposed to set it.
            raise KeyError(
                "no output variable %r was produced for sequence %s in %s. Published: %s. "
                "The module that sets it is the last one in the chain; check that it ran and "
                "that it assigned simulation._output_variables[%r]."
                % (key, self._sequence or "(unknown)", self._output_dir,
                   ", ".join(sorted(self._output_variables)) or "nothing", key))

