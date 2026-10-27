"""
Insert the peptide into every membrane in the registry and minimise each resulting system.

Generalised from insert_peptide_2mem.py. Three things were hardcoded to the two tension
states of a single POPC membrane and are now looked up per composition:

  1. the membrane .gro, chosen by an if-chain on the literal strings 'high_tension' and
     'low_tension'
  2. the membrane .top, which was one shared file because both states had identical
     composition
  3. the insertion depth, taken as UserSettings.b[membrane_types.index(membrane_type)], i.e.
     keyed on position in a list

Everything else in the original was already generic over membrane_types and is unchanged.
"""

import os

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings
import user.membranes as Membranes


def module_insert_peptide_nmem(simulation):
    """
    Takes the peptide GRO and, for each registered membrane, builds the combined system and
    runs a steepest descent minimisation on it.
    """
    # Fail before the first grompp rather than after a few hundred nanoseconds of simulation.
    Membranes.require_all_calibrated()

    for membrane_type in UserSettings.membrane_types:
        membrane = Membranes.get(membrane_type)

        class LocalFiles():
            # Per-composition structure and topology, from the registry.
            membrane_GRO                = membrane.gro
            membrane_TOP                = membrane.top

            peptide_GRO                 = simulation._share.get_filename("gen-GRO")
            peptide_TOP                 = simulation._share.get_filename("gen-TOP")

            peptide_GRO_tmp             = "peptide_tmp_" + membrane_type + ".gro"

            TOP                         = UserSettings.name_output_system + "_" + membrane_type + ".top"
            TOPOUT                      = "insertion_out_" + membrane_type + ".top"

            GRO                         = UserSettings.name_output_system + "_" + membrane_type + ".gro"
            NDX                         = UserSettings.name_output_system + "_" + membrane_type + ".ndx"
            TPR                         = UserSettings.name_output_system + "_" + membrane_type + ".tpr"
            MDOUT                       = UserSettings.name_output_system + "_" + membrane_type + ".mdout.mdp"
            TRR                         = UserSettings.name_output_system + "_" + membrane_type + ".trr"
            EDR                         = UserSettings.name_output_system + "_" + membrane_type + ".edr"
            LOG                         = UserSettings.name_output_system + "_" + membrane_type + ".log"
            CPT                         = UserSettings.name_output_system + "_" + membrane_type + ".cpt"

            MDP = UserSettings.filename_steep_system_MDP

            stdout_stderr = None
            if CoreSettings.write_stdout_stderr_to_file:
                stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

        def path(filename):
            if filename is None:
                return None
            return os.path.join(simulation._temp_dir, filename)

        def center_peptide_to_membrane():
            # Give the peptide the membrane's box vector so the two can be merged.
            box_vector = Gromacs.GRO_read_box(path(LocalFiles.membrane_GRO))
            data_peptide_GRO = System.read_text_file(path(LocalFiles.peptide_GRO))

            for i in reversed(range(0, len(data_peptide_GRO))):
                line = data_peptide_GRO[i].strip()
                if not len(line) == 0:
                    data_peptide_GRO[i] = " ".join(box_vector)
                    break

            System.write_text_file(path(LocalFiles.peptide_GRO_tmp), data_peptide_GRO, add_newline=False)

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.EDITCONF_align_X_center(
                    path(LocalFiles.peptide_GRO_tmp),
                    path(LocalFiles.peptide_GRO_tmp)
                ),
                input="0\n".encode(),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            # Per-replica initial orientation. Roll the peptide about its own (now X-aligned)
            # helix axis by an angle unique to this replica, so independent replicas start with
            # different faces toward the bilayer rather than the same pose with only a different
            # velocity seed. rep1 gets 0 (skipped); with N_REPLICAS=3 the rolls are 0, 120, 240.
            # This is deterministic (reproducible choice of starting poses) and, combined with the
            # random velocity seed in the mdps, gives genuine replica independence.
            n_replicas = max(1, int(getattr(UserSettings, "N_REPLICAS", 1)))
            roll_deg = ((simulation.replica - 1) % n_replicas) * (360.0 / n_replicas)
            if roll_deg != 0.0:
                System.run_command(
                    core.gromacs_commands.GROMACSRunCommands.EDITCONF_rotate(
                        path(LocalFiles.peptide_GRO_tmp),
                        path(LocalFiles.peptide_GRO_tmp),
                        (roll_deg, 0, 0)
                    ),
                    path_stdout=path(LocalFiles.stdout_stderr),
                    path_stderr=path(LocalFiles.stdout_stderr)
                )

            # Starting height above the bilayer midplane. Keyed on the membrane, not on a
            # position in a list: insertion depth is a property of the bilayer thickness.
            translate_vector = (0, 0, -1 * membrane.b_nm)

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.EDITCONF_on_top_of_mem(
                    path(LocalFiles.peptide_GRO_tmp),
                    path(LocalFiles.peptide_GRO_tmp),
                    translate_vector
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        def insert_into_membrane():
            base_sections = [Gromacs.TOP_HEADER_SECTION, "system", "molecules"]
            append_sections = [Gromacs.TOP_HEADER_SECTION, "molecules"]

            Gromacs.GRO_merge(
                path(LocalFiles.membrane_GRO),
                path(LocalFiles.peptide_GRO_tmp),
                path(LocalFiles.GRO)
            )

            Gromacs.TOP_merge(
                path(LocalFiles.membrane_TOP),
                path(LocalFiles.peptide_TOP),
                path(LocalFiles.TOP),
                base_sections, append_sections
            )

            # Index groups must be rebuilt per composition. The group names are the same
            # ("Lipids", "Solvent", "Protein") but the atom indices behind them are not, which
            # is exactly why the shared prod-NDX of the original cannot survive here.
            Gromacs.NDX_from_GRO(
                path(LocalFiles.GRO),
                path(LocalFiles.NDX),
                UserSettings.resnm_to_indexgroup_dict
            )

        def steep_system():
            System.replace_in_file(
                path(LocalFiles.MDP),
                UserSettings.tags_steep_system_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX_TOPOUT(
                    path(LocalFiles.MDP),
                    path(LocalFiles.GRO),
                    path(LocalFiles.TOP),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT),
                    path(LocalFiles.NDX),
                    path(LocalFiles.TOPOUT)
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_simulation_command(
                core.gromacs_commands.GROMACSRunCommands.MDRUN(
                    path(LocalFiles.TPR),
                    path(LocalFiles.TRR),
                    path(LocalFiles.GRO),
                    path(LocalFiles.EDR),
                    path(LocalFiles.LOG),
                    path(LocalFiles.CPT),
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        System.copy_dir_files(UserSettings.path_basedir_insert_peptide, path(""))

        center_peptide_to_membrane()
        insert_into_membrane()
        steep_system()

        simulation._share.add_file("ins-%s-peptide_tmp-GRO" % membrane_type, LocalFiles.peptide_GRO_tmp, UserSettings.FILE_WRITE_DEBUGGING)
        simulation._share.add_file("ins-%s-GRO" % membrane_type,             LocalFiles.GRO,             False)
        simulation._share.add_file("ins-%s-TOP" % membrane_type,             LocalFiles.TOP,             True)
        simulation._share.add_file("ins-%s-TOPOUT" % membrane_type,          LocalFiles.TOPOUT,          UserSettings.FILE_WRITE_DEBUGGING)
        simulation._share.add_file("ins-%s-NDX" % membrane_type,             LocalFiles.NDX,             True)
