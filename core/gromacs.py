"""Provides various functions related to GROMACS. NOTE: Functions are only tested for specific applications. """
#TODO cleanup

from    math        import floor
from    contextlib  import suppress
from    bisect      import bisect_left, bisect_right
import  functools
import os
import sys

import  core.system  as system
#import  user.seq2itp_aminoacids as seq2itp

TOP_HEADER_SECTION          = "header"
TOP_SYSTEM_SECTION          = "system"
TOP_MOLECULES_SECTION       = "molecules"
ITP_MOLECULE_TYPE_SECTION   = "moleculetype"
ITP_ATOMS_SECTION           = "atoms"

class GROFileData:
    def __init__(self, header, atom_count, data, box_vector):
        self.header     = header
        self.atom_count = atom_count
        self.data       = data
        self.box_vector = box_vector

def TOP_get_system_charge(filepath_TOP, TOP_itp_directory_list=[]):
    def ITP_get_molecule_charges(filepath_ITP):
        class Molecule:
            def __init__(self, name):
                self.name = name
                self.charge = 0

        data = {}
        current_molecule = None

        def commit():
            # Record whatever molecule is open. Called on EVERY section header and once at end
            # of file, which is the whole correction: the original committed a molecule only
            # when a section that was neither moleculetype nor atoms came along, so a molecule
            # followed directly by the next [ moleculetype ] was dropped and silently
            # contributed zero charge. That is the shape of every single-bead species, which
            # means the ions: of the 18 moleculetypes in martini_v3.0.0_ions_v1.itp only 6 were
            # parsed, and the 12 lost included NA and CL, the counterions this function exists
            # to count. A system built with insane's own salt then computed as net charged, and
            # GRO_TOP_neutralize_system replaced solvent beads with ions it did not need,
            # leaving every production box charged under PME with nothing in the output to say
            # so.
            if current_molecule is not None:
                data[current_molecule.name] = current_molecule.charge

        with open(filepath_ITP, 'rt') as file:
            state = None
            for number, line in enumerate(file, start=1):
                line = line.split(";", 1)[0].strip()
                if not line:
                    continue

                if line.startswith("["):
                    if ITP_MOLECULE_TYPE_SECTION in line:
                        # A new molecule begins, so the previous one is complete. Its name
                        # arrives on the next data line.
                        commit()
                        current_molecule = None
                        state = ITP_MOLECULE_TYPE_SECTION
                    elif ITP_ATOMS_SECTION in line:
                        state = ITP_ATOMS_SECTION
                    else:
                        # [ bonds ], [ angles ] and the rest: stop reading atoms, but the
                        # molecule is still open and must not be dropped.
                        state = None
                    continue

                if state == ITP_MOLECULE_TYPE_SECTION:
                    current_molecule = Molecule(line.split()[0])
                elif state == ITP_ATOMS_SECTION:
                    if current_molecule is None:
                        raise ValueError(
                            "%s line %d: an [ atoms ] section appears before any "
                            "[ moleculetype ]" % (filepath_ITP, number))
                    fields = line.split()
                    if len(fields) < 7:
                        raise ValueError(
                            "%s line %d: an atom line needs at least 7 columns, the 7th being "
                            "the charge; found %d: %s"
                            % (filepath_ITP, number, len(fields), line))
                    current_molecule.charge += round(float(fields[6]))

            commit()
        return data

    TOP_data = TOP_read(filepath_TOP)

    header_data = TOP_data[TOP_HEADER_SECTION]
    molecule_charge_data = {}
    TOP_itp_directory = ""
    for line in header_data:
        line = line.strip()
        if "#include" in line:
            itp_filename = line.split("\"")[1].split("\"")[0]
            if len(TOP_itp_directory_list) > 0:
                for dirname in TOP_itp_directory_list:
                    full_path = os.path.join(dirname, itp_filename)
                    if os.path.isfile(full_path):
                        TOP_itp_directory = dirname
                        break
            molecule_charge_data = {**ITP_get_molecule_charges(os.path.join(TOP_itp_directory, itp_filename)), **molecule_charge_data}

    total_charge_system = 0
    unknown = []
    for line in TOP_data[TOP_MOLECULES_SECTION]:
        line = line.strip().split()
        if not line:
            continue
        if line[0] in molecule_charge_data:
            total_charge_system +=  molecule_charge_data[line[0]] * round(float(line[1]))
        else:
            # Silently contributing zero is how a charged system comes out looking neutral.
            unknown.append("%s x %s" % (line[0], line[1]))

    if unknown:
        # A system that computes as neutral while some of its molecules were never costed is the
        # dangerous case: no counterions are added, and grompp's "System has non-zero total
        # charge" note is swallowed by the max_warnings allowance in production.py. The run then
        # produces a perfectly plausible surface tension for a charged periodic box.
        message = ("could not determine the charge of: %s (no matching moleculetype in the "
                   "#include'd itp files of %s)" % ("; ".join(unknown), filepath_TOP))
        if total_charge_system == 0:
            import core.errors as Errors
            raise Errors.OperatorError(
                "the system charge cannot be verified",
                where=message,
                why=("the total came out as zero, but molecules whose charge is unknown "
                     "contribute zero to that total, so a charged system is indistinguishable "
                     "from a neutral one here and would run without counterions"),
                fix=("check that every moleculetype in [ molecules ] is defined in one of the "
                     "#include'd itp files, and that the membrane .top was copied together with "
                     "the itps it includes"))
        sys.stderr.write("warning: %s; the total charge %+d is computed without them\n"
                         % (message, total_charge_system))

    return total_charge_system

def GRO_TOP_neutralize_system(filepath_input_GRO, filepath_input_TOP, filepath_output_GRO, filepath_output_TOP, name_Solvent, name_positive_ion, name_negative_ion, TOP_itp_directory_list=[]):
    system_charge = TOP_get_system_charge(filepath_input_TOP, TOP_itp_directory_list=TOP_itp_directory_list)
    data_TOP = TOP_read(filepath_input_TOP)
    data_GRO = GRO_read(filepath_input_GRO)

    if system_charge == 0:
        pass # No need to change the files
    else:
        # The solvent line is what pays for the counterions, in both the topology and the
        # structure. If no molecule in [ molecules ] is called name_Solvent, the loop below adds
        # nothing and name_Ion is never assigned, and the failure arrives at the GRO loop as
        # "UnboundLocalError: local variable 'name_Ion' referenced before assignment", forty
        # frames below anything the operator wrote. The usual cause is a membrane built with a
        # solvent name other than UserSettings.name_solvent_GENION: W versus WF versus SOL.
        names = [line.strip().split()[0] for line in data_TOP[TOP_MOLECULES_SECTION]
                 if line.strip()]
        if name_Solvent not in names:
            import core.errors as Errors
            raise Errors.OperatorError(
                "system has a net charge of %+d but there is no %r molecule to take the "
                "counterions from" % (system_charge, name_Solvent),
                where="%s\n[ molecules ] contains: %s" % (filepath_input_TOP, ", ".join(names)),
                why=("the counterions replace solvent beads, so the solvent has to be findable "
                     "by name; UserSettings.name_solvent_GENION is currently %r" % name_Solvent),
                fix=("set name_solvent_GENION in user/usersettings.py to the name this topology "
                     "uses, or rebuild the membrane with that solvent name"))

        #TOP#
        output_data_TOP = []
        molecules = data_TOP[TOP_MOLECULES_SECTION]
        for i in range(len(molecules)):
            molecule=molecules[i].strip().split()

            if molecule[0] == name_Solvent:
                if system_charge > 0:
                    name_Ion = name_negative_ion
                elif system_charge < 0:
                    name_Ion = name_positive_ion
                output_data_TOP.append("" + name_Ion + " " + str(abs(system_charge)) + "\n")

                new_count = int(molecule[1]) - abs(system_charge)
                output_data_TOP.append("" + name_Solvent + " " + str(new_count) + "\n")
            else:
                output_data_TOP.append(molecules[i])

        data_TOP[TOP_MOLECULES_SECTION] = output_data_TOP

        #GRO#
        count_ION = abs(system_charge)

        for i in range(len(data_GRO.data)):
            if count_ION <= 0:
                break

            if name_Solvent in data_GRO.data[i][5:10].strip():
                data_GRO.data[i] = data_GRO.data[i][0:5] + "%-5s%5s" % (name_Ion, name_Ion) + data_GRO.data[i][15:]
                count_ION -= 1
    
    with open(filepath_output_TOP, 'wt') as file:
        for line in data_TOP[TOP_HEADER_SECTION]:
            file.write(line)
        file.write("\n")
        file.write("[" + TOP_SYSTEM_SECTION + "]\n")
        for line in data_TOP[TOP_SYSTEM_SECTION]:
            file.write(line)
        file.write("\n")
        file.write("[" + TOP_MOLECULES_SECTION + "]\n")
        for line in data_TOP[TOP_MOLECULES_SECTION]:
            file.write(line)

    with open(filepath_output_GRO, 'wt') as file:
        for line in data_GRO.header:
            file.write(line)
        file.write(str(data_GRO.atom_count) + "\n")
        for line in data_GRO.data:
            file.write(line)
        file.write(" ".join(data_GRO.box_vector))

def GRO_read_box(filepath):
    """
    Reads the box vector from GRO file.

    Args:
        filepath (string): Path to GRO file.
    Returns:
        (list): List of strings describing the box vector [v1(x), v2(y), v3(z), v1(y), v1(z), v2(x), v2(z), v3(x), v3(y)].
    """
    box_string = None
    data = system.read_text_file(filepath)
    
    for line in reversed(data):
        line = line.strip()
        if not len(line) == 0:
            box_string = line.split()
            break

    return box_string

def GRO_read(filepath_GRO):
    """
    Args:
        filepath_GRO (string):  Path to GRO file.
    Returns:
        GROFileData object, with members:
            header (string):    Title string of the file.
            atom_count (int):   Number of atoms in file.
            data (list):        List of strings, each string describing an atom. Strings contain line endings.
            box_vector (list):  List of strings describing the box vector.
    """
    file_data = system.read_text_file(filepath_GRO, strip=False)

    header      = file_data[0]
    atom_count  = int(file_data[1])
    data        = [file_data[i] for i in range(2, 2 + atom_count)]
    box_vector  = file_data[2 + atom_count].strip().split()

    return GROFileData(header, atom_count, data, box_vector)

def GRO_merge(filepath_GRO_base, filepath_GRO_insert, filepath_GRO_output):
    """
    Combines two GRO files (base + insert) into one, keeping the header and box vector of filepath_GRO_base.
    """
    data_base   = system.read_text_file(filepath_GRO_base)
    data_insert = system.read_text_file(filepath_GRO_insert)
    
    # Strip box from data_base
    for i in reversed(range(0, len(data_base))):
        line = data_base[i].strip()
        if not len(line) == 0:
            del(data_base[i])
            break
        del(data_base[i])
    
    data_output     = data_base + data_insert[2:]
    data_output[1]  = str(len(data_output) - 3) + "\n"
    system.write_text_file(filepath_GRO_output, data_output, add_newline=False)

def TOP_read(filepath):
    """
    Returns:
        (dict): Each key represents a section of the TOP/ITP file (e.g. [ atoms ] becomes key "atoms").
                The corresponding value is a list of strings containing the lines belonging to that section.
    """
    data_sections = {}
    with open(filepath, 'rt') as file:
        state = TOP_HEADER_SECTION
        data_sections[state] = []
        for row in file:
            row_stripped = row.strip()
            if len(row_stripped) == 0 or row_stripped.startswith(";"):
                continue
            if row_stripped.startswith("["):
                state = row_stripped.split("[", 1)[-1].split("]")[0].strip()
                if state not in data_sections:
                    data_sections[state] = []
                continue

            data_sections[state].append(row)

    return data_sections

def TOP_merge(filepath_TOP_base, filepath_TOP_insert, filepath_TOP_output, sections_list, sections_list_append):
    """
    Combines two TOP files into one (base + insert).
    
    Args:
        sections_list (list): List of strings describing which sections of 'filepath_TOP_base' will be written to the output file.
        sections_list_append (list): List of strings describing which sections of 'filepath_TOP_insert' are appended 
                                        to the sections of 'filepath_TOP_base' in the output file.
                                        Only sections also present in 'sections_list' are included.
    """
    data_sections_base      = TOP_read(filepath_TOP_base)
    data_sections_insert    = TOP_read(filepath_TOP_insert)
    header_union = None

    if TOP_HEADER_SECTION in sections_list and TOP_HEADER_SECTION in sections_list_append:
        header_combined = data_sections_base[TOP_HEADER_SECTION] + data_sections_insert[TOP_HEADER_SECTION]
        for i in range(len(header_combined)):
            include_filename    =  header_combined[i].split("\"", 1)[-1].split("\"")[0]
            header_combined[i]  = "#include \"" + include_filename + "\"\n"
        header_union = list(set(header_combined)) # Removes duplicates. Removes original order of list as well.

        #Ensure headers are in order in which they appear in header_combined
        i_sorted = []
        for i in range(len(header_union)):
            i_sorted.append(header_combined.index(header_union[i]))
        
        i_sorted = sorted(i_sorted)
        for i in range(len(i_sorted)):
            header_union[i] = header_combined[i_sorted[i]]

    with open(filepath_TOP_output, 'wt') as file:
        for section in sections_list:
            if section == TOP_HEADER_SECTION and header_union is not None:
                for line in header_union:
                    file.write(line)
                continue
            elif not section == TOP_HEADER_SECTION:
                file.write("\n[" + section + "]\n")
            else:
                continue
            for row in data_sections_base[section]:
                file.write(row)
            if section in sections_list_append:
                for row in data_sections_insert[section]:
                    file.write(row)

# def ITP_convert_to_GRO(filepath_ITP, filepath_GRO):
#     """
#     Similar to old_ITP_convert_to_GRO, except that atoms are placed according to hard coded amino acid structures. 
#     This helps prevent possible complications during steeping. Otherwise same in- and output as old_ITP_convert_to_GRO.
#     """

#     data_ITP = TOP_read(filepath_ITP)["atoms"]
#     seq2itp.seq2itp_ITP_to_GRO(data_ITP, filepath_GRO)

def ITP_constraints_to_bonds(filepath_ITP, filepath_ITP_noConstraints, constraint_forceConstant=5000):
    """
    Creates a copy of 'filepath_ITP' where constraints have been replaced with harmonic bonds.
    NOTE: Only tested on ITP files containing a single molecule.

    Args:
        filepath_ITP (string): Path to input ITP file.
        filepath_ITP_noConstraints (string): Path to output ITP file.
        constraint_forceConstant (int): Force constant of the new harmonic bonds.
    """
    data = []
    constraints = []
    with open(filepath_ITP, 'rt') as file:
        state_constraints = False
        for row in file:
            row = row.strip()
            if(row.startswith("[")):
                if "constraints" in row:
                    state_constraints = True
                    continue
                else:
                    state_constraints = False
            if state_constraints:
                constraints.append(row.split(";")[0])
            else:
                data.append(row)
    
    with open(filepath_ITP_noConstraints, 'wt') as file:
        for row in data:
            file.write(row + "\n")
            if(row.startswith("[") and "bonds" in row):
                # The newline was missing, so when the constraints block does not open with a
                # comment line the first converted bond was appended to this comment and
                # silently commented out: one bond short in every minimisation topology, and
                # nothing says so. Masked today only because the shipped itps happen to start
                # that section with a comment.
                file.write("; Constraints converted to bonds for energy minimization purposes\n")
                for row_constraints in constraints:
                    if not (len(row_constraints) == 0 or row_constraints.startswith(";")):
                        # str(): the default for this parameter is an int, so the default was
                        # unusable (TypeError on the concatenation). Production passes a
                        # string, which is why it never surfaced.
                        row_constraints = row_constraints + " " + str(constraint_forceConstant)
                    file.write(row_constraints + "\n")
                continue

def NDX_write(filepath_NDX, index_dict):
    """
    Args:
        filepath_NDX (string): Path to output NDX file.
        index_dict (dict(string, int)): Dictionary with keys describing the index group names 
                                        and values describing which atoms belong to that group.
    """
    with open(filepath_NDX, 'wt') as file:
        for group, indices in index_dict.items():
            count = 0
            file.write("[ " + group + " ]\n")
            for index in indices:
                file.write(str(index) + "\t")
                count+=1
                if count >= 15:
                    file.write("\n")
                    count = 0
            file.write("\n")

def NDX_read(filepath_NDX):
    """Returns index_dict, see NDX_write"""

    NDX_dict = {}
    state = None
    with open(filepath_NDX, 'rt') as file:
        for line in file:
            line = line.strip()
            if len(line) == 0:
                continue

            if line.startswith("["):
                state = line.split("[")[1].split("]")[0].strip()
                if state not in NDX_dict:
                    NDX_dict[state] = []
                continue
            
            if state is None:
                continue
            else:
                indices = line.split()
                for index in indices:
                    NDX_dict[state].append(int(index))
                continue
    
    return NDX_dict

def NDX_from_GRO(filepath_GRO, filepath_NDX_out, resnm_to_group_dict):
    """
    Produces an NDX file from a GRO file according to a reference dictionary.
    
    Args:
        resnm_to_group_dict (dict(string, string/tuple)): Dict with keys describing residue names, 
                                                    and values describing to which index group they belong.
    """
    data        = GRO_read(filepath_GRO).data
    NDX_dict    = {}

    for i in range(len(data)):
        resnm = data[i][5:10].strip()
        try:
            groups = resnm_to_group_dict[resnm]
        except KeyError:
            # This fork exists to run compositions the original never did, so an unlisted
            # residue is an expected operator error rather than a bug. A bare KeyError here
            # names neither the file, nor the line, nor the setting that needs editing.
            import core.errors as Errors
            raise Errors.OperatorError(
                "residue name %r is not assigned to an index group" % resnm,
                where="%s, atom line %d:\n%s" % (filepath_GRO, i + 1, data[i].rstrip()),
                why=("every residue in the structure has to be classified so that the index "
                     "groups the .mdp refers to can be built.\n"
                     "Currently known: %s" % ", ".join(sorted(resnm_to_group_dict))),
                fix=("add %r to UserSettings.resnm_to_indexgroup_dict with its group "
                     "(Lipid, Solvent or Protein). Note that Martini 2 and Martini 3 differ "
                     "on ion names: NA+/CL- against NA/CL." % resnm))

        if not isinstance(groups, tuple):
            groups = (groups,)

        for group in groups:
            if group not in NDX_dict:
                NDX_dict[group] = []
            NDX_dict[group].append(str(i+1)) # atom index counts from 1 instead of 0

    NDX_write(filepath_NDX_out, NDX_dict)

def NDX_split_group(filepath_NDX_input, filepath_NDX_output, resname, atoms_per_residue, resname_split):
    """
    Split an index group into multiple smaller groups, keeping the original group intact.

    Args:
        resname (string): Index group to split
        atoms_per_residue (int): Number of indices that go into each new group
        resname_split (string): Name for the new groups, this name will be appended by an index number (e.g. "CL" -> "CL1" "CL2" "CL3" ...)
    """
    NDX_dict = NDX_read(filepath_NDX_input)

    if resname not in NDX_dict:
        raise KeyError("Error: resname (" + resname + ") not found in index file (" + filepath_NDX_input + ")")

    indices = NDX_dict[resname]

    count = 0
    index_new_residue = 1
    new_residue = []
    for index in indices:
        new_residue.append(index)
        count += 1

        if count >= atoms_per_residue:
            NDX_dict[resname_split + str(index_new_residue)] = new_residue

            index_new_residue += 1
            count = 0
            new_residue = []
    
    NDX_write(filepath_NDX_output, NDX_dict)

def GMX_average_ENERGY_output(filepath_GMX_ENERGY_output):
    # Returns dict containing energy selections and average

    output_column_template = ["Energy", "Average", "Err.Est."]

    with open(filepath_GMX_ENERGY_output, 'rt') as file:
        state = False
        for line in file:
            line = line.strip()

            if state:
                if line.startswith('-'):
                    continue
                columns = line.split()
                return float(columns[1]), float(columns[2])
                continue

            if line.startswith("Energy"):
                columns = line.split()
                if columns[0:3] == output_column_template:
                    state = True
                    continue

def XVG_read(filepath_XVG):
    """Returns a dictionary containing:

        "x":        list of x data points

        "y":        list of y data point lists,
                    first index indicates data column

        "title":    string

        "legend":   list of strings matching the data columns in the file

        "xlabel":   string

        "ylabel":   string
    """
    XVG_dict    = {}

    x           = []
    ylist       = None
    legend      = []
    title       = ""
    xlabel      = ""
    ylabel      = ""

    with open(filepath_XVG, "rt") as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("@"):
                splt = line.split()
                if splt[1] == "title":
                    title = " ".join(splt[2:])
                elif splt[1] == "xaxis":
                    xlabel = " ".join(splt[3:])
                elif splt[1] == "yaxis":
                    ylabel = " ".join(splt[3:])
                elif len(splt) > 3 and splt[2] == "legend":
                    legend.append(" ".join(splt[3:]))
                continue

            data = line.split()
            if ylist is None:
                ylist = [[] for i in range(len(data)-1)]

            x.append(float(data[0]))
            data = data[1:]
            for i, value in enumerate(data):
                ylist[i].append(float(value))

    XVG_dict["x"]       = x
    XVG_dict["y"]       = ylist
    XVG_dict["title"]   = title
    XVG_dict["legend"]  = legend
    XVG_dict["xlabel"]  = xlabel
    XVG_dict["ylabel"]  = ylabel

    return XVG_dict

def XVG_ymax(filepath_XVG, x_range=None, column_index=0):
    XVG_dict = XVG_read(filepath_XVG)
    x_data = XVG_dict["x"]
    y_data = XVG_dict["y"][column_index]

    if x_range:
        idx_range = (bisect_left(x_data, x_range[0]), bisect_right(x_data, x_range[1]))
        y_data = y_data[idx_range[0]:idx_range[1]]

    ymax = y_data[0]
    for value in y_data:
        if value > ymax:
            ymax = value

    return ymax

