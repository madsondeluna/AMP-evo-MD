import os
import sys
import core.system as System
import getopt
import subprocess

##############
## Settings ##
##############

# PDB generator #
pdb_gen                             = "PeptideBuilder" # "Modeller"

# Force field #
ff                                  = "martini3001"

# Input #
filename_seq                        = "input.seq"
filename_ssd                        = "input.ssd"

# Output #
name_output_score                   = "fitness"
name_output_peptide                 = "peptide"
name_output_system                  = "system"
name_output_production              = "production"

# Simulation #
timeout_generate_peptide_MDRUN      = 60*2 # In seconds
timeout_Production_relax_MDRUN      = 600 # In seconds (10 min; relax is short, generous under load)
timeout_Production_MDRUN            = 60*60*4 # In seconds (4 h; ~1.2 h/200 ns at 4 threads, margin for contention)

# Basefile #
path_basedir_generate_peptide       = System.get_file_dir(__file__) + "/base/generate_peptide/"
path_basedir_insert_peptide         = System.get_file_dir(__file__) + "/base/insert_peptide/"
path_basedir_production             = System.get_file_dir(__file__) + "/base/production/"

# generate peptide
filename_steep_peptide_TOP          = "steep_peptide.top"
if pdb_gen == "Modeller":    
    template_pdb                        = "template.pdb"

if ff == "martini22":
    template_top                        = "template_m22.top"
elif ff == "martini3001":
    template_top                        = "template_m3001.top"
else:
    # Printing and continuing left template_top undefined, and the failure then arrived much
    # later and somewhere else, as "AttributeError: module 'user.usersettings' has no attribute
    # 'template_top'" from inside the peptide build of the first candidate.
    raise SystemExit(
        "user/usersettings.py: ff = %r is not supported. Expected 'martini22' or "
        "'martini3001'. This selects the topology template the peptide build starts from, so "
        "nothing can run without it." % (ff,))
#filename_NVT_peptide_MDP            = "nvt_peptide.mdp"
filename_MDP_steep_pept             = "steep_pept.mdp"
filename_MDP_eq_pept                = "eq_small_dt_pept.mdp"

# insert peptide
#
# membrane_types used to be the two tension states of one POPC membrane. It is now the list
# of compositions defined in user/membranes.py. Every module that loops over membrane_types
# keeps working unchanged; what differs between entries is the lipid mixture, not the area.
import user.membranes as Membranes
membrane_types                      = Membranes.membrane_names()
filename_steep_system_MDP           = "steep_system.mdp"

# production
filename_steep_softcore_MDP         = "steep_softcore_system.mdp"
filename_steep_restr_MDP            = "steep_restr_system.mdp"
filename_NPT_equilibrate_system_MDP = "equilibrate_low_timestep.mdp"
filename_NPT_system_MDP             = "npt.mdp"

# tags
tags_steep_peptide_TOP              = ("$PEPTIDE_ITP",)
tags_steep_peptide_2pept_TOP        = ("$PROTEIN_NUM",)
tags_NVT_peptide_MDP                = ("$INCLUDE_CWD",)
tags_steep_system_MDP               = ("$INCLUDE_CWD",)
tags_NPT_equilibrate_system_MDP     = ("$INCLUDE_CWD",)
tags_NPT_system_MDP                 = ("$INCLUDE_CWD",)
tags_steep_softcore_MDP             = ("$INCLUDE_CWD",)
tags_steep_restr_MDP                = ("$INCLUDE_CWD",)
tags_SIM_dt                         = ("$SIM_DT",)
tags_SIM_nsteps                     = ("$SIM_NSTEPS",)


# Miscellaneous: generate_peptide.py #
STEEP_constraint_forceConstant      = "1000000" #"5000"
# append_ITP_noConstraints            = ".noConstraints.itp"
# append_ITP_seq2itp                  = ".seq2itp.itp"
# append_ITP_noConstr_Restr           = ".noConstr_Restr.itp"
# append_TOP_noConstraints            = ".noConstraints.top"
# append_TOP_noConstr_Restr           = ".noConstr_Restr.top"

# Miscellaneous: insert_peptide.py #
resnm_to_indexgroup_dict            = {
                                        "POPC"  : ("POPC", "Lipids"),
                                        "POPE"  : ("POPE", "Lipids"),
                                        "POPS"  : ("POPS", "Lipids"),
                                        "POPG"  : ("POPG", "Lipids"),
                                        "CDL2"  : ("CDL2", "Lipids"),
                                        "CDL0"  : ("CDL0", "Lipids"),
                                        "DPSM"  : ("DPSM", "Lipids"),
                                        "CHOL"  : ("CHOL", "Lipids"),
                                        "W"     : "Solvent", "NA"  : "Solvent", "CL"  : "Solvent", "ION"  : "Solvent",
                                        "GLY"   : "Protein", "ALA"  : "Protein", "ASP"  : "Protein", "ASN"  : "Protein", "GLU"  : "Protein",
                                        "GLN"   : "Protein", "VAL"  : "Protein", "LEU"  : "Protein", "ILE"  : "Protein", "MET"  : "Protein",
                                        "THR"   : "Protein", "SER"  : "Protein", "CYS"  : "Protein", "LYS"  : "Protein", "ARG"  : "Protein",
                                        "HIS"   : "Protein", "PHE"  : "Protein", "PRO"  : "Protein", "TRP"  : "Protein", "TYR"  : "Protein"      
                                      }

# Miscellaneous: production.py #
name_positive_ion_GENION            = "NA"
name_negative_ion_GENION            = "CL"
name_solvent_GENION                 = "W"
neutralize                          = True
# Keep every candidate's trajectory on disk. True for the class-validation run, where all
# trajectories are retained for inspection and reproducibility (9 peptides x 3 replicas is a
# few GB). For a STATE=0 GA campaign of thousands of candidates, run.py refuses True to avoid
# filling the disk; decide the campaign retention policy (all, top-N, or downsampled) with the
# disk numbers in hand before switching STATE and keeping this True.
WRITE_XTC                           = True

# Number of independent replicas per sequence in the fixed-sequence (STATE=1) validation run.
# run.py expands each sequence into this many replicas (dirs <SEQ>_rep1.._repN). Replicas differ
# by initial peptide orientation (a roll about the helix axis, set in the insertion module) and
# by the random velocity seed, so they are independent starting points. Single-replica 200 ns
# was dominated by per-trajectory noise; three replicas give a mean and a between-replica error.
# This governs the validation only; a STATE=0 GA campaign decides its own replica count.
N_REPLICAS                          = 3

# Keep the entire working directory (every intermediate trajectory: minimisation, relax,
# equilibration, production, plus structures and logs) next to each candidate's output on
# success, not only on failure. True for the validation, where all data is kept on disk for
# inspection and reproducibility. Disk-heavy over a campaign; revisit with WRITE_XTC for STATE=0.
RETAIN_ALL_TRAJECTORIES             = True

# Miscellaneous
FILE_WRITE_DEBUGGING                = False

# # Miscellaneous: seq2itp_aminoacids.py #
# generate_bond_length                = 0.350
# generate_sidechain_angle            = 90
# generate_box_factor                 = 50

# # Dependency: Seq2ITP #
# path_Seq2ITP                        = System.get_file_dir(__file__) + "/seq2itp/seq2itp.pl"
# #path_Seq2ITP_2pept                  = System.get_file_dir(__file__) + "/seq2itp/seq2itp_2pept.pl"

# Dependency: Porter5 #
path_Porter5                        = System.get_file_dir(__file__) + "/Porter5/Porter5.py"

# Dependency: modeller #
if pdb_gen == "Modeller":
    path_modeller                        = System.get_file_dir(__file__) + "/modeller/run_modeller.py"

# fitness
SIM_dt                    = 0.03 # ps
SIM_time                  = 200000 # ps  (200 ns per replica, class-validation)

# Fraction of the trajectory discarded before averaging. The original value (0.95) keeps
# only the last 5 percent, i.e. 20 ns out of 400 ns. That is a short window for a quantity as
# noisy as surface tension, and it is the dominant source of per-candidate error; the rerun
# averaging in the GA is what compensates for it. Exposed here so the trade-off between
# window length and equilibration can be examined rather than inherited.
ENERGY_init_fraction      = 0.95
ENERGY_init_time          = ENERGY_init_fraction*SIM_time

# Orientation correction factor (Eq. 3 of van Hilten et al. 2023). 'a' is the lower bound of
# the cosine window; the per-membrane 'b' now lives on each MembraneDef, since insertion
# depth depends on the bilayer, not on a position in a list.
a = 0.0

# other stuff
selection_sorting                   = True # goes in SelectionMethods, should be True if fitness is maximized, False if minimized (i.e. in case of free energy difference as fitness)             
SSP                                 = False
normalize_by_nr_of_beads            = False # if True, fitness values are devided by the number of beads in the peptide

rep_lengths                         = [24] # possible block lengths
repeats                             = [1] # possible repeat numbers
min_length                          = 24
max_length                          = 24 # maximal sequence length

# Net charge window, inclusive. The original value [-24, 24] is a symmetric upper bound that
# never binds for L=24, which is why the broken lower-bound branch of charge_correction was
# never exercised. A cationic window is the point of the exercise here, so the lower bound is
# now live and charge_correction has been rewritten to be able to raise charge, not only
# lower it. Keep R enabled in 'propensities' if the lower bound is positive, otherwise the
# correction has only K to work with.
charge_range                        = [3, 9]
# Dead. mutate_repeat is never called: its only call site in generate_new_population is
# commented out upstream. Independently, with rep_lengths = [24], repeats = [1] and
# min_length == max_length, any repeat mutation would produce a length the constraint rejects.
# Kept so the upstream signature still matches, but changing it has no effect.
repeat_mutations                    = [-1,1]
pept_position                       = [1/3] # position on the membrane (*box_length from the center)
rotate_vector                       = (0,0,90)
equilibration_time                  = 500 # ps
blocktime                           = 500 # ps
bootstrapping_sampling              = 10000
# Interaction-energy groups read from the production .edr. These were dead code in the
# original (energygrps is commented out in npt.mdp and compute_fitness only ever read
# SurfTen); enabling them costs one uncommented line in the mdp and no extra simulation,
# which is what lets the binding and damage objectives share the same trajectories.
ENERGY_group_terms                  = ["Coul-SR:Protein-Lipids",
                                       "LJ-SR:Protein-Lipids",
                                       "Coul-SR:Protein-Solvent",
                                       "LJ-SR:Protein-Solvent"]
ENERGY_select_groups_Protein1       = ["Coul-SR:Protein1-Lipids", "LJ-SR:Protein1-Lipids", "LJ-SR:Protein1-Solvent"]
ENERGY_select_groups_Protein2       = ["Coul-SR:Protein2-Lipids", "LJ-SR:Protein2-Lipids", "LJ-SR:Protein2-Solvent"]

# Per-membrane reference areas and free-membrane surface tensions moved to the MembraneDef
# entries in user/membranes.py: they are properties of a composition, not globals.
conversion_factor                   = 6.02214076e-2 # bar*nm to kJ*mol-1*nm-2
plotting                            = False
# Residue sampling weights, used for both initialisation and point mutation. Weights are
# multiplied by 1000 and truncated to int when the sampling list is built, so three decimals
# is the resolution and a weight of 0.0 removes the residue from the alphabet entirely.
#
# Two notes on what this dict actually was in the original code, because the comment there is
# misleading. First, the alphabet restriction lives here and not in GeneticAlgorithm.genes,
# which still holds all 20. Second, the original comment claims these are "propensities from
# ALPS-like peptides (Drin et al 2007)", but the active dict was uniform 0.1 over the ten
# helix-prone residues; the real ALPS propensities sat commented out below it and were never
# used. The behaviour matched the paper (uniform over a reduced alphabet); only the comment
# was wrong.
#
# R is enabled here, unlike the original. Excluding arginine is defensible for curvature
# sensors and indefensible for AMPs: it carries much of the cationic character of natural
# antimicrobial peptides and interacts with anionic headgroups differently from lysine.
#
# G and P remain excluded because the pipeline imposes a fixed helix; turning them on means
# the Martini secondary structure assignment is modelling a helix the sequence would not
# adopt. That premise is the main structural limitation of the whole approach and should be
# argued explicitly rather than relaxed here by accident.
#
# I stays excluded, but for a different reason than helix propensity. In the martini3001
# residue definitions shipped with vermouth, LEU and ILE carry identical bead types for both
# backbone (P2) and side chain (SC2), exactly as they do in martini22 (P5 / AC1). Bead type
# is what sets every nonbonded interaction, so the two residues are indistinguishable in
# hydrophobicity, lipid contact and partitioning: everything this pipeline measures. They
# differ only in side-chain geometry (LEU: bond, 0.363 nm, fc 7500; ILE: constraint,
# 0.341 nm), which is far below the noise of the surface-tension average. Searching over both
# spends evaluations separating two sequences the force field treats as one. Note that VAL is
# genuinely distinct in Martini 3 (SP2 / SC3), so this argument does not extend to it.
#
# This exclusion is enforced here and nowhere else: GeneticAlgorithm.genes still holds all 20
# residues, so setting propensities_mutation = False silently reintroduces I.
propensities_initial                = True # False = uniform over GeneticAlgorithm.genes (all 20)
propensities_mutation               = True # False = uniform over GeneticAlgorithm.genes (all 20)
propensities                        = {"A": 0.1,
                                       "C": 0.0,
                                       "D": 0.0,
                                       "E": 0.1,
                                       "F": 0.1,
                                       "G": 0.0,
                                       "H": 0.0,
                                       "I": 0.0,
                                       "K": 0.1,
                                       "L": 0.1,
                                       "M": 0.1,
                                       "N": 0.0,
                                       "P": 0.0,
                                       "Q": 0.1,
                                       "R": 0.1,
                                       "S": 0.1,
                                       "T": 0.0,
                                       "V": 0.0,
                                       "W": 0.1,
                                       "Y": 0.1}

# propensities_ALPS                   = {"A": 0.063,
#                                        "C": 0.000,
#                                        "D": 0.025,
#                                        "E": 0.024,
#                                        "F": 0.040,
#                                        "G": 0.132,
#                                        "H": 0.012,
#                                        "I": 0.053,
#                                        "K": 0.037,
#                                        "L": 0.111,
#                                        "M": 0.021,
#                                        "N": 0.036,
#                                        "P": 0.018,
#                                        "Q": 0.033,
#                                        "R": 0.026,
#                                        "S": 0.173,
#                                        "T": 0.101,
#                                        "V": 0.059,
#                                        "W": 0.014,
#                                        "Y": 0.022}

#################
## Objective   ##
#################
# Which scalarisation of the measured observables drives the GA. Both objectives read the
# same observables.json written by compute_observables, produced by the same pair of
# production runs, so switching between them costs nothing and running both on one dataset is
# a post-processing step rather than a second campaign.
#
#   "damage"  -> induced leaflet tension on the target membrane (approach B)
#   "binding" -> membrane partitioning proxy, target versus host (approach A)
OBJECTIVE                           = "damage"

# --- damage objective -------------------------------------------------------------------
# Delta gamma = gamma(free membrane) - gamma(membrane + peptide), at constant area. Positive
# means the peptide relieves tension in the leaflet it occupies, which is the same thing as
# saying it generates leaflet strain. That is the mechanical step upstream of toroidal pore
# formation, which is why it is used here as a damage proxy rather than a binding proxy.
#
# Caveat to carry into any claim made from these numbers: one peptide in a 128-lipid patch is
# roughly P:L 1:128, far below the ratio at which AMPs act cooperatively. Delta gamma here is
# the dilute-limit, per-peptide tension generation. Extrapolating it to lysis is an
# assumption, not a measurement.
DAMAGE_host_tolerance               = 0.0   # kJ mol-1 nm-2 of host strain tolerated unpenalised
DAMAGE_host_weight                  = 1.0   # hinge penalty weight above that tolerance

# --- binding objective ------------------------------------------------------------------
# Proxy choice is deliberately left open: it should be decided from the calibration run, not
# assumed. Short-range interaction energy is cheap and already in the .edr, but it is
# enthalpy of interaction, with no entropy and no desolvation cost, and it is unreliable when
# comparing peptides of different composition. Insertion depth is more robust and less
# informative. Thermodynamic integration is the correct quantity and does not fit in the loop;
# reserve it for re-ranking the top candidates.
#
#   "interaction_energy" -> Coul-SR + LJ-SR between peptide and lipids
#   "depth"              -> centre-of-mass insertion depth d_mem
#   "combined"           -> weighted sum of the two, z-scored within the generation
BINDING_proxy                       = "interaction_energy"
BINDING_depth_weight                = 0.5   # only used when BINDING_proxy == "combined"
BINDING_host_weight                 = 1.0   # weight of the host term in the selectivity difference

# Apply the cosine orientation factor (Eq. 3) to the fitness. Shared by both objectives:
# a peptide that has left the surface, gone transmembrane or dissolved is not measuring what
# either objective means to measure.
APPLY_orientation_factor            = True

# Below this orientation factor the candidate is rejected outright rather than scaled.
#
# Scaling alone is not enough, and the reason is worth stating because it is not obvious. The
# factor multiplies only a positive raw score, so a detached peptide with factor 0 scores
# exactly 0.0 -- and on a maximised fitness, 0.0 beats every ATTACHED peptide whose raw score
# is negative. In early generations most random sequences have a negative delta_gamma, so
# "leave the membrane" was a reachable local optimum worth more than most honest attempts, and
# the fitness trace would show only a plateau at zero.
#
# Rejecting instead of scaling is also the physically honest reading: the factor is not a
# penalty on a real measurement, it is a statement that the measurement does not mean what the
# objective assumes. Set to 0.0 to reject only a peptide that is never at the interface, or
# raise it to demand that a candidate spend most of the window there.
ORIENTATION_reject_below            = 0.05

# Score given to a rejected candidate. Below any attainable raw value -- real delta_gamma is of
# order 1 to 10 kJ mol-1 nm-2 -- so such a candidate is never selected as a parent, and well
# above run.WORST_FITNESS (-1e6) so that the report can tell a peptide that detached from an
# evaluation that failed to run at all. The two mean different things: the first is a result.
DETACHED_FITNESS                    = -1.0e3


###############
## Functions ##
###############

# # Seq2ITP #
# def command_seq2ITP(path_sequence, path_structure, path_ITP_out):
#     command = "perl "
#     command += path_Seq2ITP
#     command += " --sequence " + path_sequence
#     command += " --structure " + path_structure
#     command += " --itp " + path_ITP_out
#     command += " --noelastic"
#     return command

# Porter5 #
def command_Porter5(path_sequence):
    command = "python3 "
    command += path_Porter5
    command += " -i " + path_sequence
    command += " --cpu 1 "
    command += " --fast"
    return command

# modeller #
def command_modeller(path_sequence, path_structure, path_PDB, name, ali, outdir):
    command = "python "
    command += path_modeller + " "
    command += path_sequence + " "
    command += path_structure + " "
    command += path_PDB + " "
    command += name + " "
    command += ali + " "
    command += outdir
    return command


# Resolved from the environment rather than hardcoded. The original pointed at an absolute
# path inside one author's home directory and broke on first run anywhere else.
path_martinize2                     = os.environ.get("MARTINIZE2", "martinize2")


# Neutral peptide termini, via martinize2's -nt. Charged termini are the martinize2 default and
# were what upstream used. Two reasons to turn them off here. They add a charge at each end whose
# interaction with the membrane depends on how the helix happens to be oriented, which is
# variability in the measured surface tension that has nothing to do with the sequence being
# evaluated. And they are inconsistent with the experimental literature this pipeline is
# calibrated against, where at least one terminus is usually capped. Reference 2 made the same
# choice for the same reason.
#
# This does NOT change the net charge. Verified against martinize2 (vermouth 0.15.0): the default
# charges the N-terminus +1 and the C-terminus -1, so they cancel, and the same sequence comes out
# at the same total charge either way. Only the two terminal backbone bead types change (Q5 to P6
# in martini3001). charge_range is therefore the side-chain charge in both cases, and it does not
# need adjusting for this setting.
NEUTRAL_TERMINI                     = True


def command_martinize2(path_PDB_in, structure, path_PDB_out, path_TOP, ff=ff, elastic=False):
    command = path_martinize2
    command += " -f " + path_PDB_in
    command += " -x " + path_PDB_out
    command += " -o " + path_TOP
    command += " -ss " + structure
    command += " -ff " + ff
    if NEUTRAL_TERMINI:
        command += " -nt"
    if elastic:
        command += " -elastic"
    return command

def command_whereamI():
    command = "pwd"
    return command

# def command_seq2ITP_2pept(path_sequence, path_structure, path_ITP_out, molname):
#     command = "perl "
#     command += path_Seq2ITP_2pept
#     command += " --sequence " + path_sequence
#     command += " --structure " + path_structure
#     command += " --itp " + path_ITP_out
#     command += " --noelastic"
#     command += " --molname " + molname
#     return command

# The orientation factor lives in user/modules/compute_observables.orientation_factor, which is
# the single implementation the pipeline uses. Two variants used to sit here: one duplicating it
# including a window bound that has since been corrected there, and one implementing a different
# window shape altogether. Both were dead, and an alternative definition of the same physical
# weighting next to the settings is a trap for anyone reconstructing what was applied.
