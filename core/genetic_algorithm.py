import sys
import os
import pickle
import user.usersettings as UserSettings
import numpy as np
import math

from random         import Random
from datetime       import datetime
from core.system    import mprint

# Guarded, like the import in core/system.py. The module already tests "if MPI is not None"
# throughout, so the optional path was intended; only the import contradicted it.
MPI = None
try:
    from mpi4py import MPI
except ImportError:
    pass

class _SimulationBuffer:
    def __init__(self, ptr_scoring_function):
        self.__buffer               = set()
        self.ptr_scoring_function   = ptr_scoring_function

    def __sim_buffer(self, buffer_data):
        output_dict = {}

        for chromosome in buffer_data:
            chromosome.compute_fitness(output_dict, self.ptr_scoring_function)

        return output_dict

    def add(self, chromosome_set):
        self.__buffer = self.__buffer.union(chromosome_set)

    def simulate(self): # Returns output dict (seq string, fitness) describing all of the simulated chromosomes, without references to the object
        state_MPI = False
        if MPI is not None:
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            if size > 1:
                state_MPI = True

        output_dict = None
        
        if state_MPI:
            output_dict = self.__simulate_MPI(comm)
        else:
            output_dict = self.__sim_buffer(self.__buffer)

        self.__buffer = set()
        return output_dict

    # returns an output dict, similar to self.simulate() (i.e. individual __sim_buffer calls return dicts from different threads, combine and return combined dict)
    def __simulate_MPI(self, comm):
        chunk = lambda lst, N_chunks: [lst[i::N_chunks] for i in range(N_chunks)]

        size = comm.Get_size()
        rank = comm.Get_rank()

        number_of_ranks = size
        if number_of_ranks > len(self.__buffer):
            number_of_ranks = len(self.__buffer)

        if rank == 0:   # Master
            # 0. Send #
            list_of_buffers = chunk(list(self.__buffer), number_of_ranks)
            for i in range(1, number_of_ranks):
                comm.send(list_of_buffers[i], dest=i)

            # 2. Simulate own part #
            output_dict = self.__sim_buffer(list_of_buffers[0])

            # 4. Receive and combine #
            for i in range(1, number_of_ranks):
                slave_dict  = comm.recv(source=i)
                output_dict = {**output_dict, **slave_dict}

            # 5. Send to slaves #
            for i in range(1, number_of_ranks):
                comm.send(output_dict, dest=i)

            return output_dict

        else:           # Slave
            # 1. Receive #
            buffer_data     = comm.recv(source=0)

            # 2. Simulate own part #
            output_dict     = self.__sim_buffer(buffer_data)

            # 3. Send results to master #
            comm.send(output_dict, dest=0)

            # 6. Receive combined #
            output_dict     = comm.recv(source=0)

            return output_dict

        raise RuntimeError("Should not happen.")


#######################
## String operations ##
#######################

class _StringChromosome:
    def __init__(self, initial_sequence):
        """
        Args:
            initial_sequence (list(string)): List of strings describing the genes.
        """

        self.sequence = initial_sequence.copy()
        self.fitness  = None
        self.score    = None

    def __eq__(self, other):
        return self.get_sequence_string() == other.get_sequence_string()
    def __ne__(self, other):
        return self.get_sequence_string() != other.get_sequence_string()
    def __hash__(self):
        return hash(self.get_sequence_string())

    def copy(self):
        return _StringChromosome(self.sequence)

    def mutate(self, P_gene, genes, RNG):
        for i in range(len(self.sequence)):
            if RNG.random() < P_gene:
                self.sequence[i] = RNG.choice(genes)

    def mutate_propensities(self, P_gene, RNG, propensities): # point mutation at random position (following user defined propensities)
        for i in range(len(self.sequence)):
            if RNG.random() < P_gene:
                new_AA = RNG.choice([AA for AA in propensities for prop in range(int(propensities[AA]*1000))])
                self.sequence[i] = new_AA

    def deletion(self, repeat, RNG): # deletion at random position
        seq_length = float(len(self.sequence) * repeat)
        #P_deletion = (seq_length - min(UserSettings.rep_lengths)) / (UserSettings.max_length - min(UserSettings.rep_lengths))
        P_deletion = 0
        if RNG.random() < P_deletion and len(self.sequence) > min(UserSettings.rep_lengths):
            i = RNG.choice([*range(len(self.sequence))])
            #print("del", self.sequence,"-->", self.sequence[:i]+self.sequence[(i+1):])
            self.sequence = self.sequence[:i] + self.sequence[(i+1):]

    def insertion(self, repeat, RNG, propensities): # insertion at random position
        seq_length = float(len(self.sequence) * repeat)
        #P_insertion = 1 - ((seq_length - min(UserSettings.rep_lengths)) / (UserSettings.max_length - min(UserSettings.rep_lengths)))
        P_insertion = 0
        if RNG.random() < P_insertion and len(self.sequence) < max(UserSettings.rep_lengths) and (len(self.sequence) + 1) * repeat <= UserSettings.max_length:
            i = RNG.choice([*range(len(self.sequence))])
            new_AA = RNG.choice([AA for AA in propensities for prop in range(int(propensities[AA]*1000))])
            #print("ins", self.sequence,"-->", self.sequence[:i] + [new_AA] + self.sequence[i:])
            self.sequence = self.sequence[:i] + [new_AA] + self.sequence[i:]

    CATIONIC_RESIDUES = ("K", "R")
    ANIONIC_RESIDUES = ("D", "E")

    def net_charge(self):
        """Net charge at pH 7. Histidine is treated as neutral, as in the original."""
        pos = sum(self.sequence.count(aa) for aa in _StringChromosome.CATIONIC_RESIDUES)
        neg = sum(self.sequence.count(aa) for aa in _StringChromosome.ANIONIC_RESIDUES)
        return pos - neg

    def charge_correction(self, repeat, RNG, propensities, charge_range, max_iter=200):
        """
        Mutate residues until the net charge falls inside charge_range.

        Rewritten. The original could only move the net charge toward zero: it picked a
        charged residue and replaced it with a neutral one, whichever direction was needed.
        That satisfies a symmetric upper bound such as the shipped [-24, 24], which for L=24
        never binds anyway, and cannot satisfy a lower bound at all. It also raised NameError
        on an exactly neutral sequence outside the window, because neither of its two
        unguarded ifs assigned charge_list, and it had no iteration guard.

        A cationic window is the whole point for antimicrobial peptides, so this version can
        raise charge as well as lower it, and gives up cleanly rather than looping when the
        alphabet makes the window unreachable (for instance a positive lower bound with every
        cationic residue zeroed out in the propensities).
        """
        lo, hi = min(charge_range), max(charge_range)

        def sampling_pool(exclude=(), restrict_to=None):
            # Mirrors the weighted sampling used elsewhere: weight * 1000, truncated to int.
            pool = []
            for aa, weight in propensities.items():
                if aa in exclude:
                    continue
                if restrict_to is not None and aa not in restrict_to:
                    continue
                pool.extend([aa] * int(weight * 1000))
            return pool

        charged = set(_StringChromosome.CATIONIC_RESIDUES) | set(_StringChromosome.ANIONIC_RESIDUES)

        for _ in range(max_iter):
            charge = self.net_charge() * repeat
            if lo <= charge <= hi:
                return True

            if charge < lo:
                # Need more positive charge: either neutralise an acidic residue or turn a
                # non-charged position cationic. Both raise the net charge by one.
                donors = [n for n, aa in enumerate(self.sequence)
                          if aa in _StringChromosome.ANIONIC_RESIDUES]
                if donors:
                    i = RNG.choice(donors)
                    pool = sampling_pool(exclude=charged)
                else:
                    i_candidates = [n for n, aa in enumerate(self.sequence) if aa not in charged]
                    if not i_candidates:
                        return False
                    i = RNG.choice(i_candidates)
                    pool = sampling_pool(restrict_to=_StringChromosome.CATIONIC_RESIDUES)
            else:
                # Need less positive charge: mirror image of the branch above.
                donors = [n for n, aa in enumerate(self.sequence)
                          if aa in _StringChromosome.CATIONIC_RESIDUES]
                if donors:
                    i = RNG.choice(donors)
                    pool = sampling_pool(exclude=charged)
                else:
                    i_candidates = [n for n, aa in enumerate(self.sequence) if aa not in charged]
                    if not i_candidates:
                        return False
                    i = RNG.choice(i_candidates)
                    pool = sampling_pool(restrict_to=_StringChromosome.ANIONIC_RESIDUES)

            if not pool:
                # The alphabet cannot express the required move; the caller's constraint
                # function will reject the sequence rather than this looping forever.
                return False

            self.sequence[i] = RNG.choice(pool)

        return lo <= self.net_charge() * repeat <= hi

    def mutate_repeat(self, repeat, mutation_list, RNG): # change the repeat number
        #P_mutate_repeat = 1/len(self.sequence)
        P_mutate_repeat = 0
        if len(self.sequence) * repeat > UserSettings.max_length:
            new_repeat = math.floor(float(UserSettings.max_length)/len(self.sequence))
            return new_repeat
        elif RNG.random() < P_mutate_repeat:
            new_repeat = repeat + RNG.choice(mutation_list)
            if min(UserSettings.repeats) <= new_repeat * len(self.sequence) <= UserSettings.max_length:
                #print("repeat mutation","x"+str(repeat),"-->","x"+str(new_repeat),"\tOLD_LENGTH=",repeat*len(self.sequence),"\tNEW_LENGTH=",new_repeat*len(self.sequence))
                return new_repeat
            else:
                return repeat
        else:
            return repeat

    def scramble(self, RNG): # randomly scramble the sequence
        P_scramble = 0
        if RNG.random() < P_scramble:
            RNG.shuffle(self.sequence)

    def crossover(self, other, P_crossover, RNG):
        if len(self.sequence) > 3:
            if RNG.random() < P_crossover:
                position = RNG.randint(1, len(self.sequence) - 2)

                temp                        = self.sequence[position:]
                self.sequence[position:]    = other.sequence[position:]
                other.sequence[position:]   = temp

    def crossover_replen(self, other, P_crossover, RNG): # crossover for repeated block peptides
        if RNG.random() < P_crossover:
            A = self.sequence
            B = other.sequence

            positionA = RNG.randint(1, len(A) - 1)
            positionB = RNG.randint(1, len(B) - 1)
            A1 = A[:positionA]
            A2 = A[positionA:]
            B1 = B[:positionB]
            B2 = B[positionB:]

            #print("".join(A1),"X","".join(A2))
            #print("".join(B1),"X","".join(B2))

            # COMPOSE CHILD 1
            if len(A) <= len(A1):
                child1 = A1[-abs((len(A)-len(B2))):] + B2
            else:
                child1 = A1 + B2[:(len(A)-len(A1))]

            # COMPOSE CHILD 2
            if len(B) <= len(B1):
                child2 = B1[-abs((len(B)-len(A2))):] + A2
            else:
                child2 = B1 + A2[:(len(B)-len(B1))]

            # TRUNCATE CHILDREN
            if len(A) < len(child1):
                child1 = child1[:len(A)]
            if len(B) < len(child2):
                child2 = child2[:len(B)]

            # COMPLETE CHILDREN
            if len(A) > len(child1):
                child1 = child1 + (10*B)[:(len(A) - len(child1))]
            if len(B) > len(child2):
                child2 = child2 + (10*A)[:(len(B) - len(child2))]

            self.sequence = child1
            other.sequence = child2

    def compute_fitness(self, output_dict, ptr_scoring_function):
        if self.fitness is None:
            self.fitness = ptr_scoring_function(self.sequence)
        output_dict[self.get_sequence_string()] = self.fitness

    def get_sequence_string(self):
        return "".join(self.sequence)

    def toString(self):
        return self.get_sequence_string() + "\t" + str(self.fitness)

    def find_replen(self): # find the minimal repeated block length
        rep_len_list = UserSettings.rep_lengths
        # Initialised to the no-repeat case. Both were previously assigned only inside the loop
        # body, so a sequence whose length is not divisible by any entry in rep_lengths reached
        # the assignment below with nothing bound and died mid-campaign as "UnboundLocalError:
        # cannot access local variable 'new_seq'". That is reachable by editing rep_lengths, and
        # by restarting from a checkpoint written under different settings.
        repeat = 1
        new_seq = self.sequence
        for l in rep_len_list:
            if len(self.sequence) % l != 0:
                continue
            if self.sequence[:l] == self.sequence[-l:]:
                nlist = [*range(1,int(float(len(self.sequence))/l))]
                if len(nlist) > 0:
                    ntrue = 0
                    for n in nlist:
                        if self.sequence[:l] == self.sequence[(n*l):((n+1)*l)]:
                            ntrue+=1
                    if len(nlist) == ntrue:
                        repeat = int(float(len(self.sequence))/l)
                        new_seq = self.sequence[:l]
                        break
                else:
                    repeat = 1
                    new_seq = self.sequence
        self.sequence = new_seq
        return repeat





###########################
## Population operations ##
###########################

class _Population:
    class _iterator:
        def __init__(self, population):
            self.population_list    = list(population)
            self.index              = 0

        def __next__(self):
            if self.index >= len(self.population_list):
                raise StopIteration
            self.index += 1
            return self.population_list[self.index - 1]

    def __new_population_from_output_dict(self, output_dict):
        """Replace the current data with data from output_dict dictionary.
        NOTE: updated to include accounting for reruns. 
        If a sequence has been simulated more than once, the resulting fitness will be the average of all the runs of that sequence.
        """
        self.population = set()

        for sequence in output_dict:
            if sequence not in self.n_run_dict:
                self.n_run_dict[sequence] = 1
            else:
                self.n_run_dict[sequence] += 1
                prev_fitness = self.chromosome_pool[sequence][0]
                prev_ddF = self.chromosome_pool[sequence][1]
                output_dict[sequence][0] = (prev_fitness * (self.n_run_dict[sequence] - 1) + output_dict[sequence][0]) / self.n_run_dict[sequence] # Weighted averages
                output_dict[sequence][1] = (prev_ddF * (self.n_run_dict[sequence] - 1) + output_dict[sequence][1]) / self.n_run_dict[sequence] # Weighted averages

            chromosome          = _StringChromosome(list(sequence))
            chromosome.fitness  = output_dict[sequence][0]

            self.population.add(chromosome)

    def __generate_sequence(self, length):
        return self.random.choices(self.genes, k=length)

    def __generate_sequence_propensities(self, length, propensities):
        return self.random.choices([AA for AA in propensities for prop in range(int(propensities[AA]*1000))], k=length)


    def __init__(self, ptr_genetic_algorithm):
        self.ptr_genetic_algorithm  = ptr_genetic_algorithm
        self.random                 = self.ptr_genetic_algorithm.RNG
        self.genes                  = self.ptr_genetic_algorithm.genes
        self.population             = set()
        self.chromosome_pool        = {} # Dict containing all simulated sequences and corresponding fitnesses
        self.n_run_dict             = {} # Dict containing all simulated sequences and how often they were simulated
        self.simulation_buffer      = _SimulationBuffer(self.ptr_genetic_algorithm.ptr_scoring_function)

        population_size = ptr_genetic_algorithm.population_size

        # Both loops below used to be unbounded. Neither can make progress when the settings are
        # mutually incompatible, and the symptom is identical to a slow start: the last thing
        # printed is "Checking settings and initializing the genetic algorithm.." and then a
        # pegged core, indefinitely, with no further output. Two realistic ways in: a charge
        # window tightened past what the enabled residues can reach, so the constraint rejects
        # essentially every draw; and a rep_lengths entry larger than max_length, so no repeat
        # count satisfies the inner condition.
        max_attempts = 1000 * max(population_size, 1)
        attempts = 0
        rejected_by = {}

        while len(self.population) < population_size:
            attempts += 1
            if attempts > max_attempts:
                import core.errors as Errors
                Errors.fail(
                    "could not build the initial population: %d of %d sequences after %d draws"
                    % (len(self.population), population_size, attempts),
                    where="core/genetic_algorithm.py, population initialisation",
                    why=("every draw was rejected by the constraints, so the search space the "
                         "settings describe is empty or nearly so.\n"
                         "charge_range = %s, and the residues with a non-zero propensity are %s\n"
                         "rep_lengths = %s, repeats = %s, min_length = %s, max_length = %s"
                         % (UserSettings.charge_range,
                            "".join(sorted(a for a, p in UserSettings.propensities.items()
                                           if p > 0)) or "(none)",
                            UserSettings.rep_lengths, UserSettings.repeats,
                            UserSettings.min_length, UserSettings.max_length)),
                    fix=("widen charge_range, or enable more cationic residues (K, R) in "
                         "propensities, or relax the length settings, in user/usersettings.py"))

            repeat = 1000
            replen = self.random.choice(UserSettings.rep_lengths)
            inner = 0
            while replen * repeat > UserSettings.max_length:
                inner += 1
                if inner > 10000:
                    import core.errors as Errors
                    Errors.fail(
                        "no repeat count can fit a block of %d residues into max_length %d"
                        % (replen, UserSettings.max_length),
                        where="core/genetic_algorithm.py, population initialisation",
                        why=("sequences are built as a block of rep_lengths repeated a number "
                             "of times drawn from repeats, and every value in repeats = %s "
                             "overshoots max_length = %d for this block length"
                             % (UserSettings.repeats, UserSettings.max_length)),
                        fix=("in user/usersettings.py make max_length at least the largest "
                             "rep_lengths entry (%s), or add a smaller value to repeats"
                             % max(UserSettings.rep_lengths)))
                repeat = self.random.choice(UserSettings.repeats)

            if UserSettings.propensities_initial:
                initial_sequence = self.__generate_sequence_propensities(replen, UserSettings.propensities)
            else:
                initial_sequence = self.__generate_sequence(replen)
            if self.ptr_genetic_algorithm.ptr_constraint_function(initial_sequence*repeat):
                rejected_by["constraint"] = rejected_by.get("constraint", 0) + 1
                continue

            self.population.add(_StringChromosome(initial_sequence*repeat))

        if attempts > 2 * population_size:
            # Not an error, but the acceptance rate is the one number that says how close the
            # settings are to the empty-search-space failure above.
            mprint("Initial population: %d sequences accepted from %d draws (%.0f%%). A low rate "
                   "means the constraints and the propensities are close to incompatible."
                   % (population_size, attempts, 100.0 * population_size / attempts))

    def __iter__(self):
        return _Population._iterator(self.population)

    def __write_checkpoint_file(self):
        if MPI is not None: # Ensure only master writes
            if not MPI.COMM_WORLD.Get_rank() == 0:
                return

        # The iteration is part of the restart state. Without it a resumed run restarted its
        # numbering at 0, so the progress file gained a second "Iteration 0" block and
        # tools/ga_report.py, which keys generations by that number, merged the two into one.
        pickle_data = (self.chromosome_pool, self.n_run_dict, self.population,
                       self.ptr_genetic_algorithm.current_iteration)

        filepath_checkpoint = self.ptr_genetic_algorithm.filepath_checkpoint
        if filepath_checkpoint is not None:
            # Written to a temporary file and moved into place, because this overwrites the only
            # record of a campaign that may be days old. A kill or a full disk during the write
            # used to truncate it, and the next --restart then failed on a half-written pickle
            # with nothing to fall back to. os.replace is atomic within a filesystem, so the
            # checkpoint on disk is always either the previous generation or this one.
            temporary = filepath_checkpoint + ".tmp"
            with open(temporary, 'wb') as file:
                pickle.dump(pickle_data, file, pickle.DEFAULT_PROTOCOL)
                file.flush()
                os.fsync(file.fileno())
            os.replace(temporary, filepath_checkpoint)
    
    def read_checkpoint_file(self, filepath_checkpoint):
        import core.errors as Errors

        Errors.require_file(
            filepath_checkpoint,
            why=("--restart resumes a campaign from the checkpoint the genetic algorithm "
                 "writes after every generation"),
            fix=("point --restart at the checkpoint in the output directory of the run you "
                 "want to continue, or drop --restart to start a new campaign"))

        try:
            with open(filepath_checkpoint, 'rb') as file:
                loaded = pickle.load(file)
            # Checkpoints written before the iteration was recorded hold three elements. They
            # still load; such a run resumes with its numbering restarted, as it did before.
            if len(loaded) == 4:
                (chromosome_pool, n_run_dict, population, iteration) = loaded
            else:
                (chromosome_pool, n_run_dict, population) = loaded
                iteration = 0
        except Exception as exc:
            # Anything unpickling raises here means the file is not a usable checkpoint. The
            # distinction worth drawing for the operator is between a truncated file, which is
            # what a kill during the write produces, and the wrong file entirely.
            Errors.fail("the checkpoint could not be read (%s: %s)" % (type(exc).__name__, exc),
                        where=os.path.abspath(filepath_checkpoint),
                        why=("it is either truncated, which is what a kill or a full disk "
                             "during the write leaves behind, or it is not a checkpoint"),
                        fix=("check the file size against the other campaigns' checkpoints; if "
                             "it is truncated the generation it was written for is lost and the "
                             "campaign has to restart from the last good checkpoint"))

        self.chromosome_pool = chromosome_pool
        self.n_run_dict = n_run_dict
        self.population = population
        self.ptr_genetic_algorithm.restart_from_iteration = iteration
        mprint("Restarted from %s: %d sequences already evaluated, population of %d."
               % (os.path.abspath(filepath_checkpoint), len(chromosome_pool), len(population)))

    def evaluate_fitness(self):
        # Simulate #
        self.simulation_buffer.add(self.population)
        output_dict = self.simulation_buffer.simulate()

        # Post-process #
        self.__new_population_from_output_dict(output_dict)
        self.chromosome_pool = {**self.chromosome_pool, **output_dict}

        self.__write_checkpoint_file()

        return output_dict

    def _compute_score(self):
        for individual in self.population:
            individual.score = individual.fitness + self.ptr_genetic_algorithm.ptr_diversity_function(
                                                        individual.sequence, 
                                                        self.population, 
                                                        self.ptr_genetic_algorithm.current_iteration/self.ptr_genetic_algorithm.iterations
                                                    )

    def generate_new_population(self, N_parents):
        self._compute_score()

        selected_population = self.ptr_genetic_algorithm.ptr_selection_function(self.population, N_parents, self.ptr_genetic_algorithm)

        new_population      = set()
        N_elites            = self.ptr_genetic_algorithm.N_elites
        N_rerun_elites      = self.ptr_genetic_algorithm.N_rerun_elites

        if not N_rerun_elites == 0:
            elites = SelectionMethods._n_rerun_elites(self.population, self.n_run_dict, self.ptr_genetic_algorithm)
            for elite in elites:
                new_population.add(elite.copy())

        if not N_elites == 0:
            elites = SelectionMethods._fitness_elites(self.population, N_elites, self.ptr_genetic_algorithm)
            for i in range(N_elites):
                new_population.add(elites[i].copy())

        P_crossover     = self.ptr_genetic_algorithm.P_crossover
        P_mutate        = self.ptr_genetic_algorithm.ptr_P_mutate_function(float(UserSettings.max_length), self.ptr_genetic_algorithm.current_iteration/self.ptr_genetic_algorithm.iterations) # second argument does not do anything (is not used in function)
        random          = self.random

        # The same unbounded spin the initialisation loop was hardened against, reachable here
        # on generation 50 of a 100-generation campaign rather than at startup, where nobody is
        # watching. Two ways in: settings edited between a run and its --restart so that the
        # constraints now reject every child, and a population converged enough that the
        # duplicate and reverse-duplicate filtering below can no longer fill the quota.
        target = len(self.population)
        max_attempts = 1000 * max(target, 1)
        attempts = 0

        while len(new_population) < target:
            attempts += 1
            if attempts > max_attempts:
                import core.errors as Errors
                Errors.fail(
                    "could not build generation %d: %d of %d sequences after %d attempts"
                    % (self.ptr_genetic_algorithm.current_iteration, len(new_population),
                       target, attempts),
                    where="core/genetic_algorithm.py, generate_new_population",
                    why=("every crossover and mutation was rejected, either by the constraints "
                         "or as a duplicate of a sequence already in the generation.\n"
                         "charge_range = %s, residues with a non-zero propensity: %s\n"
                         "min_length = %s, max_length = %s, parents = %d"
                         % (UserSettings.charge_range,
                            "".join(sorted(a for a, p in UserSettings.propensities.items()
                                           if p > 0)) or "(none)",
                            UserSettings.min_length, UserSettings.max_length,
                            len(selected_population))),
                    fix=("if the settings were changed since this campaign started, restore "
                         "them; a constraint the current population cannot satisfy makes every "
                         "child illegal. Otherwise the population has converged: raise the "
                         "mutation rate or the parent count."))
            A = random.choice(selected_population).copy()
            B = random.choice(selected_population).copy()
            #print("PARENTS\t"+"".join(A.sequence)+"\t"+"".join(B.sequence))
            # translate to block
            repeatA = A.find_replen()
            repeatB = B.find_replen()
            #print("\nBEFORE\t","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # crossover
            if len(UserSettings.rep_lengths) > 1:
                if len(A.sequence) == len(B.sequence):
                    A.crossover(B, P_crossover, self.random)
                else:
                    A.crossover_replen(B, P_crossover, self.random)
            if len(UserSettings.rep_lengths) == 1:
                A.crossover(B, P_crossover, self.random)
            #print("AFTER CROSSOVER","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # point mutations
            if UserSettings.propensities_mutation == True:
                A.mutate_propensities(P_mutate, self.random, UserSettings.propensities)
                B.mutate_propensities(P_mutate, self.random, UserSettings.propensities)
            if UserSettings.propensities_mutation == False:
                A.mutate(P_mutate, self.genes, self.random)
                B.mutate(P_mutate, self.genes, self.random)
            #print("AFTER MUTATION","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # repeat mutation
            #if len(UserSettings.rep_lengths) > 1 and len(UserSettings.repeat_mutations) > 1:
            #    repeatA = A.mutate_repeat(repeatA, UserSettings.repeat_mutations, self.random)
            #    repeatB = B.mutate_repeat(repeatB, UserSettings.repeat_mutations, self.random)

            # deletion
            #if len(UserSettings.rep_lengths) > 1:
            #    A.deletion(repeatA, self.random)
            #    B.deletion(repeatB, self.random)
            #print("AFTER DELETION","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # insertion
            #if len(UserSettings.rep_lengths) > 1:
            #    A.insertion(repeatA, self.random, UserSettings.propensities)
            #    B.insertion(repeatB, self.random, UserSettings.propensities)
            #print("AFTER INSERTION","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # charge correction
            A.charge_correction(repeatA, self.random, UserSettings.propensities, UserSettings.charge_range)
            B.charge_correction(repeatB, self.random, UserSettings.propensities, UserSettings.charge_range)
            #print("AFTER CHARGE CORR.","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # scramble
            #A.scramble(self.random)
            #B.scramble(self.random)
            #print("AFTER SCRAMBLE","".join(A.sequence),"x"+str(repeatA),"\t","".join(B.sequence),"x"+str(repeatB))

            # translate back to full sequence
            #print("RESULT\t", "".join(A.sequence),"x"+str(repeatA), "\t", "".join(B.sequence),"x"+str(repeatB))
            A.sequence = A.sequence * repeatA
            B.sequence = B.sequence * repeatB
            #print("CHILDREN\t"+"".join(A.sequence)+"\t"+"".join(B.sequence)+"\n")

            if (self.ptr_genetic_algorithm.ptr_constraint_function(A.sequence)
                or self.ptr_genetic_algorithm.ptr_constraint_function(B.sequence)):
                continue


            # Fill new population, and check and remove duplicates and reverse duplicates
            ind_list = []
            for ind in new_population:
                ind_list.append("".join(ind.sequence))

            if "".join(A.sequence) in ind_list or "".join(A.sequence[::-1]) in ind_list:
                pass
            else:
                new_population.add(A)
            if len(new_population) == len(self.population):
                break

            if "".join(B.sequence) in ind_list or "".join(B.sequence[::-1]) in ind_list:
                pass
            else:
                new_population.add(B)
            if len(new_population) == len(self.population):
                break

        population_size = self.ptr_genetic_algorithm.population_size
        while len(new_population) > population_size: # Old code to avoid population overflow. Likely unnecessary but kept here to ensure no overflow.
            new_population.pop() #NOTE: Randomly removes individual from population when called.
        self.population = new_population

    def select_final(self, N):
        """
        Retrieve N highest fitness chromosomes from all iterations.

        NOTE: Not always the best result.
        Usually you want a high fitness result that has been rerun several times,
        resulting in better sampling.
        """
        sorted_population = sorted(self.chromosome_pool.items(), key=lambda x: x[1][0], reverse=UserSettings.selection_sorting)
        return sorted_population[:N]




#######################
## Selection methods ##
#######################

class SelectionMethods:
    #format: func(population (list), N_parents (int), ptr_genetic_algorithm(, others))

    @staticmethod
    def _n_rerun_elites(population, N_rerun_dict, ptr_genetic_algorithm):
        result = []
        sorted_population = sorted(population, key=lambda x: x.fitness, reverse=UserSettings.selection_sorting)

        count = 0
        for individual in sorted_population:
            if count >= ptr_genetic_algorithm.N_rerun_elites:
                break
            if N_rerun_dict["".join(individual.sequence)] >= ptr_genetic_algorithm.N_rerun_threshold:
                result.append(individual)
                count += 1

        return result

    @staticmethod
    def _fitness_elites(population, N_elites, ptr_genetic_algorithm):
        sorted_population = sorted(population, key=lambda x: x.fitness, reverse=UserSettings.selection_sorting)
        return sorted_population[:N_elites]

    @staticmethod
    def best_N(population, N_parents, ptr_genetic_algorithm): 
        """Retrieve N highest fitness chromosomes from current iteration."""

        sorted_population = sorted(population, key=lambda x: x.score, reverse=UserSettings.selection_sorting)
        return sorted_population[:N_parents]

    @staticmethod
    def tournament(population, N_parents, ptr_genetic_algorithm): # No replacement
        pop         = list(population)
        selected    = []

        for _ in range(N_parents):
            tour    = sorted(ptr_genetic_algorithm.RNG.sample(pop, ptr_genetic_algorithm.selfun_tournament_size), key=lambda x: x.score, reverse=UserSettings.selection_sorting)
            winner  = tour[0]

            pop.remove(winner)
            selected.append(winner) #TODO use probability instead to select winner? (1: p; 2:p*(1-p); 3:p*(1-p)^2; etc...)?

        return selected

    @staticmethod
    def exponential_rank(population, N_parents, ptr_genetic_algorithm): # No replacement
        c           = ptr_genetic_algorithm.selfun_exponential_rank_C
        p           = lambda i, N: (c - 1)/(c**N - 1) * (c**(N - i))
        selected    = []
        pop         = sorted(list(population), key=lambda x: x.score)

        while len(selected) < N_parents:
            N_pop   = len(pop)
            s       = [0] * (N_pop+1)

            for i in range(1, N_pop+1):
                s[i] = s[i-1] + p(i, N_pop)

            r       = ptr_genetic_algorithm.RNG.random()
            index   = next(i-1 for i in range(len(s)) if r < s[i])

            selected.append(pop[index])
            del pop[index]

        return selected

class P_mutateFunctions:
    @staticmethod
    def default_P_mutate_function(sequence_length, cur_iter_percent):
        return 1.0/float(sequence_length)

    @staticmethod
    def exponential_decay(sequence_length, cur_iter_percent):
        """Very large mutation rate early, will decay to zero eventually."""
        return (1.0 - cur_iter_percent)**2



################################
## Core and default functions ##
################################

class GeneticAlgorithm:
    @staticmethod
    def __default_scoring_function(sequence): # Use this function signature
        """Replace this function with your own fitness function."""

        test_string = "PANTSERAFWEERGESCHUT"

        score = 0.0
        for i in range(len(sequence)):
            if sequence[i] == test_string[i]:
                score += 1
        return score

    __default_genes = (
        'A', 'C', 'D', 'E', 'F',
        'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R',
        'S', 'T', 'V', 'W', 'Y'
    )

    @staticmethod
    def __default_diversity_function(sequence, population, cur_iter_percent): # Use this function signature
        """Replace this function with your own custom function."""
        #NOTE: It is possible to include a scaling function as well, using cur_iter_percent (current_iteration/total_iterations). 
        # This allows for large diversity in early iterations while still allowing the algorithm to converge later on.
        # Default function -> No effect (i.e. score == fitness)
        return 0.0

    @staticmethod
    def __default_constraint_function(sequence): # Use this function signature
        """
        Replace this function with your own custom constraint function.
        
        Should return True if the sequence is not allowed, 
        preventing the sequence from being added to the population.
        
        Should return False otherwise.
        Default=no constraints
        """
        if len(sequence) < UserSettings.min_length:
            return True
        elif len(sequence) > UserSettings.max_length:
            return True

        # martinize2 does not know how to handle P at the termini
        # this constraint should always be active
        if sequence[0] == "P":
            #print("caught Proline")
            return True
        if sequence[-1] == "P":
            #print("caught Proline")
            return True

        return False

    def __init__(self):
        # Default values, set manually afterwards #
        self.iterations                 = 1000
        self.population_size            = 20
        self.sequence_length            = 20
        self.P_crossover                = 1
        self.N_parents                  = 5
        self.N_elites                   = 1 # N best invididuals that get copied directly to next generation
        self.N_rerun_elites             = 0 # N best individuals (that have #reruns > N_rerun_threshold) that get copied directly to next generation
        self.N_rerun_threshold          = 2
        self.ptr_scoring_function       = GeneticAlgorithm.__default_scoring_function
        self.ptr_selection_function     = SelectionMethods.best_N
        self.ptr_diversity_function     = GeneticAlgorithm.__default_diversity_function
        self.ptr_constraint_function    = GeneticAlgorithm.__default_constraint_function
        self.selfun_tournament_size     = 10
        self.selfun_exponential_rank_C  = 0.5
        self.ptr_P_mutate_function      = P_mutateFunctions.default_P_mutate_function
        self.genes                      = GeneticAlgorithm.__default_genes
        self.RNG                        = Random()
        self.filepath_progress_output   = None
        self.filepath_checkpoint        = None
        self.current_iteration          = 0
        # Set by read_checkpoint_file; the run loop continues the numbering from it.
        self.restart_from_iteration     = 0
        self.additional_output          = {}

        # Other values
        self._bRestartRun               = False
        self.description                = "Genetic Algorithm"
    
    def __check_settings(self): #TODO cleanup + include newer settings.
        def error(name, value, msg):
            mprint("Error:\t" + name + " (" + str(value) + ") " + msg.strip())
            sys.exit(1)
        
        if self.N_parents > self.population_size:
            error("N_parents", self.N_parents, "should be between 2 and population_size.")
        elif self.N_parents < 2:
            error("N_parents", self.N_parents, "should be equal to or larger than 2.")
        if self.population_size < 2:
            error("population_size", self.population_size, "should be equal to or larger than 2.")
        if self.ptr_selection_function is SelectionMethods.tournament:
            mprint("Tournament Selection enabled.")
            if self.selfun_tournament_size <= 0 or not self.selfun_tournament_size < (self.population_size - self.N_parents):
                error("selfun_tournament_size", self.selfun_tournament_size, "should be positive and smaller than population_size - N_parents.")

    def __write_progress_file(self, current_iteration, output_dict, n_run_dict):
        if MPI is not None: # Ensure only master writes
            if not MPI.COMM_WORLD.Get_rank() == 0:
                return

        def header_string(header):
            return header + ":"
        def body_string(key, value):
            return "\t" + str(key) + "\t" + str(value) + ""
        def footer_string():
            return ""

        def write_to_new_file(data):
            if self.filepath_progress_output is None:
                return

            # 'at', not 'wt'. run() always writes this header, including on a restart, and the
            # documented restart invocation points --output_dir and --restart at the same
            # directory. Truncating here destroyed the entire per-generation history of the
            # campaign being continued: the checkpoint survived, so the run resumed correctly
            # and silently lost the only human-readable record of how it got there, which is
            # also what tools/ga_report.py reads.
            resuming = os.path.isfile(self.filepath_progress_output)
            with open(self.filepath_progress_output, 'at' if resuming else 'wt') as file:
                if resuming:
                    # One marker line, not the whole settings banner again. The banner is parsed
                    # as settings only before the first generation; written a second time it
                    # lands in the middle of the run, where tools/ga_report.py counts every one
                    # of its lines as a row that failed to parse and reports that the format and
                    # the tool have diverged. The line below ends in a colon, which that parser
                    # already treats as a header rather than a row.
                    file.write("\n")
                    file.write("### resumed from checkpoint at iteration %d:\n"
                               % getattr(self, "restart_from_iteration", 0))
                    mprint("Appending to the existing progress file %s; the previous run's "
                           "history is kept above this point."
                           % os.path.abspath(self.filepath_progress_output))
                    return
                for line in data:
                    print(line)
                    file.write(line + "\n")

        def append_to_file(data):
            if self.filepath_progress_output is None:
                return
                
            with open(self.filepath_progress_output, 'at') as file:
                for line in data:
                    file.write(line + "\n")

        def create_progress_file():
            data = []
            data.append("###########")
            data.append("## EVOMD ##")
            data.append("###########")
            data.append("")

            data.append(header_string("Date"))
            data.append(datetime.now().strftime("%Y/%m/%d %H:%M:%S"))
            data.append("")

            data.append(header_string("Description"))
            description = self.description.split("\n")
            for line in description:
                data.append("\t" + line)
            data.append("")
            
            data.append(header_string("Settings"))

            if not UserSettings.propensities_initial:
                data.append(body_string("Genes:\t\t\t\t", self.genes))
            else:
                data.append(body_string("Genes:\t\t\t\t", [x for x in UserSettings.propensities.keys() if UserSettings.propensities[x] > 0.0]))
            data.append(body_string("Iterations:\t\t\t", self.iterations))
            data.append(body_string("Population Size:\t\t", self.population_size))
            data.append(body_string("Number of parents:\t\t", self.N_parents))
            data.append(body_string("Number of elites:\t\t", self.N_elites))
            data.append(body_string("Number of rerun elites:\t\t", self.N_rerun_elites))
            data.append(body_string("Rerun elites threshold:\t\t", self.N_rerun_threshold))
            data.append(body_string("Scoring function:\t\t", self.ptr_scoring_function.__name__))
            data.append(body_string("Selection method:\t\t", self.ptr_selection_function.__name__))
            data.append(body_string("Diversity function:\t\t", self.ptr_diversity_function.__name__))
            data.append(body_string("Constraint function:\t\t", self.ptr_constraint_function.__name__))
            data.append(body_string("P(crossover):\t\t\t", self.P_crossover))
            data.append(body_string("P_mutate function:\t\t", self.ptr_P_mutate_function.__name__))
            data.append(body_string("Tournament size (if applicable):", self.selfun_tournament_size))
            data.append(body_string("Exp_rank param C (if applicable):", self.selfun_exponential_rank_C))
            if self.filepath_checkpoint is not None:
                data.append(body_string("Filepath checkpoint:\t\t", self.filepath_checkpoint))
            for key, val in self.additional_output.items():
                data.append(body_string(key, val))
            
            data.append(footer_string())

            write_to_new_file(data)

        def write_iteration():
            sorted_data = sorted(output_dict.items(), key=lambda x: x[1][0], reverse=UserSettings.selection_sorting)
            data        = []

            data.append(header_string("Iteration " + str(current_iteration)))
            data.append(header_string("finished"))
            data.append(datetime.now().strftime("%Y/%m/%d %H:%M:%S"))

            for index in sorted_data:
                fitness, ddF = index[1][0], index[1][1]
                # ddF is exactly 0.0 for a failed evaluation (run.py returns [WORST_FITNESS, 0.0])
                # and can be 0.0 for a candidate that simply did not perturb the membrane. This
                # division used to be unguarded, which turned one failed candidate into a
                # ZeroDivisionError here, in the progress writer, after the whole generation had
                # already been simulated. That defeated the entire failure-tolerance design in
                # run.py and discarded a generation of GROMACS time, with a traceback pointing at
                # logging code that looks unrelated to the failure.
                penalty = (fitness / ddF) if ddF else float("nan")
                data.append(body_string(index[0], str(fitness) + "\tddF=" + str(ddF) + "\tpenalty=" + str(penalty) + "\tn_run=" + str(n_run_dict[index[0]]) ))

            data.append(footer_string())

            append_to_file(data)

        if current_iteration == -1:
            create_progress_file()
        else:
            write_iteration()

    def initialize(self, path_restart_file=None):
        """Confirms settings. Prepares genetic algorithm"""

        mprint("Checking settings and initializing the genetic algorithm..")

        self.__check_settings()

        self._population    = _Population(self)
        if path_restart_file is not None:
            self._bRestartRun = True
            self._population.read_checkpoint_file(path_restart_file)
    
    def run(self):
        mprint("Starting run:\n")

        self.__write_progress_file(-1, None, None)

        size = 1
        if MPI is not None:
            size = MPI.COMM_WORLD.Get_size()
            if size > 1:
                mprint("Simulations will be distributed over %d MPI threads.\n" % (size))

        # Short fix when running a restart/continuation using a checkpoint file.
        if self._bRestartRun: 
            self._population.generate_new_population(self.N_parents)

        # Continue the numbering rather than restarting it. The progress file is appended to on
        # a restart, and tools/ga_report.py keys generations by this number, so repeating it
        # would merge the resumed generations into the ones they continue from.
        offset = getattr(self, "restart_from_iteration", 0) + 1 if self._bRestartRun else 0

        for i in range(self.iterations):
            self.current_iteration = offset + i
            mprint("Status: Iteration %d/%d, %d simulations per iteration." % (
                self.current_iteration, 
                offset + self.iterations, 
                self.population_size, 
                )
            )

            output_dict = self._population.evaluate_fitness()
            # current_iteration, not the loop variable: on a restart they differ by the
            # offset, and the loop variable would relabel the resumed generations 0, 1, 2
            # again, which is what collided in the progress file.
            self.__write_progress_file(self.current_iteration, output_dict,
                                       self._population.n_run_dict)
            self._population.generate_new_population(self.N_parents)

    def get_result(self, N_output=1):
        """
        Returns:
            (list(tuple(string, float))): #'N_output' tuples containing the chromosome in string format and corresponding fitness.
        """

        return self._population.select_final(N_output)

