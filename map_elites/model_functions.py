#! /usr/bin/env python
#| This file is a part of the pymap_elites framework.
#| Copyright 2019, INRIA
#| Main contributor(s):
#| Jean-Baptiste Mouret, jean-baptiste.mouret@inria.fr
#| Eloise Dalin , eloise.dalin@inria.fr
#| Pierre Desreumaux , pierre.desreumaux@inria.fr
#|
#|
#| **Main paper**: Mouret JB, Clune J. Illuminating search spaces by
#| mapping elites. arXiv preprint arXiv:1504.04909. 2015 Apr 20.
#|
#| This software is governed by the CeCILL license under French law
#| and abiding by the rules of distribution of free software.  You
#| can use, modify and/ or redistribute the software under the terms
#| of the CeCILL license as circulated by CEA, CNRS and INRIA at the
#| following URL "http://www.cecill.info".
#|
#| As a counterpart to the access to the source code and rights to
#| copy, modify and redistribute granted by the license, users are
#| provided only with a limited warranty and the software's author,
#| the holder of the economic rights, and the successive licensors
#| have only limited liability.
#|
#| In this respect, the user's attention is drawn to the risks
#| associated with loading, using, modifying and/or developing or
#| reproducing the software by the user in light of its specific
#| status of free software, that may mean that it is complicated to
#| manipulate, and that also therefore means that it is reserved for
#| developers and experienced professionals having in-depth computer
#| knowledge. Users are therefore encouraged to load and test the
#| software's suitability as regards their requirements in conditions
#| enabling the security of their systems and/or data to be ensured
#| and, more generally, to use and operate it in the same conditions
#| as regards security.
#|
#| The fact that you are presently reading this means that you have
#| had knowledge of the CeCILL license and that you accept its terms.
#
# from scipy.spatial import cKDTree : TODO

import math
import numpy as np
from numpy.random import choice
import statistics
from collections import defaultdict

from map_elites import common_invasion as cm
from examples import lsq_f

import copy



def generate_genome(sequences, k):
    """Generates genome with certain K size, the L depends on the length of the sequences"""
    assert(len(sequences) >= k)
    l = len(sequences[0])
    g = np.empty([k, l], dtype=bool)
    sel_id = np.random.choice(len(sequences), k, replace=False)
    j = 0
    for seq in sel_id:
        g[j] = sequences[seq]
        j += 1
    return g


def generate_all_mutants(genome):
    """This generate all the mutation of the original genome, all the K_i is included"""

    all_mut_e = []

    for i in range(0, len(genome)):
        s_genome = genome[i]
        all_mut = np.tile(s_genome, (len(s_genome) - 1, 1))
        for t in range(0, len(s_genome) - 1):
            all_mut[t][t] ^= 1
        g_mut = np.repeat(genome[np.newaxis, ...], len(all_mut), axis=0)
        for j in range(0, len(all_mut)):
            g_mut[j][i] = all_mut[j]
        all_mut_e.append(g_mut)

    all_mutants_array = np.vstack(all_mut_e)

    print()
    return all_mutants_array


def gen_lucky_mut(s_genome, all_g, env):
    """This generates a lucky mutation that is beneficial in at least one environment"""

    # Which env?
    ran = np.random.binomial(1, 0.5, 1)
    env = env[ran]
    env = env.flatten()
    # print(s_genome)

    # Calculate Fitness of the starting genome
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env)
    fit_s = -math.sqrt(resnorm)


    l = []
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env)
        fit_m = -math.sqrt(resnorm)
        fit_diff = fit_m - fit_s
        if fit_diff > 0:
            l.append([fit_diff, i])

    if l != []:
        # Fitness increase in env becomes probability that sum up to one
        p = np.array([fitness[0] for fitness in l])
        p /= p.sum()

        # Genomes trait value
        genomes = np.array([genome[1] for genome in l])

        # Here lucky mutation is selected
        draw_n = choice(len(genomes), 1, p=p)
        mut_genome = genomes[draw_n]
        mut_genome = mut_genome.reshape(s_genome.shape)
    else:
        mut_genome = s_genome
        print("No mutation ")

    return mut_genome



def make_ind(t):
    g, trj, f, f1, f2, pos = t
    return cm.Ind(g, trj, f, f1, f2, pos)

def make_ind_inv(t):
    g, trj, f, f1, f2, inv_pot, pos = t
    return cm.Ind(g, trj, f, f1, f2, inv_pot, pos)

def add_to_archive(ind, archive):
    archive[ind.position] = ind
    return 1

def env_pair_fitness(genome, env_pair):
        """Returns average fitness of the two environment"""

        genome = genome.T

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
        fit_m_1 = -math.sqrt(resnorm)

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
        fit_m_2 = -math.sqrt(resnorm)

        average_fitness = (fit_m_1 + fit_m_2)/2

        return average_fitness, fit_m_1, fit_m_2


def scoreTradeOff(s_genome, all_g, env):
    """Calculate trade-off (How much a mutation affects the fitness in the environment pair)"""

    # Calculate Fitness of the starting genome in A
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env[0])
    fit_A = -math.sqrt(resnorm)

    # Calculate Fitness of the starting genome in B
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env[1])
    fit_B = -math.sqrt(resnorm)

    # Calculate fitness of all mutant genomes
    # fit_vec = np.empty([len(all_g)], dtype = float)

    l_a = []
    l_b = []
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env[0])
        fit_m_A = -math.sqrt(resnorm)

        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env[1])
        fit_m_B = -math.sqrt(resnorm)

        fit_diff_A = fit_m_A - fit_A
        fit_diff_B = fit_m_B - fit_B
        if fit_diff_A > 0 or fit_diff_B > 0:
            l_a.append([fit_diff_A])
            l_b.append([fit_diff_B])


    if l_a == []:
        print("Runned out of beneficial mutations score")
        #exit()
        mut_tradeOff = 0
    else:
        sum_fitA = np.array([item for sublist in l_a for item in sublist])
        sum_fitB = np.array([item for sublist in l_b for item in sublist])

        mut_tradeOff = - statistics.mean(sum_fitA*sum_fitB)/math.sqrt(statistics.mean(sum_fitA**2) * statistics.mean(sum_fitB**2))

    #print(mut_tradeOff)
    return mut_tradeOff



def compute(max_evals=10,
            k=0,
            env_pair_dict=[],
            seq_list=[],
            sim=[],
            params=cm.default_params):

    print(params)

    for sim_n in range(0, sim):
        assert (len(seq_list) >= k)



        # init archive (empty)
        archive = {}
        init = 0

        # main loop
        n_evals = 0 # number of evaluations
        # TOD count the successes
        while (n_evals < max_evals):
            # If environment is empty fill with random individuals
            if init == 0:
                # Starting genome is the same for all env_pair
                g = generate_genome(seq_list, k)
                for i in env_pair_dict.keys():
                    # Trajectory is initialized with first environment
                    # Later new positions are appended
                    trj = i
                    # Same for position but this is the actual position
                    pos = i
                    # Generate fitness
                    fit, fitA, fitB = env_pair_fitness(g, env_pair_dict[i])
                    # Pack traits and feature, make Ind and add to archive
                    *to_gen_ind, = g, trj, fit, fitA, fitB, pos
                    ind = make_ind(to_gen_ind)
                    add_to_archive(ind, archive)
                    init = 1
            else:
                for i in archive.keys():
                    start_g = archive[i].genome
                    all_mut = generate_all_mutants(start_g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                    mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                    if mutated_genome != []:
                        archive[i].genome = mutated_genome
                    else:
                        print("Runned out of beneficial mutation")
                        print()
                        print(archive[i].genome)
                    archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[i])
                    cm.__save_file(score_tradeOff, i, n_evals, sim_n)
                n_evals += 1
                print("Evaluation: ")
                print(n_evals)
                cm.__save_archive(archive, n_evals, sim_n)
    return archive






def compute_mut(archive, env_pair_dict, steps = 20):

    #We test starting fitness in mid environement before starting evo
    for i in archive.keys():
        archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[0.9])

    start_fit = {}

    for i in env_pair_dict.keys():
        start_fit[i] = [archive[i].fitness]

    vec_fit = start_fit.copy()

    for i in vec_fit.keys():
        vec_fit[i] = []

    vec_toff = start_fit.copy()

    for i in vec_toff.keys():
        vec_toff[i] = []


    for s in range(0, steps):
        for i in archive.keys():
            fit_difference = 0
            start_g = archive[i].genome
            all_mut = generate_all_mutants(start_g)
            env = env_pair_dict[0.9]
            score_tradeOff = scoreTradeOff(start_g, all_mut, env)
            mutated_genome = gen_lucky_mut(start_g, all_mut, env)
            if mutated_genome != []:
                archive[i].genome = mutated_genome
            else:
                print("Runned out of beneficial mutation")
                print()
            archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[0.9])
            print("Step: ")
            #fit_difference = - (np.asarray(start_fit[i]) - archive[i].fitness)
            fit_difference = archive[i].fitness
            vec_fit[i].append(float(fit_difference))
            vec_toff[i].append(score_tradeOff)
        s = s + 1
        print(s)
        print()

    print("Fitness of the mutants:")
    print(vec_fit)
    cm.__save_file_mut(vec_fit,vec_toff)
    return archive


def compute_invasion(max_evals=10,
            k=0,
            env_pair_dict=[],
            seq_list=[],
            sim=[],
            params=cm.default_params):

    print(params)

    for sim_n in range(0, sim):
        assert (len(seq_list) >= k)

        np.random.seed(sim_n + params["seed"] * 3)

        # init archive (empty)
        archive = {}
        init = 0

        # main loop
        n_evals = 0 # number of evaluations
        # TOD count the successes
        while (n_evals < max_evals):
            # If environment is empty fill with random individuals
            if init == 0:
                # Starting genome is the same for all env_pair
                g = generate_genome(seq_list, k)
                for i in env_pair_dict.keys():
                    # Trajectory is initialized with first environment
                    # Later new positions are appended
                    trj = [i]
                    # Same for position but this is the actual position
                    pos = i
                    # Generate fitness
                    fit, fitA, fitB = env_pair_fitness(g, env_pair_dict[i])
                    #Genrate invasion_potential
                    inv_pot = {}
                    for e in env_pair_dict:
                        fit, fitA, fitB = env_pair_fitness(g, env_pair_dict[e])
                        inv_pot[e] = [fit, fitA, fitB]
                    # Pack traits and feature, make Ind and add to archive
                    *to_gen_ind, = g, trj, fit, fitA, fitB, inv_pot, pos
                    ind = make_ind_inv(to_gen_ind)
                    add_to_archive(ind, archive)
                    init = 1
            else:
                invasion = np.random.binomial(1, 0.05, 1)
                if invasion == True and n_evals > 20:
                    keys = list(archive.keys())
                    print()
                    coordinates = np.random.choice(keys, 2, replace=False)
                    invader = archive[coordinates[0]]
                    wild_type = archive[coordinates[1]]

                    # During which environmental exposure is the invasion happening?
                    env = [1, 2]
                    env_inv = int(np.random.choice(env, 1, replace=False))

                    inv_f = invader.invasion_potential[coordinates[0]]
                    w_f = wild_type.invasion_potential[coordinates[0]]
                    print()

                    env_inv = 0

                    if inv_f[env_inv] > w_f[env_inv]:
                        archive[coordinates[1]] = copy.deepcopy(archive[coordinates[0]])
                        archive[coordinates[1]].trajectory.append(coordinates[1])
                        cm.__save_file_mig(coordinates[0], coordinates[1], n_evals, sim_n)
                else:
                    for i in archive.keys():
                        start_g = archive[i].genome
                        all_mut = generate_all_mutants(start_g)
                        env = env_pair_dict[i]
                        score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                        mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                        if mutated_genome != []:
                            archive[i].genome = mutated_genome
                        else:
                            print("Runned out of beneficial mutation")
                            print()
                            print(archive[i].genome)
                        archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[i])
                        for e in env_pair_dict.keys():
                            fit, fitA, fitB = env_pair_fitness(archive[i].genome, env_pair_dict[e])

                            archive[i].invasion_potential[e] = [fit, fitA, fitB]
                        cm.__save_file(score_tradeOff, i, n_evals, sim_n)
                n_evals += 1
                print("Evaluation: ")
                print(n_evals)
                cm.__save_archive(archive, n_evals, sim_n)
    return archive

# TODO problem because then I do not know where they came from

def compute_inv(max_evals=10,
            k=0,
            env_pair_dict=[],
            seq_list=[],
            sim=[],
            params=cm.default_params):

    print(params)

    for sim_n in range(0, sim):
        assert (len(seq_list) >= k)

        # init archive (empty)
        archive = {}
        init = 0

        # main loop
        n_evals = 0 # number of evaluations
        # TOD count the successes
        while (n_evals < max_evals):
            # If environment is empty fill with random individuals
            if init == 0:
                # Starting genome is the same for all env_pair
                g = generate_genome(seq_list, k)
                for i in env_pair_dict.keys():
                    # Trajectory is initialized with first environment
                    # Later new positions are appended
                    trj = i
                    # Same for position but this is the actual position
                    pos = i
                    # Generate fitness
                    fit, fitA, fitB = env_pair_fitness(g, env_pair_dict[i])
                    #Genrate invasion_potential
                    inv_pot = {}
                    for e in env_pair_dict:
                        fit, fitA, fitB = env_pair_fitness(g, env_pair_dict[e])
                        inv_pot[e] = [fit, fitA, fitB]
                    # Pack traits and feature, make Ind and add to archive
                    *to_gen_ind, = g, trj, fit, fitA, fitB, inv_pot, pos
                    ind = make_ind_inv(to_gen_ind)
                    add_to_archive(ind, archive)
                    init = 1
            else:
                for i in archive.keys():
                    start_g = archive[i].genome
                    all_mut = generate_all_mutants(start_g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                    mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                    if mutated_genome != []:
                        archive[i].genome = mutated_genome
                    else:
                        print("Runned out of beneficial mutation")
                        print()
                        print(archive[i].genome)
                    archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[i])
                    cm.__save_file(score_tradeOff, i, n_evals, sim_n)
                n_evals += 1
                print("Evaluation: ")
                print(n_evals)
                cm.__save_archive(archive, n_evals, sim_n)
    return archive
