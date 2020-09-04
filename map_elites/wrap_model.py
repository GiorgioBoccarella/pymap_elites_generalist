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

from map_elites import common_new as cm
from examples import lsq_f


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

    return all_mutants_array


def gen_lucky_mut(s_genome, all_g, env):
    """This generates a lucky mutation that is beneficial in at least one environment"""

    # Which env?
    env = env[np.random.binomial(1, 0.5, 1)]
    env = env.flatten()
    # print(s_genome)

    # Calculate Fitness of the starting genome
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env)
    fit_s = -math.sqrt(resnorm)

    print("Fitness before mutation:")
    print(fit_s)
    print()

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
        print("Runned out of beneficial mutations")

    return mut_genome



def make_ind(t):
    g, trj, f, f1, f2, pos = t
    return cm.Ind(g, trj, f, f1, f2, pos)

def add_to_archive(ind, archive):
    archive[ind.position] = ind
    return 1

def env_pair_fitness(genome, env_pair):
        genome = genome.T

        print(env_pair)

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
        fit_m_1 = -math.sqrt(resnorm)

        print("Fit in E1:")
        print(fit_m_1)

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
        fit_m_2 = -math.sqrt(resnorm)

        print("Fit in E2:")
        print(fit_m_2)
        print()

        #print("Modularity:")
        #mod = print(fit_m_1/fit_m_2)
        #print(mod)
        #print()

        average_fitness = (fit_m_1 + fit_m_2)/2

        return average_fitness, fit_m_1, fit_m_2


def scoreTradeOff(s_genome, all_g, env):

    # Calculate Fitness of the starting genome in A
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env[0])
    fit_A = -math.sqrt(resnorm)

    print("Fit A:")
    print(fit_A)
    print()

    # Calculate Fitness of the starting genome in B
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env[1])
    fit_B = -math.sqrt(resnorm)

    print("Fit B:")
    print(fit_B)
    print()

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
            print()


    if l_a == []:
        exit()
        print("Runned out of beneficial mutations")

    sum_fitA = np.array([item for sublist in l_a for item in sublist])
    sum_fitB = np.array([item for sublist in l_b for item in sublist])

    mut_tradeOff = - statistics.mean(sum_fitA*sum_fitB)/math.sqrt(statistics.mean(sum_fitA**2) * statistics.mean(sum_fitB**2))

    print(mut_tradeOff)
    return mut_tradeOff

def compute(max_evals=1e3,
            k=0,
            env_pair_dict=[],
            seq_list=[],
            params=cm.default_params):
     
    print(params)

    assert (len(seq_list) >= k)

    # init archive (empty)
    archive = {}

    init = 0

    # main loop
    n_evals = 0 # number of evaluations
    # TOD count the successes
    while (n_evals < 400):
        #If environment is empty fill with random individuals
        if init == 0:
            for i in env_pair_dict.keys():
                g = generate_genome(seq_list, k)
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
                print("scoreTradeOff")
                print(score_tradeOff)
                mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                if mutated_genome != []:
                    archive[i].genome = mutated_genome
                else:
                    print("Runned out of beneficial mutation")
                    print()
                    print(archive[i].genome)
                archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[i])
                print("Average fitness: ")
                print(archive[i].fitness)
                print()
            n_evals += 1
            print("Evaluation: ")
            print(n_evals, )
            print()
            cm.__save_archive(archive, n_evals)
    return archive






def compute_mut(archive, env_pair_dict, steps = 20):

    start_fit = env_pair_dict.copy()

    for i in start_fit.keys():
        start_fit[i] = -0.0001

    vec_fit = start_fit.copy()

    for i in vec_fit.keys():
        vec_fit[i] = []



    for s in range(0, steps):
        for i in archive.keys():
            fit_difference = 0
            start_g = archive[i].genome
            all_mut = generate_all_mutants(start_g)
            env = env_pair_dict[i]
            mutated_genome = gen_lucky_mut(start_g, all_mut, env)
            if mutated_genome != []:
                archive[i].genome = mutated_genome
            else:
                print("Runned out of beneficial mutation")
                print()
            archive[i].fitness, archive[i].fit1, archive[i].fit2 = env_pair_fitness(archive[i].genome, env_pair_dict[1])
            print("Step: ")
            fit_difference = -(start_fit[i] - archive[i].fitness)
            vec_fit[i].append(fit_difference)
            i
        s = s + 1
        print(s)
        print()

    print(vec_fit)
    return archive

