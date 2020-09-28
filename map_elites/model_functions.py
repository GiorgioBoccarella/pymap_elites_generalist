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
    #assert(len(sequences) >= k)
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
        #print("No mutation ")

    return mut_genome


def make_ind(t):
    g, trj, f, f1, f2, pos = t
    return cm.Ind(g, trj, f, f1, f2, pos)


def make_ind_inv(t):
    pos, g,  mod, trade, f, inv_pot, sum_p0, sum_p1 = t
    return cm.Ind(pos, g, mod, trade, f, inv_pot, sum_p0, sum_p1)


def add_to_archive(pos, ind, archive):
    archive[pos] = ind
    return 1


def env_pair_fitness(genome, env_pair):
        """Returns average fitness of the two environment"""

        genome = genome.T

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
        fit_m_1 = -math.sqrt(resnorm)

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
        fit_m_2 = -math.sqrt(resnorm)

        average_fitness = (fit_m_1 + fit_m_2)/2

        return average_fitness


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


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)


def score_Modularity(s_genome, all_g, env):
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
        l_a.append([fit_diff_A])
        l_b.append([fit_diff_B])


    sum_fitA = np.array([item for sublist in l_a for item in sublist])
    sum_fitB = np.array([item for sublist in l_b for item in sublist])

    rho, theta = cart2pol(sum_fitA, sum_fitB)

    modularity = np.mean(1 - np.divide(rho*abs(np.sin(2*theta)), rho, out=np.zeros_like(rho*abs(np.sin(2*theta))), where=rho != 0))

    return modularity


def sum_p_genome(genome):

    sum_p_0 = sum(genome[0])
    sum_p_1 = sum(genome[1])

    return sum_p_0, sum_p_1


def compute_versatility(max_evals=10,
            k=0,
            env_pair_dict_l=[],
            seq_list=[],
            sim=[],
            params=[]):

    print(params)
    cm.save_params(params)

    seq_list = [i for i in seq_list if sum(i) == params["p1"]]

    for sim_n in range(0, sim):
        print("Sim: ", sim_n)
        #The environment pairs are assigned from outside compute function
        env_pair_dict = env_pair_dict_l[sim_n]
        cm.save_env(env_pair_dict, sim_n)

        #This is the random seed for mutations
        np.random.seed(sim_n + params["seed"] + 4)

        # init archive (empty)
        archive = {}
        init = 0

        # main loop
        n_evals = 0  # number of evaluations
        while (n_evals < max_evals + 1):
            # If environment is empty fill with random individuals
            if init == 0:
                # Starting genome is the same for all env_pair
                g = generate_genome(seq_list, k)
                for i in env_pair_dict.keys():
                    # Generate fitness
                    fit = env_pair_fitness(g, env_pair_dict[i])
                    #Genrate invasion_potential
                    inv_pot = {}
                    for e in env_pair_dict:
                        fit = env_pair_fitness(g, env_pair_dict[e])
                        inv_pot[e] = [fit]
                    #Score starting feature of the genome gen 0
                    all_mut = generate_all_mutants(g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(g, all_mut, env)
                    score_mod = score_Modularity(g, all_mut, env)
                    pos = i
                    # Pack traits and feature, make Ind and add to archive
                    sum_p0, sum_p1 = sum_p_genome(g)
                    *to_gen_ind, = pos, g, score_mod, score_tradeOff, fit, inv_pot, sum_p0, sum_p1,
                    ind = make_ind_inv(to_gen_ind)
                    add_to_archive(pos, ind, archive)
                    init = 1
            else:
                invasion = np.random.binomial(1, params['invasion_rate'], 1)
                if invasion == True and sim_n > (sim/2) - 1:
                    for ind in archive.keys():
                        keys = list(archive.keys())
                        coordinates = np.random.choice(keys, 1, replace=False)
                        invader = archive[ind]
                        wild_type = archive[coordinates[0]]


                        inv_f = invader.invasion_potential[coordinates[0]]
                        w_f = wild_type.invasion_potential[coordinates[0]]
                        env_inv = 0

                        if inv_f[env_inv] > w_f[env_inv]:
                            archive[coordinates[0]] = copy.deepcopy(archive[ind])
                            archive[coordinates[0]].position = coordinates[0]
                            #If invasion is succesfull then save it
                            cm.__save_file_mig(ind, coordinates[0], n_evals, sim_n, params['invasion_rate'], params["p1"], invader.sum_p0, invader.sum_p1, wild_type.sum_p0, wild_type.sum_p1)
            if n_evals < max_evals + 1:
                for i in archive.keys():
                    start_g = archive[i].genome
                    all_mut = generate_all_mutants(start_g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                    score_mod = score_Modularity(start_g, all_mut, env)
                    mutation = np.random.binomial(1, params['mutation_rate'], 1)
                    if mutation == True:
                        mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                        if mutated_genome != []:
                            archive[i].genome = mutated_genome
                    archive[i].fitness = env_pair_fitness(archive[i].genome, env)
                    #Here calculate invasion potential
                    for e in env_pair_dict.keys():
                        fit = env_pair_fitness(archive[i].genome, env_pair_dict[e])
                        archive[i].invasion_potential[e] = [fit]
                    #Save the trade_off and modularity score
                    archive[i].trade_off = score_tradeOff
                    archive[i].modularity = score_mod
                    archive[i].sum_p0, archive[i].sum_p1 = sum_p_genome(archive[i].genome)
                    inv = 0
                    if sim_n > (sim / 2) - 1:
                        inv = 1
            cm.__save_archive(archive, n_evals, sim_n, 0, params['invasion_rate'], inv, params["p1"])
            n_evals += 1
            print("Sim: ", sim_n, "Step: ", n_evals)

        archive_c = archive.copy()

        for transf in params["env_transfer"]:
            archive = archive_c.copy()
            n_evals = max_evals + 1
            max = n_evals + 50
            while n_evals < max:
                for i in archive.keys():
                    start_g = archive[i].genome
                    all_mut = generate_all_mutants(start_g)
                    env = env_pair_dict[transf]
                    score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                    score_mod = score_Modularity(start_g, all_mut, env)
                    mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                    if mutated_genome != []:
                        archive[i].genome = mutated_genome
                    archive[i].fitness = env_pair_fitness(archive[i].genome, env)
                    # Here calculate invasion potential
                    for e in env_pair_dict.keys():
                        fit = env_pair_fitness(archive[i].genome, env_pair_dict[e])
                        archive[i].invasion_potential[e] = [fit]
                    # Save the trade_off and modularity score
                    archive[i].trade_off = score_tradeOff
                    archive[i].modularity = score_mod
                    archive[i].sum_p0, archive[i].sum_p1 = sum_p_genome(archive[i].genome)
                cm.__save_archive(archive, n_evals, sim_n, transf, params['invasion_rate'], inv, params["p1"])
                n_evals += 1
                print("Sim: ", sim_n, "Step_replicate: ", n_evals)
    return archive







