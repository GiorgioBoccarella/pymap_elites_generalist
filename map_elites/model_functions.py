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
import random
from collections import defaultdict

from map_elites import common_invasion as cm
from examples import lsq_f

import copy

def generate_genome(k, l, p, seed):
    """Generates genome with certain K size, the L depends on the length of the sequences"""
    random.seed(seed)
    a = (k, l)
    g = np.zeros(a, dtype=bool)

    for i in range(0, len(g)):
        p1 = random.sample(range(l), p)
        for j in range(0, len(p1)):
            g[i][p1[j]] ^= 1

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


def make_ind_inv(t):
    pos, g,  mod, trade, f, inv_pot, sum_p0, sum_p1, ps = t
    return cm.Ind(pos, g, mod, trade, f, inv_pot, sum_p0, sum_p1, ps)


def add_to_archive(pos, ind, archive):
    archive[pos] = ind
    return 1


def sum_p_genome(genome):

    sum_p_0 = sum(genome[0])
    sum_p_1 = sum(genome[1])

    return sum_p_0, sum_p_1


def env_pair_fitness(genome, env_pair):
        """Returns average fitness of the two environment"""

        genome = genome.T

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
        fit_m_1 = -math.sqrt(resnorm)

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
        fit_m_2 = -math.sqrt(resnorm)

        average_fitness = (fit_m_1 + fit_m_2)/2

        return average_fitness, fit_m_1, fit_m_2


def env_pair_fitness_average(genome, env_pair):
    """Returns average fitness of the two environment"""

    genome = genome.T

    [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
    fit_m_1 = -math.sqrt(resnorm)

    [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
    fit_m_2 = -math.sqrt(resnorm)

    average_fitness = (fit_m_1 + fit_m_2) / 2

    return average_fitness


def env_pair_fitness_AB(genome, env_pair):
    """Returns average fitness of the two environment"""

    genome = genome.T

    [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
    fit_m_1 = -math.sqrt(resnorm)

    [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
    fit_m_2 = -math.sqrt(resnorm)

    average_fitness = (fit_m_1 + fit_m_2) / 2

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


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    return rho, phi


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


def compute_invasion_transfer(max_evals=10,
            k=0,
            env_pair_dict_l=[],
            seq_list=[],
            sim=[],
            params=[]):

    print(params)
    cm.save_params(params)

    for sim_n in range(0, sim):
        print("Sim: ", sim_n)
        # The environment pairs are assigned from outside compute function
        env_pair_dict = env_pair_dict_l[sim_n]
        env_pair_dict_t = copy.deepcopy(env_pair_dict)
        del env_pair_dict[params["del"]]

        cm.save_env(env_pair_dict_t, sim_n)

        # This is the random seed for mutations
        np.random.seed(sim_n + params["seed"])

        # init archive (empty)
        archive = {}
        init = 0

        transf_n = 0

        # main loop
        n_evals = 0  # number of evaluations
        while n_evals < max_evals + 1:
            # If environment is empty fill with random individuals
            if init == 0:
                j = 0
                # p is the initial P1 in the genome for each lineage
                # if p is the same then modify to int ands not array
                for i in env_pair_dict.keys():
                    g = generate_genome(params['k'], params['l'], params['p'], sim_n + params["seed"])
                    # Generate fitness
                    fit = env_pair_fitness_average(g, env_pair_dict[i])
                    #Genrate invasion_potential
                    inv_pot = {}
                    for e in env_pair_dict:
                        fit, fit_1, fit_2 = env_pair_fitness(g, env_pair_dict[e])
                        inv_pot[e] = [fit, fit_1, fit_2]
                    #Score starting feature of the genome g
                    all_mut = generate_all_mutants(g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(g, all_mut, env)
                    score_mod = score_Modularity(g, all_mut, env)
                    pos = i
                    # Pack traits and feature, make Ind and add to archive
                    sum_p0, sum_p1 = sum_p_genome(g)
                    ps = sum_p0 + sum_p1
                    *to_gen_ind, = pos, g, score_mod, score_tradeOff, fit, inv_pot, sum_p0, sum_p1, ps
                    ind = make_ind_inv(to_gen_ind)
                    add_to_archive(pos, ind, archive)
                    j += 1
                init = 1
            # If we want to look at the effect of the mutations then params['invasion_rate'] > 0
            # Generally if 20 simulations I do 10 with invasion and 10 without to have a control
            invasion = np.random.binomial(1, params['invasion_rate'], 1)
            if invasion == True and sim_n > (sim / 2) - 1 and n_evals > 30:
                for ind in archive.keys():
                    keys = list(archive.keys())
                    coordinates = np.random.choice(keys, 1, replace=False)
                    w = coordinates[0]
                    invader = archive[ind]
                    wild_type = archive[w]

                    inv_f = invader.invasion_potential[w]
                    w_f = wild_type.invasion_potential[w]

                    e = [1, 2]
                    coordinates = np.random.choice(e, 1, replace=False)
                    env_inv = int(coordinates)

                    if inv_f[env_inv] > w_f[env_inv]:
                        archive[w] = copy.deepcopy(archive[ind])
                        archive[w].position = w
                        #If invasion is succesfull then save it
                        cm.__save_file_mig(ind, w, n_evals, sim_n, params['s_invasion'],
                                            invader.ps, wild_type.ps, invader.sum_p0, invader.sum_p1,
                                            wild_type.sum_p0, wild_type.sum_p1)
            if n_evals < max_evals + 1:
                for i in archive.keys():
                    start_g = archive[i].genome
                    all_mut = generate_all_mutants(start_g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                    score_mod = score_Modularity(start_g, all_mut, env)
                    mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                    if mutated_genome != []:
                        archive[i].genome = mutated_genome
                    archive[i].fitness = env_pair_fitness_average(archive[i].genome, env)
                    # Here calculate invasion potential (only if invasion can happen, otherwise is not used)
                    for e in env_pair_dict.keys():
                        fit, fit_1, fit_2 = env_pair_fitness(archive[i].genome, env_pair_dict[e])
                        archive[i].invasion_potential[e] = [fit, fit_1, fit_2]
                    #Save the trade_off and modularity score
                    archive[i].trade_off = score_tradeOff
                    archive[i].modularity = score_mod
                    archive[i].sum_p0, archive[i].sum_p1 = sum_p_genome(archive[i].genome)
                    inv = 0
                    if sim_n > (sim / 2) - 1:
                        inv = 1
            cm.__save_archive(archive, n_evals, sim_n, 0, 0, params['s_invasion'], inv, params["p"])
            n_evals += 1
            # At specific time point specified in params['transfer'] the linage are tested outside the environment
            # in which they are evolving in. (Useful for arrow plot)
            if n_evals in params['transfer']:
                transf_n += 1
                for transf in params["env_transfer"]:
                    archive_t = copy.deepcopy(archive)
                    n_evals_t = n_evals
                    max = n_evals_t + 100
                    n_evals_t = n_evals_t + 1
                    while n_evals_t < max:
                        for i in archive_t.keys():
                            start_g = archive_t[i].genome
                            all_mut = generate_all_mutants(start_g)
                            env = env_pair_dict_t[transf]
                            score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                            score_mod = score_Modularity(start_g, all_mut, env)
                            mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                            if mutated_genome != []:
                                archive_t[i].genome = mutated_genome
                            archive_t[i].fitness = env_pair_fitness_average(archive_t[i].genome, env)
                            # Here calculate invasion potential
                            for e in env_pair_dict.keys():
                                fit, fit_1, fit_2 = env_pair_fitness(archive_t[i].genome, env_pair_dict[e])
                                archive_t[i].invasion_potential[e] = [fit, fit_1, fit_2]
                            # Save the trade_off and modularity score
                            archive_t[i].trade_off = score_tradeOff
                            archive_t[i].modularity = score_mod
                            archive_t[i].sum_p0, archive_t[i].sum_p1 = sum_p_genome(archive_t[i].genome)
                        cm.__save_archive(archive_t, n_evals_t, sim_n, transf, transf_n, params['s_invasion'],
                                          inv, params["p"])
                        n_evals_t += 1
                        #print("Sim: ", sim_n, "Step_replicate: ", n_evals_t, " env_tras: ", transf)
    return archive


def compute_invasion_transfer_p_different(max_evals=10,
            k=0,
            env_pair_dict_l=[],
            seq_list=[],
            sim=[],
            params=[]):

    print(params)
    cm.save_params(params)

    for sim_n in range(0, sim):
        print("Sim: ", sim_n)
        # The environment pairs are assigned from outside compute function
        env_pair_dict = env_pair_dict_l[sim_n]
        cm.save_env(env_pair_dict, sim_n)

        # This is the random seed for mutations
        np.random.seed(sim_n + params["seed"] + 4)

        # init archive (empty)
        archive = {}
        init = 0

        # main loop
        n_evals = 0  # number of evaluations
        while n_evals < max_evals + 1:
            # If environment is empty fill with random individuals
            if init == 0:
                j = 0
                # p is the initial P1 in the genome for each lineage
                # if p is the same then modify to int ands not array
                p = np.array(params['p_list'])
                for i in env_pair_dict.keys():
                    g = generate_genome(seq_list, k, p[j])
                    # Generate fitness
                    fit = env_pair_fitness(g, env_pair_dict[i])
                    #Genrate invasion_potential
                    inv_pot = {}
                    for e in env_pair_dict:
                        fit = env_pair_fitness(g, env_pair_dict[e])
                        inv_pot[e] = [fit]
                    #Score starting feature of the genome g
                    all_mut = generate_all_mutants(g)
                    env = env_pair_dict[i]
                    score_tradeOff = scoreTradeOff(g, all_mut, env)
                    score_mod = score_Modularity(g, all_mut, env)
                    pos = i
                    # Pack traits and feature, make Ind and add to archive
                    sum_p0, sum_p1 = sum_p_genome(g)
                    ps = sum_p0
                    *to_gen_ind, = pos, g, score_mod, score_tradeOff, fit, inv_pot, sum_p0, sum_p1, ps
                    ind = make_ind_inv(to_gen_ind)
                    add_to_archive(pos, ind, archive)
                    j += 1
                init = 1
            # If we want to look at the effect of the mutations then params['invasion_rate'] > 0
            # Generally if 20 simulations I do 10 with invasion and 10 without to have a control
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

                        if inv_f[env_inv] > w_f[env_inv] and ind != 0.9 and coordinates[0] != 0.9:
                            archive[coordinates[0]] = copy.deepcopy(archive[ind])
                            archive[coordinates[0]].position = coordinates[0]
                            #If invasion is succesfull then save it
                            cm.__save_file_mig(ind, coordinates[0], n_evals, sim_n, params['invasion_rate'], invader.ps,
                            wild_type.ps, invader.sum_p0, invader.sum_p1, wild_type.sum_p0, wild_type.sum_p1)
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
                    # Here calculate invasion potential (only if invasion is activated, otherwise is not used)
                    if sim_n > (sim / 2) - 1:
                        for e in env_pair_dict.keys():
                            fit = env_pair_fitness(archive[i].genome, env_pair_dict[e])
                            archive[i].invasion_potential[e] = [fit]
                    #Save the trade_off and modularity score
                    archive[i].trade_off = score_tradeOff
                    archive[i].modularity = score_mod
                    archive[i].sum_p0, archive[i].sum_p1 = sum_p_genome(archive[i].genome)
                    inv = 0
                    if sim_n > (sim / 2) - 1:
                        # This is useful just for the data output and plotting
                        inv = 1
            cm.__save_archive(archive, n_evals, sim_n, 0, params['invasion_rate'], inv, params["p1"])
            n_evals += 1
            print("Sim: ", sim_n, "Step: ", n_evals)
            # At specific time point specified in params['transfer'] the linage are tested outside the environment
            # in which they are evolving in. (Useful for arrow plot)
            if n_evals in params['transfer']:
                for transf in params["env_transfer"]:
                    archive_t = copy.deepcopy(archive)
                    n_evals_t = n_evals
                    max = n_evals_t + 200
                    cm.__save_archive(archive_t, n_evals_t, sim_n, transf, params['invasion_rate'], inv, params["p1"])
                    n_evals_t = n_evals_t + 1
                    while n_evals_t < max:
                        for i in archive_t.keys():
                            start_g = archive_t[i].genome
                            all_mut = generate_all_mutants(start_g)
                            env = env_pair_dict[transf]
                            score_tradeOff = scoreTradeOff(start_g, all_mut, env)
                            score_mod = score_Modularity(start_g, all_mut, env)
                            mutated_genome = gen_lucky_mut(start_g, all_mut, env)
                            if mutated_genome != []:
                                archive_t[i].genome = mutated_genome
                            archive_t[i].fitness = env_pair_fitness(archive_t[i].genome, env)
                            # Here calculate invasion potential
                            for e in env_pair_dict.keys():
                                fit = env_pair_fitness(archive_t[i].genome, env_pair_dict[e])
                                archive_t[i].invasion_potential[e] = [fit]
                            # Save the trade_off and modularity score
                            archive_t[i].trade_off = score_tradeOff
                            archive_t[i].modularity = score_mod
                            archive_t[i].sum_p0, archive_t[i].sum_p1 = sum_p_genome(archive_t[i].genome)
                        cm.__save_archive(archive_t, n_evals_t, sim_n, transf, params['invasion_rate'], inv, params["p1"])
                        n_evals_t += 1
                        print("Sim: ", sim_n, "Step_replicate: ", n_evals_t, " env_tras: ", transf)
    return archive







