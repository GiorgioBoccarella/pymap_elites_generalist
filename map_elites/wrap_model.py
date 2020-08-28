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
import multiprocessing
from pathlib import Path
import sys
import random
from collections import defaultdict
from scipy.spatial import distance

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
    """This takes one of the K_i and returns all the mutation of K_i to the original genome"""
    np.random.seed()
    ran1 = np.random.randint(0, len(genome))
    #Here it's decided which K is mutated
    s_genome = genome[ran1]
    all_mut = np.tile(s_genome, (len(s_genome) - 1, 1))

    print(s_genome)

    for t in range(0, len(s_genome) - 1):
        all_mut[t][t] ^= 1

    g_mut = np.repeat(genome[np.newaxis, ...], len(all_mut), axis=0)

    for i in range(0, len(all_mut)):
        g_mut[i][ran1] = all_mut[i]

    print()
    return g_mut


def gen_lucky_mut(s_genome, all_g, env):

    # Calculate Fitness of the starting genome
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env)
    fit_s = -math.sqrt(resnorm)

    print("Fitness of wild genome:")
    print(fit_s)
    print()

    # Calculate fitness of all mutant genomes
    # fit_vec = np.empty([len(all_g)], dtype = float)

    l = []
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env)
        fit_m = -math.sqrt(resnorm)
        fit_diff = fit_m - fit_s
        if fit_diff > 0:
            l.append([fit_diff, i])

    # Fitness increase in env becomes probability that sum up to one
    p = np.array([fitness[0] for fitness in l])
    p /= p.sum()

    # Genomes trait value
    genomes = np.array([genome[0] for genome in l])

    # Here lucky mutation is selected
    draw = choice(genomes, 1, p=p)

    return draw

def make_ind(t):
    g, trj, f, pos = t
    return cm.Ind(g, trj, f, pos)



def add_to_archive(ind, archive):
    env_pair = cm.make_hashable(ind.position)
    archive[env_pair] = ind
    return 1

def env_pair_fitness(genome, env_pair):
        genome = genome.T

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[0])
        fit_m_1 = -math.sqrt(resnorm)

        [x, resnorm, residual] = lsq_f.lsqnonneg(genome, env_pair[1])
        fit_m_2 = -math.sqrt(resnorm)

        average_fitness = (fit_m_1 + fit_m_2)/2

        return average_fitness


def compute(max_evals=1e3,
            k=2,
            env_pair_list=[],
            seq_list=[],
            params=cm.default_params,
            log_file=None):
     
    print(params)

    assert (len(seq_list) >= k)

    # This I have to check
    n_env_pair = len(env_pair_list)

    # init archive (empty)
    archive = {}

    init_count = 0

    # main loop
    n_evals = 0 # number of evaluations
    successes = defaultdict(list) # count the successes
    while (n_evals < 2):
        #If environment is empty fill with random individuals
        if len(archive) < n_env_pair:
            for i in range(0, len(env_pair_list)):
                g = generate_genome(seq_list, k)
                env_c = env_pair_list[i]
                # Trajectory is initialized with first environment
                # Later new positions are appended
                trj = env_c.env_distance
                # Same for position but this is the actual position
                pos = env_c.env_distance
                # Generate fitness
                fit = env_pair_fitness(g, env_c.env)
                # Pack traits and feature, make Ind and add to archive
                *to_gen_ind, = g, trj, fit, pos
                ind = make_ind(to_gen_ind)
                print(ind.position)
                add_to_archive(ind, archive)
        else:
            for i in archive.keys():
                start_g = archive[i].genome
                all_mut = generate_all_mutants(start_g)
                env = env_pair_list
                gen_lucky_mut(start_g, all_mut, )

        else:
        for i in archive:
            generate_all_mutants(i.)

        for i in archive.keys():
            print(archive[i].genome)
        n_evals += 1

    return archive
