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
import multiprocessing
from pathlib import Path
import sys
import random
from collections import defaultdict
from sklearn.neighbors import KDTree
from scipy.spatial import distance
import matplotlib.pyplot as plt
import collections

def makehash():
    return collections.defaultdict(makehash)

from map_elites import common_new as cm


def distance(env1, env2):
    d = np.linalg.norm(env1 - env2)
    return d


# evaluate a single vector (z) with a function f and return a species
# t = vector, function
def make_ind(t):
    g, x, y, f, position = t
    return cm.Ind(g, x, y, f, position)


def mutate(ind):
    z = ind.copy()
    # select a random trait
    for i in range(0, 100):
        T_x = random.randint(0, len(ind) - 1)
    # select a random trait different from previous one
        T_y = random.choice([i for i in range(0, len(ind)) if i != T_x])
        step = min(ind[T_x], (-ind[T_y] + 1))
        if step > 0.1:
            step = 0.1
            break
        else:
            continue

    z[T_x] -= step
    z[T_y] += step
    return z


def weighted_random_choice(w_env):
    sum_f = 0
    for idx in w_env.values():
        sum_f += idx.fitness
    pick = random.uniform(0, sum_f)
    current = 0
    for key, value in w_env.items():
        pick -= value.fitness
        if pick <= 0:
            return key


def compute(dim_map=-1,
            dim_x=-1,
            f=None,
            end_sim=1e3,
            tasks=[],
            params=cm.default_params,
            sim=1,
            N=100,
            log_file=None):
    """Multi-task MAP-Elites
    - if there is no centroid : random assignation of niches
    - if there is no task: use the centroids as tasks
    - if there is a centroid list: use the centroids to compute distances
    when using the distance, use the bandit to select the tournament size (cf paper):

    Format of the logfile: evals archive_size max mean 5%_percentile, 95%_percentile

    Reference:
    Mouret and Maguire (2020). Quality Diversity for Multitask Optimization
    Proceedings of ACM GECCO.
    """
    #print(params)
    assert(f != None)
    assert(dim_x != -1)

    env_list = tasks
    env_list_d = makehash()

    N = 500

    # main loop
    gen = 0  # number of generations
    initialize = 0

    while (gen < end_sim):
        # Initialize the 4 environments with 10 individuals
        if initialize < 1:
            for i in range(len(env_list)):
                #for j in range(params['random_init'] * N):
                for j in range(N):
                    g = np.random.rand(dim_x)
                    g = (g / sum(g)) * (dim_x / 2)
                    x, y, fit = 0, 0, 0
                    position = i
                    ind_feature = [g, x, y, fit, position]
                    id = make_ind(ind_feature)
                    env_id = "env_{0}".format(i)
                    env_list_d[str(env_id)][j] = id
                    initialize = 1
        else:
        # If some individual exist then:
        # Calculate fitness for every id in every environment
            for i in env_list_d:
                for j in env_list_d[i]:
                    ind = env_list_d[i][j]
                    ind.fitness = f(ind.genome, env_list[int(ind.position)])
            # Selection and reproduction loop
            tmp = env_list_d.copy()
            for i in env_list_d:
                for j in range(len(env_list_d[i])):
                    tmp[i][j] = env_list_d[i][weighted_random_choice(env_list_d[i])]
            env_list_d = tmp
            # Variation loop
            # 10% of population mutate
            p_m = 0.1
            for i in env_list_d:
                for j in env_list_d[i].values():
                    mu = np.random.binomial(1, p_m)
                    if mu > 0:
                        j.genome = mutate(j.genome)
                    else:
                        continue

            filename = '/home/giorg/Documents/results/sim_' + str(sim) + '.dat'
            for i in env_list_d:
                average_f = 0
                for idx in env_list_d[i].values():
                    average_f += idx.fitness
                with open(filename, 'a+') as fl:
                    fl.write(str(i))
                    fl.write(" ,")
                    fl.write(str(average_f))
                    fl.write(" ,")
                    fl.write(str(gen))
                    fl.write("\n")

        gen += 1
