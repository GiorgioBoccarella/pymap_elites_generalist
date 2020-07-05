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



def fit(ind, env):
    #np.linalg.norm(ind - env) is the mismatch
    # The 10 is arbitrary, worst fitness possible is 6.83772..
    if env[len(env)-1] == 0:
        s_side = int(len(ind)/2)
        f_e = env[0:s_side]
        f_i = ind[0:s_side]
        f = 10 - np.linalg.norm(f_i - f_e)
    elif env[len(env)-1] == 1:
        s_side = int(len(ind) / 2)
        f_e = env[int(s_side):int(len(env) - 1)]
        f_i = ind[0:s_side]
        f = 10 - np.linalg.norm(f_i - f_e)
    return f

# evaluate a single vector (z) with a function f and return a species
# t = vector, function
def make_ind(t):
    g, x, y, f, position = t
    return cm.Ind(g, x, y, f, position)



# select the niche according to
def select_niche(x, z, f, centroids, tasks, t_size, params, use_distance=False):
    to_evaluate = []
    if not use_distance:
        # No distance: evaluate on a random niche
        niche = np.random.randint(len(tasks))
        s_niche = tasks[niche]
        d = distance(s_niche[0:int(len(s_niche) - 1)], z)
        to_evaluate += [(z, f, tasks[niche], centroids[niche, :], params)]
    else:
        # we select the parent (a single one), then we select the niche
        # with a tournament based on the task distance
        # the size of the tournament depends on the bandit algorithm
        niches_centroids = []
        niches_tasks = [] # TODO : use a kd-tree
        rand = np.random.randint(centroids.shape[0], size=t_size)
        for p in range(0, t_size):
            n = rand[p]
            niches_centroids += [centroids[n, :]]
            niches_tasks += [tasks[n]]
        cd = distance.cdist(niches_centroids, [x.centroid], 'euclidean')
        cd_min = np.argmin(cd)
        to_evaluate += [(z, f, niches_tasks[cd_min], niches_centroids[cd_min], params)]
    return to_evaluate, d


def mutate(ind):
    z = ind.copy()
    # select a random trait
    t_x = random.randint(0, len(ind) - 1)

    mu, sigma = 0, 0.05
    step = np.random.normal(mu, sigma, 1)
    z[t_x] += step
    return z

def weighted_random_choice(w_env):
    sum_f = 0
    for idx in w_env.values():
        sum_f += idx.fitness
    pick = random.uniform(0, sum_f)
    current = 0
    for key, value in w_env.items():
        current += value.fitness
        if current > pick:
            return key


def compute(dim_map=-1,
            dim_x=-1,
            f=None,
            end_sim=1e4,
            env_list=[],
            tasks = [],
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
    print(params)
    assert(f != None)
    assert(dim_x != -1)

    env_list = tasks
    env_list_d = makehash()


    # main loop
    gen = 0  # number of generations
    initialize = 0
    dim_x = 8

    while (gen < end_sim):
        # Initialize the 4 environments with 10 individuals
        if initialize < 1:
            for i in range(len(env_list)):
                #for j in range(params['random_init'] * N):
                for j in range(100):
                    g = np.random.rand(dim_x)
                    g = (g / sum(g)) * (dim_x / 2)
                    x, y, f = 0, 0, 0
                    position = i
                    ind_feature = [g, x, y, f, position]
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
                    ind.fitness = fit(ind.genome, env_list[int(ind.position)])
            # Selection and reproduction loop
            tmp = env_list_d.copy()
            for i in env_list_d:
                for j in range(len(env_list_d[i])):
                    tmp[i][j] = env_list_d[i][weighted_random_choice(env_list_d[i])]
            env_list_d = tmp
            # Variation loop
            # 10% of population mutate
            p_m = 0.001
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
                with open(filename, 'a+') as f:
                    f.write(str(i))
                    f.write(" ,")
                    f.write(str(average_f))
                    f.write(" ,")
                    f.write(str(gen))
                    f.write("\n")

            gen += 1

