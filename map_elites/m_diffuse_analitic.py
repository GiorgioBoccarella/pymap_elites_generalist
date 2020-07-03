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
import seaborn as sns

from map_elites import common_source as cm


def distance(env1, env2):
    d = np.linalg.norm(env1 - env2)
    return d


def add_to_archive(ind, env_id, env_list):
    env_list[env_id] = ind

def fit(ind, env):
    #np.linalg.norm(ind - env) is the mismatch
    # The 10 is arbitrary, worst fitness possible is 6.83772..
    if env[len(env)-1] == 0:
        s_side = int(len(ind)/2)
        f_e = env[0:s_side]
        f_i = ind[0:s_side]
        f = len(ind) - np.linalg.norm(f_e - f_i)
    elif env[len(env)-1] == 1:
        s_side = int(len(ind) / 2)
        f_e = env[int(s_side):int(len(env) - 1)]
        f_i = ind[0:s_side]
        f = len(ind) - np.linalg.norm(f_e - f_i)
    return f

# evaluate a single vector (z) with a function f and return a species
# t = vector, function
def __evaluate(t):
    g, x, y, f, position = t
    # from position extract environment
    #fitnes must be += because multiple env
    fitness = fit(g, env)
    return cm.Ind(g, x, y, fitness, position)

# bandit opt for optimizing tournament size
# probability matching / Adaptive pursuit Thierens GECCO 2005
# UCB: schoenauer / Sebag
# TODO : params for values, and params for window


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
    for i in range(0, 100):
        t_x = random.randint(0, len(ind) - 1)

    mu, sigma = 0, 0.1
    step = np.random.normal(mu, sigma, 1)
    z[t_x] += step
    return z


def compute(dim_map=-1,
            dim_x=-1,
            f=None,
            end_sim=1e5,
            envs_list=[],
            params=cm.default_params,
            sim= 1,
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

    # init archive (empty). Archive here is dictionary of dictionary
    # Each dictionary is an environment and and in each environments there are several individual
    envs_list_d = {}

    for x in range(len(envs_list)):
        y = "env_{0}".format(x)
        envs_list_d[str(y)] = {}

    # init multiprocessing
    num_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(num_cores)

    # main loop
    gen = 0  # number of generations
    while (gen < end_sim):
        to_evaluate = []
        for e in envs_list_d:
        if len(e.) <= params['random_init'] * envs_list:
            # initialize the map with random individuals
            for i in range(0, params['random_init_batch']):
                # create a random individual perfectly specialized to one of the task
                x = np.random.randint(0, envs_list)
                x = np.repeat(0.5, int(len(envs_list[x]) - 1))
                x = x.astype(float)
                # we take a random task
                n = np.random.randint(0, envs_list)
                to_evaluate += [(x, f, envs_list[n], params)]
            e_list = cm.parallel_eval(__evaluate, to_evaluate, pool, params)
            for i in range(0, len(list(e_list))):
                add_to_archive(e_list[i], archive)
        else:
            # main variation/selection loop
            keys = list(archive.keys())
            # we do all the randint together because randint is slow
            rand1 = np.random.randint(len(keys), size=params['batch_size'])
            for n in range(0, params['batch_size']):
                # ind selection
                x = archive[keys[rand1[n]]]
                # add variation
                z = mutate(x.x)
                # different modes for multi-task (to select the niche)
                to_evaluate += select_niche(x, z, f, centroids, tasks, t_size, params, use_distance)[0]
                d = select_niche(x, z, f, centroids, tasks, t_size, params, use_distance)[1]
            # parallel evaluation of the fitness
            s_list = cm.parallel_eval(__evaluate, to_evaluate, pool, params)
            n_evals += len(to_evaluate)
            b_evals += len(to_evaluate)
            # natural selection
            suc = 0
            for i in range(0, len(list(s_list))):
                suc += add_to_archive(s_list[i], archive, d)[0]
                d += add_to_archive(s_list[i], archive, d)[1]
            if use_distance:
                successes[t_size] += [(suc, n_evals)]
        if use_distance: # call the bandit to optimize t_size
            t_size = bandit(successes, n_tasks)

        # write archive
        if params['dump_period'] != -1 and b_evals > params['dump_period']:
            cm.__save_archive(archive, n_evals, sim)
            b_evals = 0
            n_e = [len(v) for v in successes.values()]
            print(n_evals, n_e)
            env_l = list(archive.keys())
            env_a = np.array(env_l)
            keys_a = archive.keys()
            second_word = [archive[x] for x in keys_a]
            #for i in second_word:
            #   print(i.x)

            #vals = np.fromiter(archive.items().x, dtype=float)
            ###################################################### My mod
            env_l = list(archive.keys())
            env_a = np.array(env_l)
            nrows, ncols = 6, 6
            image = np.zeros(nrows * ncols)

            # Set every other cell to a random number (this would be your data)
            j = 0
            for i in second_word:
                image[int(env_a[j])] = np.std(i.x)
                j += 1

            for i in range(len(image)):
                if image[i] >= 0.25:
                    image[i] = 1
                elif image[i] < 0.25 and image[i] > 0:
                    image[i] = 2


            # Reshape things into a 9x9 grid.
            image = image.reshape((nrows, ncols))

            dct = {1: 0., 2: 5., 0: 20.}
            n = [[dct[i] for i in j] for j in image]
            plt.imshow(n, cmap='brg', vmin=1, vmax=10)
            plt.savefig('/home/giorg/Documents/plots/map_plot_%i.png' %(n_evals), dpi=300)
            plt.close
            ######################################
            np.savetxt('t_size.dat', np.array(n_e))
        if log_file != None:
            fit_list = np.array([x.fitness for x in archive.values()])
            log_file.write("{} {} {} {} {} {} {}\n".format(n_evals, len(archive.keys()), fit_list.max(), np.mean(fit_list), np.median(fit_list), np.percentile(fit_list, 5), np.percentile(fit_list, 95)))
            log_file.flush()
    cm.__save_archive(archive, n_evals, sim)
    filename_dtot = '/home/giorg/Documents/results/archive_sim_dtot_s_' + '.dat'
    with open(filename_dtot, 'a+') as f:
        f.write(str(sim) + ' ')
        f.write("\n")
    return archive



# a small test
if __name__ == "__main__":
    def rastrigin(xx):
        x = xx * 10.0 - 5.0
        f = 10 * x.shape[0]
        for i in range(0, x.shape[0]):
            f += x[i] * x[i] - 10 * math.cos(2 * math.pi * x[i])
        return -f, np.array([xx[0], xx[1]])
    # CVT-based version
    my_map = compute(dim_map=2, dim_x = 10, n_niches=1500, f=rastrigin)

