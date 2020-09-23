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

import math
import numpy as np
import multiprocessing
from pathlib import Path
import sys
import random
from collections import defaultdict
#from sklearn.cluster import KMeans

default_params = \
    {
        "seed": 85,
        "l_n": 4,
        "env_list": [0.2, 0.5, 0.9, 1.1, 1.3],
        "k": 2,
        "sim": 10,
        "max_evals": 450
    }

class Ind:
    def __init__(self, genome, trajectory, fitness, fit1, fit2, position=None):
        self.genome = genome
        self.trajectory = trajectory
        self.fitness = fitness
        self.fit1 = fit1
        self.fit2 = fit2
        self.position = position



def make_hashable(array):
    return tuple(map(float, array))

def parallel_eval(evaluate_function, to_evaluate, pool, params):
    if params['parallel'] == True:
        s_list = pool.map(evaluate_function, to_evaluate)
    else:
        s_list = map(evaluate_function, to_evaluate)
    return list(s_list)


# format: fitness, centroid, desc, genome \n
# fitness, centroid, desc and x are vectors
def __save_archive(archive, gen, sim):
    def write_array(a, f):
        for i in a:
            f.write(str(i) + ' ')
    filename = '/home/giorg/Documents/results/archive_sim_' + '.dat'
    with open(filename, 'a+') as f:
        for k in archive.values():
            f.write(str(k.fitness) + ' ')
            f.write(str(k.fit1) + ' ')
            f.write(str(k.fit2) + ' ')
            #write_array(k.genome, f)
            #f.write(str(k.trajectory))
            f.write(str(k.position) + ' ')
            f.write(str(gen) + ' ')
            f.write("\n")

def __save_file(tradeoff, env,  gen, sim):
    def write_array(a, f):
        for i in a:
            f.write(str(i) + ' ')
    filename = '/home/giorg/Documents/results/tradeoff_sim_' + '.dat'
    with open(filename, 'a+') as f:
        f.write(str(tradeoff) + ' ')
        f.write(str(env) + ' ')
        f.write(str(gen) + ' ')
        f.write(str(sim) + " ")
        f.write("\n")

def __save_file_mut(vec_mut, vec_off):
    def write_array(a, f):
        for i in a:
            f.write(str(i) + ' ')
    filename = '/home/giorg/Documents/results/mut_traj_sim_' + '.dat'
    with open(filename, 'a+') as f:
        for k in vec_mut:
            f.write(str(k) + " ")
            write_array(vec_mut[k], f)
            f.write("\n")
            f.write(str(k) + "t ")
            write_array(vec_off[k], f)
            f.write("\n")