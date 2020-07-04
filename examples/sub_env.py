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
import sys

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import map_elites.m_diffuse_analitic as mt_map_elites
import map_elites.common_source as cm_map_elites


def fitness(ind, env):
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

def distance(env1, env2):
    d = np.linalg.norm(env1 - env2)
    return d

def half_sum(sequence):
    sequence = np.array(sequence)

    n = len(sequence)
    x = sequence[:int((n/2))]
    y = sequence[int(-(n/2)):]
    diff = sum(x - y)

    return diff


def subset_env(env):

    half = env[0:int(len(env)/2)]
    s_half = env[int(len(env)/2):int(len(env))]

    sum_env_1 = np.concatenate((half, s_half), axis=1)
    sum_env_2 = np.concatenate((half, s_half), axis=1)

    f = np.repeat(1, len(sum_env_1))
    s = np.repeat(0, len(sum_env_1))

    sum_env_1 = np.column_stack((sum_env_1, f))
    sum_env_2 = np.column_stack((sum_env_2, s))

    env_list_n = np.vstack((sum_env_1, sum_env_2))

    return env_list_n

# dim_map, dim_x, function
px = cm_map_elites.default_params.copy()
px["dump_period"] = 1000
px["parallel"] = True
px["batch_size"] = 10


#Generate environements
# 10 bits = 1024 env
n = 4
env_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]

#From string to binary
for i in range(len(env_list)):
    env_list[i] = [int(numeric_string) for numeric_string in env_list[i]]
#Every environment sum is == constant n/2
env_list = [i for i in env_list if sum(i) == n/2]

#if we split the sequence in half we have the same sum
env_list = [i for i in env_list if half_sum(i) == 0]

env_list = subset_env(env_list)

dim_x = n

n_sim = 1

for s in range(0, n_sim):
    archive = mt_map_elites.compute(dim_x=dim_x, f=fitness, tasks=env_list, max_evals=1e5, params=px, sim=n_sim,
                                    log_file=open('mt_no_dist.dat', 'w'))
    print(s)



def add_to_env(ind_feature, envs_list_d):
     current_env = list(envs_list_d.items())[0][0]
     envs_list_d[current_env] = {int(len(list(envs_list_d.items())[0][1])), ind_feature}

for i in range(10):
    g = np.random.rand(10)
    g = (g / sum(g)) * ( 10 / 2)
    x = 0
    y = 0
    f = 0
    position = 1
    ind_feature = [g, x, y, f, position]
    add_to_env(ind_feature, envs_list_d)