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
import random

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import map_elites.model_functions as mt_map_elites

from examples import generate_env
from map_elites import common_invasion as cm


params = cm.default_params

# Seed MUST BE different from 0 (see gen_env)
# For each sim generate random seed
seed = params["seed"]
l_n = params["l_n"]
env_list = params["env_list"]



def environment_from_params(env_list_v, l_n, seed):
    example_env = env_list_v[len(env_list_v) - 1]
    envPair = generate_env.environmentPair(l_n, seed)
    env = envPair(example_env)

    envPair_c = generate_env.environmentPair(env, 0)

    envList = []

    for i in range(0, len(env_list_v)):
        envList.append(envPair_c(env_list_v[i]))

    envList = np.array(envList)
    envList = envList.real

    return envList

all_env_sim = []

for i in range(0, 100 + 2):
    new_seed = seed + i
    envList = environment_from_params(env_list, l_n, new_seed)
    env_pair_d = {}
    for d, s in zip(env_list, envList):
        env_pair_d[d] = s
    all_env_sim.append(env_pair_d)


#Generate all possible combination
# 10 bits = 1024 env
#With N = 4 => 16 sequences etc..
seq_list = [bin(x)[2:].rjust(l_n, "0") for x in range(2**l_n)]
#From string to binary
for i in range(len(seq_list)):
    seq_list[i] = [int(numeric_string) for numeric_string in seq_list[i]]

#Is the n divisible by 4? 25% are 1
#assert(l_n % n_one == 0)
#seq_list = [i for i in seq_list if sum(i) == n_one]

print("Sequences Generated")
#archive = mt_map_elites.compute_invasion(max_evals=params["max_evals"], k=params["k"], env_pair_dict_l=all_env_sim, seq_list=seq_list, sim=params["sim"],
#           params=cm.default_params)

archive = mt_map_elites.compute_versatility(max_evals=params["max_evals"], k=params["k"], env_pair_dict_l=all_env_sim, seq_list=seq_list, sim=params["sim"],
            params=cm.default_params_1)


archive = mt_map_elites.compute_versatility(max_evals=params["max_evals"], k=params["k"], env_pair_dict_l=all_env_sim, seq_list=seq_list, sim=params["sim"],
            params=cm.default_params_2)

print("DONE")