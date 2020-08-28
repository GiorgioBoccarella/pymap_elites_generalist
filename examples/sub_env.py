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
import sys

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import map_elites.wrap_model as mt_map_elites

from examples import generate_env
from map_elites import common_new as cm

#Seed MUST BE different from 0 (see gen_env)
#For each sim generate random seed
seed = 1
n = 6
#Generate one environment Pair
first_env = 1.3
envPair = generate_env.environmentPair(n, seed)
env = envPair(first_env)

#Fuction that creates a family of environment starting from env
envPair_c = generate_env.environmentPair(env, 0)

#Add all environments in list with respective distance
env_dist = [[0], [0.9], [1.3]]
envList = [] # This would be implemented as parameter
dist_env_add = np.array([0, 0.9])
print(dist_env_add)
for i in range(0, len(dist_env_add)):
    envList.append(envPair_c(dist_env_add[i]))

#The starting env is added here
envList.append(env)

envList = np.array(envList)
envList = envList.real


env_pair_l = {}
for i in range(0, len(envList)):
    env_pair_l.append(cm.Env(env_dist[i], envList[i]))

for i in range(0, len(env_pair_l)):
    print(env_pair_l[i].env_distance)
    print(env_pair_l[i].env)

#Generate all possible combination
# 10 bits = 1024 env
#With N = 4 => 16 sequences etc..
seq_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]
#From string to binary
for i in range(len(seq_list)):
    seq_list[i] = [int(numeric_string) for numeric_string in seq_list[i]]

archive = mt_map_elites.compute(max_evals=1e3, k=6, env_pair_list = env_pair_l, seq_list= seq_list,
            params=cm.default_params, log_file=None)

print(archive)
print("yee")