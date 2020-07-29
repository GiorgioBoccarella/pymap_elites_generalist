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

import map_elites.mod_model as mt_map_elites
import map_elites.common_new as cm_map_elites
import map_elites.l

#TODO Fintess function to fix on all plastic trait
def fit(ind, env):
    #np.linalg.norm(ind - env) is the mismatch
    # The 10 is arbitrary, worst fitness possible is 6.83772..
    f = 4 - np.linalg.norm(ind - env)
    f = np.exp(f)
    return f

# dim_map, dim_x, function
px = cm_map_elites.default_params.copy()


#Generate all possible combination
# 10 bits = 1024 env
n = 4
seq_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]
#From string to binary
for i in range(len(seq_list)):
    seq_list[i] = [int(numeric_string) for numeric_string in seq_list[i]]


#TODO generate environment pair with proper distance
env_dist = np.array([0, 0.5, 1.0, 1.5])


def generate_genome(sequences, k):
    assert(len(sequences) >= k)
    l = len(sequences[0])
    g = np.empty([k, l], dtype=bool)
    sel_id = np.random.choice(len(sequences), k, replace=False)
    j = 0
    for seq in sel_id:
        g[j] = sequences[seq]
        j += 1
    return g


def mutate_g(genome):
    ran1 = np.random.randint(0, len(genome))
    s_genome = genome[ran1]
    all_mut = np.tile(s_genome, (len(s_genome), 1))

    for t in range(0, len(s_genome) - 1):
        all_mut[t][t] ^= 1

    print(all_mut)


n_sim = 1
#Env params is the centroid
for s in range(0, n_sim):
    archive = mt_map_elites.compute(dim_x=4, f=fit, env=env_list, k=2, env_params=[],
                                    seq_list=seq_list, end_sim=2e4, params=px, sim=n_sim)
    print(s)

