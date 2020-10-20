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
import numpy as np
import os
# Default params are general for all simulation

initial_params = \
    {
        "seed": 7423,
        "l_n": 100,
        "env_list": [0.3, 0.9, 1.1, 1.3],
        "max_evals": 10,
        "sim": 10
    }


default_params_1 = \
    {
        "seed": 4326,
        'del': 0,
        'transfer': 0,
        "env_transfer": [1.1],
        "p": 0,
        'l': 100,
        'k': 2,
        'invasion_rate': 0,
        's_invasion': 1.1,
        'length_transfer': 100
    }



class Ind:
    def __init__(self, position, genome, modularity, trade_off, fitness, invasion_potential, sum_p0, sum_p1, ps):
        self.position = position
        self.genome = genome
        self.modularity = modularity
        self.trade_off = trade_off
        self.fitness = fitness
        self.invasion_potential = invasion_potential
        self.sum_p0 = sum_p0
        self.sum_p1 = sum_p1
        self.ps = ps



def make_hashable(array):
    return tuple(map(float, array))

def parallel_eval(evaluate_function, to_evaluate, pool, params):
    if params['parallel'] == True:
        s_list = pool.map(evaluate_function, to_evaluate)
    else:
        s_list = map(evaluate_function, to_evaluate)
    return list(s_list)

# define the name of the directory to be created
# Usually make the folder in advance so they can concatenate
folder = "/home/giorg/Documents/trial/#1/"
os.makedirs(folder)

print(folder)

# format: fitness, centroid, desc, genome \n
def __save_archive(archive, gen, sim, transfer_in, transf_n, inv_rate, inv, p):
    filename = str(folder) + 'archive_sim_' + '.dat'
    def write_array(a, f):
        for i in a:
            f.write(str(i) + ' ')
    with open(filename, 'a+') as f:
        for k in archive.values():
            f.write(str(k.position) + ' ')
            f.write(str(k.fitness) + ' ')
            f.write(str(k.modularity) + ' ')
            f.write(str(k.trade_off) + ' ')
            for kk in k.invasion_potential.values():
                for l in range(0, len(kk)):
                    f.write(str(kk[l]) + ' ')
            f.write(str(k.sum_p0) + ' ')
            f.write(str(k.sum_p1) + ' ')
            f.write(str(k.ps) + ' ')
            f.write(str(gen) + ' ')
            f.write(str(sim) + ' ')
            f.write(str(transfer_in) + " ")
            f.write(str(transf_n) + " ")
            f.write(str(inv_rate) + " ")
            f.write(str(inv) + " ")
            f.write(str(p) + " ")
            f.write("\n")




def __save_file_mig(invader, wild, epoch, sim, inv_rate, ips, wps, sum_p0_invader, sum_p1_invader, sum_p0_wild, sum_p1_wild):
    filename = str(folder) + 'invasion_id' + '.dat'
    with open(filename, 'a+') as f:
        f.write(str(invader) + " ")
        f.write(str(wild) + " ")
        f.write(str(epoch) + " ")
        f.write(str(sim) + " ")
        f.write(str(inv_rate) + " ")
        f.write(str(ips) + " ")
        f.write(str(wps) + " ")
        f.write(str(sum_p0_invader) + " ")
        f.write(str(sum_p1_invader) + " ")
        f.write(str(sum_p0_wild) + " ")
        f.write(str(sum_p1_wild) + " ")
        f.write("\n")


def save_params(params):
    filename = str(folder) + 'params' + '.dat'
    with open(filename, 'a+') as f:
        f.write(str(params) + " ")
        f.write('\n')


def save_env(env, sim):
    filename = str(folder) + 'env_list' + '.dat'
    def write_array(a, f):
        for i in a:
            f.write(str(i) + ' ')
    with open(filename, 'a+') as f:
        for e in env.values():
            e_1 = np.array(e[0])
            write_array(e_1, f)
            e_2 = np.array(e[1])
            write_array(e_2, f)
        f.write(str(sim) + " ")
        f.write('\n')

