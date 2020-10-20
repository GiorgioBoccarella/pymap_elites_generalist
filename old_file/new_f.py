import math
import numpy as np
import random
from numpy.random import choice
import sys

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sim_script import lsq_f
from sim_script import generate_env
from map_elites import common as cm

#Seed MUST BE different from 0 (see gen_env)
#For each sim generate random seed
seed = 1
n = 6
#Generate one environment Pair
first_env = 1.3
envPair = generate_env.environmentPair(n, seed)
env = envPair(first_env)

#Fuction that creates a family of environemt starting from env
envPair_c = generate_env.environmentPair(env, 0)

#Add all environements in list with respective distance
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


l = []
for i in range(0, len(envList)):
    l.append(cm.Env(env_dist[i], envList[i]))

for i in range(0, len(l)):
    print(l[i].env_distance)
    print(l[i].env)



#Generate all possible combination
# 10 bits = 1024 env
#With N = 4 => 16 sequences and so o
seq_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]
#From string to binary
for i in range(len(seq_list)):
    seq_list[i] = [int(numeric_string) for numeric_string in seq_list[i]]


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

g = generate_genome(seq_list, 2)

print("This is G:")
g = g*1
print(g)
print()

print("This is E:")
print(env)
print()

[x, resnorm, residual] = lsq_f.lsqnonneg(g.T, env)

print("This is the fitness: ")
print(-math.sqrt(resnorm))

print()

genome = g.copy()

def mutate_g(genome):
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

m = mutate_g(genome)

print(m)

print()

#TODO print fitness for example sequence

def gen_lucky_mut(s_genome, all_g, env):

    #Calculate Fitness of the starting genome
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env)
    fit_s = -math.sqrt(resnorm)

    print("startfit")
    print(fit_s)

    #Calculate fitness of all mutant genomes
    #fit_vec = np.empty([len(all_g)], dtype=float)
    l = []
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env)
        fit_m = -math.sqrt(resnorm)
        fit_diff = fit_m - fit_s
        if fit_diff > 0:
            l.append([fit_diff, i])


    p = np.array([a_tuple[0] for a_tuple in l])
    p /= p.sum()
    genomes = np.array([a_tuple[0] for a_tuple in l])


    draw = choice(genomes, 1,
                  p=p)

    return draw


l_c = gen_lucky_mut(g, m, env)

print(l_c)



