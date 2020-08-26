import math
import numpy as np
import random
import sys

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from examples import lsq_f
from examples import generate_env


#TODO generate different environment for every simulation (different seed)
#In generate_env(x, y)
# x = L
# y = seed
envPair = generate_env.environmentPair(4, 3)
env = envPair(1.0).real


#Generate all possible combination
# 10 bits = 1024 env
#With N = 4 => 16 sequences and so on
n = 4
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
env = env[0]
print(env)
print()

[x, resnorm, residual] = lsq_f.lsqnonneg(g.T, env)

print("This is the fitness: ")
print(-math.sqrt(resnorm))

print()

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

    g_mut = np.repeat(g[np.newaxis, ...], len(all_mut), axis=0)

    for i in range(0, len(all_mut)):
        g_mut[i][ran1] = all_mut[i]

    print()


    return g_mut

m = mutate_g(g)

print()

#TODO print fitness for example sequence
#Problem is that you can not calculate the fitness if you do not have the full genome
def gen_lucky_mut(s_genome, all_g, env):

    #Calculate Fitness of the starting genome
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env)
    fit_s = -math.sqrt(resnorm)

    print(fit_s)

    #Calculate fitness of all mutant genome
    fit_vec = np.empty([len(all_g)], dtype=float)
    l = []
    j = 0
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env)
        fit_vec[i] = -math.sqrt(resnorm)
        l.append([fit_vec[j], i])
        j += 1

    l_c = []
    for j in range(0, len(l)):
        if l[j][0] > fit_s:
            l_c.append(l[j])

    return l_c , l

l_c, l = gen_lucky_mut(g, m, env)

print(l)
print(l_c)

exit()




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

def lucky_mut(all_mut, env):
    fit_vec = np.empty([len(all_mut)], dtype=float)
    for i in all_mut:
        fit_vec[i] = lsq.lsqnonneg(i.T, env)

    #fit as a probability that sums to one

    draw = choice(list_of_candidates, number_of_items_to_pick,
    p = probability_distribution)
    return mut





print()
mutate_g(g)
print()
print(g)