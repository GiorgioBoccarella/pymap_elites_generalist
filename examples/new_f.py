import math
import numpy as np
import random
import sys

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from map_elites import lsq

#Generate all possible combination
# 10 bits = 1024 env
n = 4
seq_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]
#From string to binary
for i in range(len(seq_list)):
    seq_list[i] = [int(numeric_string) for numeric_string in seq_list[i]]


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

g = generate_genome(seq_list, 3)

print(g)

def mutate_g(genome):
    ran1 = np.random.randint(0, len(genome))
    s_genome = genome[ran1]
    all_mut = np.tile(s_genome, (len(s_genome) - 1, 1))

    print()
    print(s_genome)
    print()

    for t in range(0, len(s_genome) - 1):
        all_mut[t][t] ^= 1

    return all_mut


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