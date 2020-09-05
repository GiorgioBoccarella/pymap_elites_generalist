import numpy as np
import collections
import random


def makehash():
    return collections.defaultdict(makehash)

class Ind:
     def __init__(self, genome, trait_x, trait_y, fitness, position=None):
        self.genome = genome
        self.trait_x = trait_x
        self.trait_y = trait_y
        self.fitness = fitness
        self.position = position



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


def add_to_env(ind_feature, envs_list_d):
    current_env = list(envs_list_d.items())[0][0]
    n_ind = int(len(list(envs_list_d.items())[0][1]))
    envs_list_d[current_env][n_ind + 1] = {n_ind, ind_feature}


def fit(ind, env):
    #np.linalg.norm(ind - env) is the mismatch
    # The 10 is arbitrary, worst fitness possible is 6.83772..
    if env[-1] == 0:
        s_side = int(len(ind)/2)
        f_e = env[0:s_side]
        f_i = ind[0:s_side]
        f = 2.01 - np.linalg.norm(f_i - f_e)
        f = np.exp(f)
    elif env[-1] == 1:
        s_side = int(len(ind) / 2)
        f_e = env[int(s_side):int(len(env) - 1)]
        f_i = ind[s_side:]
        f = 2.01 - np.linalg.norm(f_i - f_e)
        f = np.exp(f)
    return f


# Generate environments
# 10 bits = 1024 env


n = 4
env_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]

# From string to binary
for i in range(len(env_list)):
    env_list[i] = [int(numeric_string) for numeric_string in env_list[i]]

# Every environment sum is == constant n/2
env_list = [i for i in env_list if sum(i) == n/2]

# If we split the sequence in half we have the same sum
env_list = [i for i in env_list if half_sum(i) == 0]

env_list = subset_env(env_list)


env_list_d = makehash()


dim_x = 8

def make_ind(t):
    g, x, y, f, position = t
    return Ind(g, x, y, f, position)

for i in range(len(env_list)):
    # for j in range(params['random_init'] * N):
    for j in range(10):
        g = np.random.rand(dim_x)
        g = (g / sum(g)) * (dim_x / 2)
        x = 0
        y = 0
        f = 0
        position = i
        ind_feature = [g, x, y, f, position]
        id = make_ind(ind_feature)
        env_id = "env_{0}".format(i)
        env_list_d[str(env_id)][j] = id
        initialize = 1
    for i in env_list_d:
        for j in env_list_d[i]:
            ind = env_list_d[i][j]
            ind.fitness = fit(ind.genome, env_list[int(ind.position)])



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



g = np.empty([1, 10], dtype=bool)


def mutate_g(genome):
    ran1 = np.random.randint(0, len(genome))
    s_genome = genome[ran1]
    all_mut = np.tile(s_genome, (len(s_genome),1))

    for t in range(0, len(s_genome) - 1):
        all_mut[t][t] ^= 1

    return all_mut


def gen_lucky_mut(s_genome, all_g, env):
    """This generates a lucky mutation that is beneficial in at least one environment"""

    # Which env?
    env = env[np.random.binomial(1, 0.5, 1)]
    env = env.flatten()
    # print(s_genome)

    # Calculate Fitness of the starting genome
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env)
    fit_s = -math.sqrt(resnorm)

    print("Fitness before mutation:")
    print(fit_s)
    print()

    # Calculate fitness of all mutant genomes
    # fit_vec = np.empty([len(all_g)], dtype = float)

    l = []
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env)
        fit_m = -math.sqrt(resnorm)
        fit_diff = fit_m - fit_s
        if fit_diff > 0:
            l.append([fit_diff, i])

    if l != []:
        # Fitness increase in env becomes probability that sum up to one
        p = np.array([fitness[0] for fitness in l])
        p /= p.sum()

        # Genomes trait value
        genomes = np.array([genome[1] for genome in l])

        # Here lucky mutation is selected
        draw_n = choice(len(genomes), 1, p=p)
        mut_genome = genomes[draw_n]
        mut_genome = mut_genome.reshape(s_genome.shape)
    else:
        mut_genome = s_genome
        print("Runned out of beneficial mutations")

    return mut_genome

# Calculate Fitness of the starting genome in A
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env[0])
    fit_A = -math.sqrt(resnorm)

    print("Fit A:")
    print(fit_A)
    print()

    # Calculate Fitness of the starting genome in B
    [x, resnorm, residual] = lsq_f.lsqnonneg(s_genome.T, env[1])
    fit_B = -math.sqrt(resnorm)

    print("Fit B:")
    print(fit_B)
    print()

    # Calculate fitness of all mutant genomes
    # fit_vec = np.empty([len(all_g)], dtype = float)

    l_a = []
    l_b = []
    for i in all_g:
        t = i.T
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env[0])
        fit_m_A = -math.sqrt(resnorm)
        [x, resnorm, residual] = lsq_f.lsqnonneg(t, env[1])
        fit_m_B = -math.sqrt(resnorm)
        fit_diff_A = fit_m_A - fit_A
        fit_diff_B = fit_m_B - fit_B
        if fit_diff_A > 0 or fit_diff_B > 0:
            l_a.append([fit_diff_A])
            l_b.append([fit_diff_B])
            print()

    print(i_d)

    if l_a == []:
        exit()

    print("L_a")
    print(l_a)
    print()

    print("L_b")
    print(l_b)
    print()