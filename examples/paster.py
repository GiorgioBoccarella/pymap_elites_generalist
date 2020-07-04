import numpy as np
import collections

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


envs_list_d = makehash()


def make_ind(t):
    g, x, y, f, position = t
    return Ind(g, x, y, f, position)


for i in range(4):
    print(i)
    for j in range(10):
        g = np.random.rand(10)
        g = (g / sum(g)) * (10 / 2)
        x = 0
        y = 0
        f = 0
        position = 1
        ind_feature = [g, x, y, f, position]
        id = make_ind(ind_feature)
        env_id = "env_{0}".format(i)
        envs_list_d[str(env_id)][str(j)] = id
