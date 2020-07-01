import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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


def subset_env(env_list):

    half = env_list[0:int(len(env_list)/2)]
    s_half = env_list[int(len(env_list)/2):int(len(env_list))]

    sum_env_1 = np.concatenate((half, s_half), axis=1)
    sum_env_2 = np.concatenate((half, s_half), axis=1)

    f = np.repeat(1, len(sum_env_1))
    s = np.repeat(0, len(sum_env_1))

    sum_env_1 = np.column_stack((sum_env_1, f))
    sum_env_2 = np.column_stack((sum_env_2, s))

    env_list_n = np.stack((sum_env_1, sum_env_2))

    return env_list_n


env_list = subset_env(env_list)

print(f_env)
print(s_env)








