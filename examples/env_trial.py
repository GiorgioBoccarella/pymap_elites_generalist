import math
import numpy as np


def distance(env1, env2):
    d = np.linalg.norm(env1 - env2)
    return d

#Generate environements
# 10 bits = 1024 env
n = 4
env_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]

#From string to binary
for i in range(len(env_list)):
    env_list[i] = [int(numeric_string) for numeric_string in env_list[i]]
#Every environment sum is == 5
env_list = [i for i in env_list if sum(i) == n/2]
env_list = np.array(env_list)