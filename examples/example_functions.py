import numpy as np
import random

#Generate environements
# 10 bits = 1024 env
n = 14
env_list = [bin(x)[2:].rjust(n, "0") for x in range(2**n)]


#From string to binary
for i in range(len(env_list)):
    env_list[i] = [int(numeric_string) for numeric_string in env_list[i]]

#Every environment sum is == 5
env_list = [i for i in env_list if sum(i) == 7]


#10 is also the number of traits that the individual can "invest" in
# T = (T1, T2, T3...T10)
# Sum of T_i is always 1 (trade-off)
# Mutation is for example T2 =- 0.1 and respectively T6 =+ 0.1

#Initialize one individual
ind1 = [np.random.dirichlet(np.ones(10),size=1)] #sum = 1
#ind2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  dtype=float)
print(np.sum(ind1)) #sum close to one but not exactly 1 #TODO

def mutate(ind):
    # select a random trait
    T_x = random.randint(0, 9)
    # select a random trait different from previous one
    T_y = random.choice([i for i in range(0, 10) if i != T_x])

    #Avoid ind[T_x] to go lower then 0 and ind[T_y] to go higher then 1
    step = np.random.uniform(0, min(ind[T_x], (-ind[T_y] + 1)))

    ind[T_x] -= step
    ind[T_y] += step
    return ind

def mutate(ind):
    # select a random trait
    T_x = random.randint(0, 9)
    # select a random trait different from previous one
    T_y = random.choice([i for i in range(0, 10) if i != T_x])

    #Avoid ind[T_x] to go lower then 0 and ind[T_y] to go higher then 1
    step = np.random.uniform(0, min(ind[T_x], (-ind[T_y] + 1)))

    ind[T_x] -= step
    ind[T_y] += step
    return ind


def mutateEnv(ind):
    z = ind.copy()
    # select a random trait
    for i in range(0, 100):
        T_x = random.randint(0, len(ind) - 2)
    # select a random trait different from previous one
        T_y = random.choice([i for i in range(0, 10) if i != T_x])
        step = min(ind[T_x], (-ind[T_y] + 1))
        if step >= 0.1:
            step = 1
            break
        else:
            continue

    z[T_x] -= step
    z[T_y] += step
    return z

fit = 0
#fitness(ind1, env_list[834])
def fitness(ind, env):
    #np.linalg.norm(ind - env) is the mismatch
    # The 7 is arbitrary, worst fitness possible is 6.83772..
    fit = 7 - np.linalg.norm(ind - env)
    return fit


#Print fitness of ind in every environment
#Computational costly but if we store fitness in every environement then we can see if they are "specializing"
#Fitness in a set of enviornment much greater then in a different subset
for i in range(len(env_list)):
    print(env_list[i], fitness(ind1, env_list[i]))




m = 5
n = 10
a = np.zeros((m, n), dtype=int)
cols = np.random.binomial(1, 0.7, size=n)
a[np.random.randint(0, m, size=cols.sum()), np.nonzero(cols)[0]] = 1

to_evaluate = []
to_evaluate_centroid = []
if len(archive) <= 2:
    # create a random individual perfectly specialized to one of the task
    x = np.random.randint(0, n_tasks)
    x = np.asarray(tasks[x])
    x = x.astype(float)
    type = 1
    # we take a random task
    n = np.random.randint(0, n_tasks)
    to_evaluate += [(x, f, tasks[n], type, centroids[n], params)]
    s_list = cm.parallel_eval(__evaluate, to_evaluate, pool, params)
    n_evals += len(to_evaluate)
    b_evals += len(to_evaluate)
    add_to_archive(s_list[0], archive)
    # create a random generalist
    y = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    y = y.astype(float)
    type = 0
    # we take a random task
    n = np.random.randint(0, n_tasks)
    to_evaluate += [(y, f, tasks[n], type, centroids[n], params)]
    s_list = cm.parallel_eval(__evaluate, to_evaluate, pool, params)
    n_evals += len(to_evaluate)
    b_evals += len(to_evaluate)
    add_to_archive(s_list[1], archive)





