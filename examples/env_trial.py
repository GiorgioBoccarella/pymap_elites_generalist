import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random

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
#Every environment sum is == 5
env_list = [i for i in env_list if sum(i) == n/2]

#if we split the sequence in half we have the same sum
env_list = [i for i in env_list if half_sum(i) == 0]






















env_list = np.array(env_list)

n_e = len(env_list)

map = np.zeros((n_e, n_e))
map1 = np.zeros((n_e, n_e))

for i in range(n_e):
    for j in range(n_e):
        map[i][j] = distance(env_list[i], env_list[j])

#order = map[0][:]
#env_list1 = env_list[order.argsort()]

#plt.imshow(map, cmap='hot', interpolation='nearest')
#plt.show()

for i in range(n_e):
    for j in range(n_e):
        map1[i][j] = distance(env_list1[i], env_list1[j])

mask = np.triu(map1)
ax = sns.heatmap(map1, mask=mask)
plt.show()

#TODO
# I could append the coordinates to a vector and then plot a line
# Implementing in main model

for key, value in dict.iteritems():
    temp = [key,value]
    dictlist.append(temp)


def plot_best(the_map, occ_n, iteration_number):
    mask = np.triu(map)
    ax = sns.heatmap(map, mask=mask)

    x = [0.5] + [x + 0.5 for x in occ_n[0:len(occ_n) - 1]] + [len(map) - 0.5]
    y = [0.5] + [x + 0.5 for x in occ_n[1:len(occ_n)]] + [len(map) - 0.5]

    plt.plot(x, y, marker='o', linewidth=4, markersize=12, linestyle="-", color='white')
    #plt.savefig('images/new1000plot_%i.png' % (iteration_number), dpi=300)
    plt.show()




import matplotlib.pyplot as plt


xx = range(5)
yy = range(5)
[plt.plot([x,x],[min(yy),max(yy)],color='k') for x in xx]
[plt.plot([min(xx),max(xx)],[y,y],color='k') for y in yy]

ind = np.array([0.5, 0.5])
time = 0


def initialize_map(N):

    the_map = np.zeros((N, N))

    for i in range(0, N):
        for j in range(0, i):
                the_map[j][i] = the_map[i][j]

    return the_map

def plot_best(the_map, id_vec):

    x = [0.5] + [x + 0.5 for x in id_vec[0:len(id_vec) - 1]] + [len(id_vec) - 0.5]
    y = [0.5] + [x + 0.5 for x in id_vec[1:len(id_vec)]] + [len(id_vec) - 0.5]

    plt.plot(x, y, marker='o', linewidth=4, markersize=12, linestyle="-", color='white')
    plt.show()


while time < 500:

    plt.plot()


	plt.pause(0.01)

	time += dt


def plot_best(id_vec):

    for i in id_vec:
        x = id_vec/2 - 0.5
        y = id_vec/2 - 0.5
    plt.scatter(x, y, color='red')
    plt.show()


import numpy as np
... import matplotlib.pyplot as plt
... r = np.arange(5)
... p = np.arange(5)
... R,P = np.meshgrid(r,p)
... data = np.random.random((5,5))
... xx = range(5)
... yy = range(5)
... [plt.plot([x,x],[min(yy),max(yy)],color='k') for x in xx]
... [plt.plot([min(xx),max(xx)],[y,y],color='k') for y in yy]
... plt.scatter(R[:-1,:-1]+0.5,P[:-1,:-1]+0.5, color = 'blue')