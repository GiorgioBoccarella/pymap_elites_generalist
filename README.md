# Python3 Map-Elites 
This repository contains "reference implementations" of:
- CVT-MAP-Elites (Vassiliades, Chatzilygeroudis, Mouret, 2017)
- Multitask-MAP-Elites (Mouret and Maguire, 2020)
- Tikhonov, M., Kachru, S., & Fisher, D. S. (2020).

CVT-MAP-Elites can be used instead of the standard MAP-Elites described in:
Mouret JB, Clune J. Illuminating search spaces by mapping elites. arXiv preprint arXiv:1504.04909. 2015 Apr 20.


## Dependencies

- python3
- numpy


## References:
If you use this code in a scientific paper, please cite:

**Main paper**: Mouret JB, Clune J. Illuminating search spaces by mapping elites. arXiv preprint arXiv:1504.04909. 2015 Apr 20.

**CVT Map-Elites**: Vassiliades V, Chatzilygeroudis K, Mouret JB. Using centroidal voronoi tessellations to scale up the multi-dimensional archive of phenotypic elites algorithm. IEEE Transactions on Evolutionary Computation. 2017 Aug 3.

**Variation operator**: Vassiliades V, Mouret JB. Discovering the Elite Hypervolume by Leveraging Interspecies Correlation. Proc. of GECCO. 2018.

**Multitask-MAP-Elites**: Mouret JB, Maguire G. Quality Diversity for Multi-task Optimization. Proc of GECCO. 2020.


## Basic usage
(you need to have the map_elites module map_elites in your Python path)

```python

mt_map_elites.compute_invasion_transfer_new(params_sim=cm.params_sim_t)

```
params_sim_t = \
    {
        "seed_e": 800, # seed environment
        "seed_s": 3450, # seed mutation
        "l_n": 100, # L lenght of environment
        "env_list": [0.1, 1.0, 1.3],# generate list of environment with Delta E
        "max_evals": 250, # max mutation steps
        "sim": 40, # number of simulation
        "del": [1.3], # Delete environment from list
        "p": 0, # Initial p regulator
        'l': 100, # lenght of genome k-basis
        'k': 4, # number of k
        'invasion': True, # Invasion true or false
        'invasion_rate': 1, # How often invasion (0-1)
        'avg': True, # Invade base on average fitness
        'env_transfer': [1.3], # Transfer in this env
        'length_transfer': 100, # Lenght Transfer
        "folder": "test_transfer_1.3_k4_l100/"
    }
```

See the `examples` directory for a few examples.
