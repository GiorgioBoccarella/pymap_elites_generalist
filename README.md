# Python3 Map-Elites and the evolution of trade-offs

## Informal abstract:

Specialist and Generalist often coexist in nature and the classic explanation for
this coexistence is the presence of trade-offs. Theory often assumes that a
‘jack-of-all-trades is a master of none’ and generalist always encounter genetic
trade-off that stops them from displacing specialist. The presence of these fixed
trade-offs has been hard to detect empirically and when present it may be
ameliorated over time. For these reasons in this work, I relaxed this assumption
and considered the genome architecture with features as trade-off and
modularity as emergent and evolving properties. The framework here used
does not aim to be specific for a particular species but aims to uncover deeper
fundamental question of evolution in a pair of similar or dissimilar
environments which is the scenario respectively encountered by specialist and
generalist. Briefly, individuals have a simple bacteria metabolic network and fitness depend on matching a certain environment stoichiometry. 
I hereby show that evolved genome property has a strong influence
on evolutionary trajectory that may help or hinder evolution of a generalist or
specialist lineage. In particular, I demonstrated that generalists suffer from
specialist invasions that cause the fixation of nonoptimal genome architecture
and ultimately hampers the capacity of a generalist to be well adapted to their
niche. Hence demonstrating that no cost generalism is possible but hard to
reach due to specialist competition. In this case, successful invasions can have
important consequences on genealogical structures and potentially diminish
the capacity for adaptation to multiple environments. One factor that is relevant
for this effect is the initial state of the genome at the beginning of evolution.
Previous evolutionary history is shown to be a key factor for predicting the
evolutionary trajectory of a lineage.


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


