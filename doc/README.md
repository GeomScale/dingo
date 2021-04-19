# Metabolic networks and dingo package

Metabolic network reconstruction has allowed for an in-depth insight into the molecular mechanisms by providing models that correlate the genome with molecular physiology.  The analysis of such  reconstructions  allow  the  identification  of  key  features  of  metabolism,  fundamental  for a great range of fields;  from the study of ecosystems resilience to this of complex diseases (e.g., neuro-degenerative diseases) and advanced precision medicine.  

dingo is a python package for metabolic network analysis. dingo performs the following methods (operations) on a given metabolic network:  

1. Samples from the flux space of a metabolic network using the  [Multiphase Monte Carlo Sampling algorithm](https://arxiv.org/abs/2012.05503) according to:  

- the uniform distribution,
- the multivariate exponential distribution,
- the multivariate Gaussian distribution.

2. Applies the [FVA method](https://www.sciencedirect.com/science/article/abs/pii/S1096717603000582) .

3. Applies the [FBA method](https://www.nature.com/articles/nbt.1614).

To perform high dimensional sampling, dingo relies on the C++ package [volesti](https://github.com/GeomScale/volume_approximation), which provides several Markov Chain Monte Carlo (MCMC) algorithms for sampling high dimensional convex polytopes. dingo is part of [GeomScale](https://geomscale.github.io/) project.  

# How to use dingo

There are two ways to run dingo. The first way is to run dingo from terminal using the main function and the second is to use dingo as a library by importing its routines in your code.

You can find more details in the following links,

1. [Use dingo from terminal](https://github.com/GeomScale/dingo/blob/develop/doc/documentation_main.md)
2. [Use dingo as a python package](https://github.com/GeomScale/dingo/blob/develop/doc/documentation_package.md)