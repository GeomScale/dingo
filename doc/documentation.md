# Why to analyze metabolic networks?

Systems biology is an approach in biological and biomedical research aiming at deciphering theunderlying mechanisms and understand the full picture of the phenomena under study.  By dictatingall the biological levels of organization of living entities, the study of metabolism is the Holy Grail for biologists.  Metabolic network reconstruction has allowed for an in-depth insight into the molecularmechanisms by providing models that correlate the genome with molecular physiology.  The analysis of such  reconstructions  allow  the  identification  of  key  features  of  metabolism,  fundamental  for a great range of fields;  from the study of ecosystems resilience to this of complex diseases (e.g.,neuro-degenerative diseases) and advanced precision medicine.  

dingo is a python package for metabolic networks sampling and
analysis. To perform high dimensional sampling, dingo relies on the C++ package [volesti](https://github.com/GeomScale/volume_approximation), which provides several Markov Chain Monte Carlo (MCMC) algorithms for sampling high dimensional convex polytopes. dingo is part of [GeomScale](https://geomscale.github.io/) project.  

# Metabolic netwrks and dingo package

Systems Biology expands in all the different levels of living entities, from the
molecular, to the organismal and ecological level. The notion that
penetrates all  levels horizontally is *metabolism*; the
process that modifies molecules and  maintains the living state of a
cell or an organism through a set of chemical reactions. The reactions begin with a particular molecule
which they convert into some other molecule(s), while they are catalyzed by
enzymes in a key-lock relationship.
The quantitative relationships between the components of a reaction  is called *stoichiometry*.
Linked reactions, where the product of the first acts as the substrate for the
next, build up metabolic pathways. Each pathway is responsible for a certain
function. We can link together the aggregation of all the pathways that take
place in an organism (and their corresponding reactions)
and represent them mathematically using  the reactions' stoichiometry.
Therefore, at the species level, metabolism is a network of its metabolic pathways and we call
these representations *metabolic networks*.

Stoichiometric coefficients are the number of molecules a biochemical reaction
consumes and produces. The coefficients of all the reactions in a network,
with $m$ metabolites and $n$ reactions ($m \le n$), form
the stoichiometric matrix $S\in \RR^{m\times n$. 
The nullspace of $S$ corresponds to the steady states of the network:
<img src="https://render.githubusercontent.com/render/math?math=S \cdot v=0"> ,
where <img src="https://render.githubusercontent.com/render/math?math=v_{lb}\leq v\leq v_{ub}"> is the flux vector that contains  the fluxes
of each chemical reaction of the network while $v_{lb}$ and $v_{ub}$ denote lower and upper bound for each reaction flux respectively.

dingo performs the following methods (operations) on a given metabolic network:  

1. Samples from the flux space of a metabolic network using the  [Multiphase Monte Carlo Sampling algorithm](https://arxiv.org/abs/2012.05503) according to:  
- the uniform distribution,
- the multivariate exponential distribution,
- the multivariate Gaussian distribution.

2. Applies the FVA method.

3. Applies the [FBA method](https://www.nature.com/articles/nbt.1614).
 
# How to use dingo

## Run dingo from terminal

## Use dingo as a library