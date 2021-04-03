# Why to analyze metabolic networks?

Systems biology is an approach in biological and biomedical research aiming at deciphering the underlying mechanisms and understand the full picture of the phenomena under study.  By dictating all the biological levels of organization of living entities, the study of metabolism is the Holy Grail for biologists.  Metabolic network reconstruction has allowed for an in-depth insight into the molecular mechanisms by providing models that correlate the genome with molecular physiology.  The analysis of such  reconstructions  allow  the  identification  of  key  features  of  metabolism,  fundamental  for a great range of fields;  from the study of ecosystems resilience to this of complex diseases (e.g.,neuro-degenerative diseases) and advanced precision medicine.  

dingo is a python package for metabolic networks sampling and
analysis. To perform high dimensional sampling, dingo relies on the C++ package [volesti](https://github.com/GeomScale/volume_approximation), which provides several Markov Chain Monte Carlo (MCMC) algorithms for sampling high dimensional convex polytopes. dingo is part of [GeomScale](https://geomscale.github.io/) project.  

# Metabolic networks and dingo package

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
the stoichiometric matrix <img src="https://render.githubusercontent.com/render/math?math=S\in \mathbb{R}^{m\times n}">. 
The nullspace of <img src="https://render.githubusercontent.com/render/math?math=S"> corresponds to the steady states of the network:
<img src="https://render.githubusercontent.com/render/math?math=S \cdot v=0"> ,
where <img src="https://render.githubusercontent.com/render/math?math=v_{lb}\leq v\leq v_{ub}"> is the flux vector that contains  the fluxes
of each chemical reaction of the network while <img src="https://render.githubusercontent.com/render/math?math=v_{lb}">
and <img src="https://render.githubusercontent.com/render/math?math=v_{ub}">
denote lower and upper bound for each reaction flux respectively.

dingo performs the following methods (operations) on a given metabolic network:  

1. Samples from the flux space of a metabolic network using the  [Multiphase Monte Carlo Sampling algorithm](https://arxiv.org/abs/2012.05503) according to:  
- the uniform distribution,
- the multivariate exponential distribution,
- the multivariate Gaussian distribution.

2. Applies the FVA method.

3. Applies the [FBA method](https://www.nature.com/articles/nbt.1614).

# How to use dingo

There are two ways to run dingo. The first way is to run dingo from terminal using the main function and the second is to use dingo as a library by importing its routines in your code.

## Run dingo from terminal

You can call dingo by providing the path to the `json` file (as those in [BiGG](http://bigg.ucsd.edu/models) dataset) of a model:

```
python -m dingo -i model.json
```
Otherwise, you can use the `matlab` script in `./ext_data` folder to transform a `.mat` file (also as those in BiGG dataset) into a `.mat` file that dingo can read,

```
python -m dingo -i dingo_model.mat
```

By default, dingo generates uniformly distributed steady states of the given metabolic network. In particular, it computes the full dimensional polytope implied by <img src="https://render.githubusercontent.com/render/math?math=S \cdot v=0,\ v_{lb}\leq v\leq v_{ub}"> and then samples from it using the [Multiphase Monte Carlo Sampling algorithm](https://arxiv.org/abs/2012.05503) (MMCS) with the target distribution being the uniform distribution. Finally, dingo saves to the current path --using `pickle` package--,  

(a) a file with the generated steady states and two vectors that contain the minimum and maximum values of each flux,
(b) a file that contains the full dimensional polytope and the matrices of the linear transformations that map the polytope to the initial space.

You could also specify the output directory,

```
python -m dingo -i model.json -o output_directory
```

Since for high dimensional networks the complete pipeline of dingo is time consuming, you could split it into separate steps.
You can ask dingo to complete only the preprocessing of a model; that is computing only the full dimensional polytope,

```
python -m dingo -i model.json -o output_directory -preprocess True
```

Then, you can use the output polytope to sample from it,

```
python -m dingo -poly output_polytope
```

However, the output polytope after a complete run and the termination of MMCS algorithm is much more rounded than the polytope just after the preprocessing. Thus, the sampling from that polytope is more efficient. You should use that polytope to sample additional steady states.

### Statistical guarantees

dingo provides two statistical guarantees for the quality of the generated sample based on two MCMC diagnostics, (a) the *effective sample size* (ESS)  and (b) the *potential scale reduction factor* (PSRF).  

dingo sets by default the target effective sample size (ESS) for each flux marginal equal to 1000. You could ask for a different value of ESS,

```
python -m dingo -i model.json -n 2000
```

You can also ask for an additional statistical guarantee by setting an upper bound on the values of the PSRF of each flux marginal,

```
python -m dingo -i model.json -psrf 1.1
```

Then, dingo samples until it achieves the target values of ESS and PSRF for each flux marginal.  

### Fast and robust computations with gurobi library

The default Linear Program (LP) solver that dingo uses is the `linprog` from package `scipy`. However, this solver might fail to perform the required computations, or would be too slow. dingo also uses the LP solver of [`gurobi`](https://www.gurobi.com/downloads/?campaignid=2027425882&adgroupid=77414946611&creative=406173548636&keyword=gurobi&matchtype=e&gclid=CjwKCAjwpKCDBhBPEiwAFgBzj6dqIVtvMCfl2VQEDb3t4azHbvHN5vKoqcKNtAgHC9iBTlt3kqK2aBoC7CYQAvD_BwE) library which much faster and accurate than `linprog`. You could ask for `gurobi` by the following flag,

```
python -m dingo -i model.json -s gurobi
```

The default choice is `linprog` since `gurobi` needs a license to work.



### Handling the objective function  

Given the metabolic network, dingo --by default-- restricts the flux space to the set of flux vectors that have an objective value equal to the optimal value of the function. dingo allows for a more  relaxed option where you could ask to sample flux vectors that have an objective value equal to at least a percentage of the optimal value,

```unix
python -m dingo -i model -opt 90
```

which ask for flux vectors with a value equal to the 90% of the optimal value.  



Moreover, you can ask dingo to ignore the objective function,

```
python -m dingo -i model -unbiased True
```

In this case dingo does not impose any additional restriction to the flux space.



### Alternative distributions

dingo can sample from the flux space according to two more alternative distributions compare to the uniform distribution (i.e. the spherical Gaussian distribution and the exponential distribution). To ask for an alternative target distribution use the flag `-d`,

```
python -m dingo -i model -d exponential
```

or,

```
python -m dingo -i model -d gaussian
```

When the choice is the exponential distribution, dingo ignores the objective function and does not impose a certain objective value to the flux vectors. When the choice is the gaussian distribution, dingo sets the mode to be the center of the maximum inscribed ball of the full dimensional polytope.



### FVA method

To call FVA method run the following command,

```
python -m dingo -i model.json -fva True -opt opt_percentage
```

while the `-opt` flag is optional with the default percentage equal to `100%`.

### FBA method

To call FBA method run the following command,

```
python -m dingo -i model.json -fba True
```



## Use dingo as a package

dingo provides several `python` modules and functions to integrate into your project. The following script provides the full set of imports that are necessary to use all the cmputational options that dingo provides.

```python
# external imports
import numpy as np
import pickle

# dingo imports
from dingo.fva import slow_fva
from dingo.fba import slow_fba
from dingo.loading_models import read_json_file
from dingo.inner_ball import slow_inner_ball
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import (
    gmscale,
    apply_scaling,
    remove_almost_redundant_facets,
    map_samples_to_steady_states,
)

# import the fast implementations if gurobi is available
try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError:
    pass

# import the C++ class representing a polytope exposed by Cython
from dingo import HPolytope
```

In the following scripts we assume that we have imported all the above.

### Read a model and sample steady states

```python
input_file_json = "path_to/ext_data/e_coli_core.json"

e_coli_network = read_json_file(input_file_json)

lb = e_coli_network[0]
ub = e_coli_network[1]
S = e_coli_network[2]
biomass_function = e_coli_network[6]

A, b, Aeq, beq, min_fluxes, max_fluxes = fast_fva(lb, ub, S, biomass_function)

N, N_shift = nullspace_sparse(Aeq, beq)
b = np.subtract(b, np.dot(A, N_shift))
A = np.dot(A, N)

A, b = remove_almost_redundant_facets(A, b)
res = gmscale(A, 0.99)
A, b = apply_scaling(A, b, res[0], res[1])
A, b = remove_almost_redundant_facets(A, b)

p = HPolytope(A,b)
p.fast_mmcs(1000, True)
```





