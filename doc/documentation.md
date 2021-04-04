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
with <img src="https://render.githubusercontent.com/render/math?math=m"> metabolites and <img src="https://render.githubusercontent.com/render/math?math=n"> reactions (<img src="https://render.githubusercontent.com/render/math?math=m \le n">), form
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

2. Applies the [FVA method](https://www.sciencedirect.com/science/article/abs/pii/S1096717603000582) .

3. Applies the [FBA method](https://www.nature.com/articles/nbt.1614).

# How to use dingo

There are two ways to run dingo. The first way is to run dingo from terminal using the main function and the second is to use dingo as a library by importing its routines in your code.

## Run dingo from terminal

You can read analytical instructions of how to use dingo from terminal by typing,

```
python -m dingo -h
```

In the sequel, we present the main options for a user.  

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

### Nullspace computation

During the preprocessing, dingo computes the matrix of the right nullspace of the augmented stoichiometric matrix. The default method is based on the QR decomposition of the matrix with  library [SuitSparse](https://people.engr.tamu.edu/davis/suitesparse.html). In particular, dingo uses the python wrapper [PySPQR](https://github.com/yig/PySPQR) of `SuiteSparse` which exploits the sparsity of the stoichiometric matrix. An alternative option is to ask for the function `nullspace()` in `scipy` library by using the flag `-null`,

```
python -m dingo -i model.json -null scipy
```

### Parallel computations

dingo provides also a parallel implementation of MMCS algorithm. This implementation leads to faster computations. To sample steady state with parallel MMCS,

```
python -m dingo -i model.json -pmmcs True --num_threads 2
```

The flag `--num_threads` requests the number of threads to be used and it can be replaced with flag `-nt` for simplicity.

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

### Load a model and perform computations

The following example samples uniformly distributed steady states from the `e_coli` model in the folder `./ext_data`. It presents the main pipeline that dingo uses to generate steady states.  

```python
input_file_json = "path_to/ext_data/e_coli_core.json"

# read the model
e_coli_network = read_json_file(input_file_json)
# for a .mat file use the function read_mat_file()

# extract the information about the model
lb = e_coli_network[0]
ub = e_coli_network[1]
S = e_coli_network[2]
metabolites = e_coli_network[3]
reactions = e_coli_network[4]
biomass_index = e_coli_network[5]
biomass_function = e_coli_network[6]

# call FVA method to restrict the flux space in the region of optimal solutions
# and to compute the minimum and maximum values of the fluxes in that region
A, b, Aeq, beq, min_fluxes, max_fluxes = slow_fva(lb, ub, S, biomass_function)

# compute the nullspace of the augmented stoichiometric matrix Aeq
# and apply the linear trasformation to the flux space to derive a
# full dimensional polytope
N, N_shift = nullspace_sparse(Aeq, beq)
b = np.subtract(b, np.dot(A, N_shift))
A = np.dot(A, N)

# apply two heuristics to improve numerical accuracy and speed
A, b = remove_almost_redundant_facets(A, b)
res = gmscale(A, 0.99)
A, b = apply_scaling(A, b, res[0], res[1])
A, b = remove_almost_redundant_facets(A, b)

# initialize the class of polytopes and call mmcs algorithm to sample 
# from the full dimensional polytope
P = HPolytope(A, b)
A_rounded, b_rounded, T, T_shift, samples = P.slow_mmcs(1000, True)

# map the samples back to the initial space to obtain the steady states
steady_states = map_samples_to_steady_states(samples, T, T_shift, N, N_shift)
```

To exploit fast computations with `gurobi` library you have to replace `slow` with `fast` in the above script, i.e. call `fast_fva()` and `p.fast_mmcs()`.  The function `nullspace_sparse()`uses the `suitesparse` library to compute the nullspace. To use the nullspace computation of `scipy` library call the function `nullspace_dense()`.  

The  functions `read_json_file()`and `read_mat_file()` return the lower and upper bounds for each reaction flux, the stoichiometric matrix, a list that contains the metabolites, a list that contains the reactions, the objective function of the biomass and the index of the biomass pseudo-reaction.  

The matrix `Aeq` that `fva` returns, is equal to the stoichiometric matrix augmented by some rows of matrix `A` that define redundant facets in the initial space. Then, we compute the matrix of the right nullspace of `Aeq` to restrict the initial polytope onto that space to obtain the full dimensional polytope <img src="https://render.githubusercontent.com/render/math?math=P = \{ x\in\mathbb{R}^n\ |\ Ax\leq b \}">.

The function `gmscale()` computes a scaling for the full dimensional polytope <img src="https://render.githubusercontent.com/render/math?math=P"> to improve the numerical accuracy and the function `remove_almost_redundant_facets()`removes the facets of <img src="https://render.githubusercontent.com/render/math?math=P"> that are redundant after the preprocessing. dingo uses those two functions to improve the runtime.  

The algorithm MMCS unifies sampling with rounding and thus, after termination it has computed a linear transformation that puts the initial full dimensional polytope into an approximate well-rounded position. Consequently, the member functions `slow_mmcs()` and `fast_mmcs()` of the polytope class, return the matrix <img src="https://render.githubusercontent.com/render/math?math=A\in\mathbb{R}^{m\times n}">and the vector <img src="https://render.githubusercontent.com/render/math?math=b\in\mathbb{R}^m"> that define the rounded polytope and the linear transformation, defined by the matrix `T` and the vector `T_shift`, that maps the samples from the rounded polytope to the initial full dimensional polytope <img src="https://render.githubusercontent.com/render/math?math=P">. 

### Sample additional steady states faster

To sample an additional set of steady states one should use the matrix `A_rounded` and the vector `b_rounded` to define a new polytope and sample from it,

```python
P = HPolytope(A, b)
A_rounded, b_rounded, T, T_shift, samples = P.slow_mmcs(1000, True)

P_rounded = HPolytope(A_rounded, b_rounded)
A_rounded_new, b_rounded_new, T_new, T_shift_new, samples2 = P_rounded.slow_mmcs(1000, True)

T = np.dot(T, T_new)
T_shift = np.add(T_shift, T_shift_new)
steady_states_new = map_samples_to_steady_states(samples2, T, T_shift, N, N_shift)
```

Similarly, let the file `output_polytope` that the following command saves in your working directory,

```
python -m dingo -i model.json
```

Then, to load the rounded polytope and sample from it you should use the following commands,

```python
file = open(args.polytope, "rb")
polytope_matrices = pickle.load(file)
file.close()

A_rounded = polytope_matrices[0]
b_rounded = polytope_matrices[1]
N = polytope_matrices[2]
N_shift = polytope_matrices[3]
T = polytope_matrices[4]
T_shift = polytope_matrices[5]

P_rounded = HPolytope(A_rounded, b_rounded)
A_rounded_new, b_rounded_new, T_new, T_shift_new, samples = P_rounded.slow_mmcs(1000, True)

T = np.dot(T, T_new)
T_shift = np.add(T_shift, T_shift_new)
steady_states_new = map_samples_to_steady_states(samples, T, T_shift, N, N_shift)
```

### Perform FBA and FVA methods

To use dingo for FVA or FBA you could use the following script,

```python
e_coli_network = read_json_file(input_file_json)

lb = e_coli_network[0]
ub = e_coli_network[1]
S = e_coli_network[2]
biomass_function = e_coli_network[6]

A, b, Aeq, beq, min_fluxes, max_fluxes = slow_fva(lb, ub, S, biomass_function)

optimum_solution, optimum_value = slow_fba(lb, ub, S, c)
```

Of course you could give as input any linear objective function.

