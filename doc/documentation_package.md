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