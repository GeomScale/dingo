## Use dingo as a package

It quite simple to use dingo in your code.  In general, dingo provides two classes:  

- `metabolic_network` represents a metabolic network
- `polytope_sampler` can be used to sample from the flux space of a metabolic network or from a general convex polytope.

 The following script shows how you could sample steady states of a metabolic network with dingo. To initialize a metabolic network object you have to provide the path to the `json` file as those in [BiGG](http://bigg.ucsd.edu/models) dataset or the `mat` file (using the `matlab` wrapper in folder `/ext_data` to modify a standard `mat` file of a model as those in BiGG dataset):

```python
from dingo import MetabolicNetwork, PolytopeSampler

model = MetabolicNetwork.from_json('path/to/model_file.json')
sampler = PolytopeSampler(model)
steady_states = sampler.generate_steady_states()
```

The output variable `steady_states` is a `numpy` array that contains the steady states of the model column-wise. You could ask from the `sampler` for more statistical guarantees on sampling,  

```python
steady_states = sampler.generate_steady_states(ess=2000, psrf = True)
```

The `ess` stands for the effective sample size (ESS) (default value is `1000`) and the `psrf` is a flag to request an upper bound equal to 1.1 for the value of the  *potential scale reduction factor* of each marginal flux (default option is `False`).  

You could also ask for parallel MMCS algorithm,

```python
steady_states = sampler.generate_steady_states(ess=2000, psrf = True, 
                                               parallel_mmcs = True, num_threads = 2)
```

The default option is to run the sequential [Multiphase Monte Carlo Sampling algorithm](https://arxiv.org/abs/2012.05503) (MMCS) algorithm.  

**Tip**: After the first run of MMCS algorithm the polytope stored in object `sampler` is usually more rounded than the initial one. Thus, the function `generate_steady_states()` becomes more efficient from run to run.  

#### Fast and slow mode

If you have installed successfully the `gurobi` library, dingo turns to the *fast mode* by default. To set a certain mode you could use the following member functions,

```python
sampler = polytope_sampler(model)

#set fast mode to use gurobi library
sampler.set_fast_mode()
#set slow mode to use scipy functions 
sampler.set_slow_mode()
```



### Apply FBA and FVA methods

To apply FVA and FBA methods you have to use the class `metabolic_network`,

```python
from dingo import MetabolicNetwork

model = MetabolicNetwork('path/to/model_file')
fva_output = model.fva()

min_fluxes = fva_output[0]
max_fluxes = fva_output[1]
max_biomass_flux_vector = fva_output[2]
max_biomass_objective = fva_output[3]
```

The output of FVA method is tuple that contains `numpy` arrays. The vectors `min_fluxes` and `max_fluxes` contains the minimum and the maximum values of each flux. The vector `max_biomass_flux_vector` is the optimal flux vector according to the biomass objective function and `max_biomass_objective` is the value of that optimal solution.  

To apply FBA method,

```python
fba_output = model.fba()

max_biomass_flux_vector = fba_output[0]
max_biomass_objective = fba_output[1]
```

while the output vectors are the same with the previous example.   



### Set the restriction in the flux space

FVA and FBA,  restrict the flux space to the set of flux vectors that have an objective value equal to the optimal value of the function. dingo allows for a more  relaxed option where you could ask for flux vectors that have an objective value equal to at least a percentage of the optimal value,

```python
model.set_opt_percentage(90)
fva_output = model.fva()

# the same restriction in the flux space holds for the sampler
sampler = polytope_sampler(model)
steady_states = sampler.generate_steady_states()
```

The default percentage is `100%`.



### Change the objective function

You could also set an alternative objective function. For example, to maximize the 1st reaction of the model,

```python
n = model.num_of_reactions()
obj_fun = np.zeros(n)
obj_fun[0] = 1
model.biomass_function(obj_fun)

# apply FVA using the new objective function
fva_output = model.fva()
# sample from the flux space by restricting 
# the fluxes according to the new objective function
sampler = polytope_sampler(model)
steady_states = sampler.generate_steady_states()
```



### Plot flux marginals

The generated steady states can be used to estimate the marginal density function of each flux. You can plot the histogram using the samples,

```python
from dingo import plot_histogram

model = MetabolicNetwork('path/to/e_coli_core.json')
sampler = PolytopeSampler(model)
steady_states = sampler.generate_steady_states(ess = 3000)

# plot the histogram for the 14th reaction in e-coli (ACONTa)
reactions = model.reactions
plot_histogram(
        steady_states[13],
        reactions[13],
        n_bins = 60,
        )
```

The default number of bins is 60. dingo uses the package `matplotlib` for plotting.

![histogram](../doc/e_coli_aconta.png)

