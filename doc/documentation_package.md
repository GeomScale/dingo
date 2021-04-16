## Use dingo as a package

It quite simple to use dingo in your code. The following script shows how you could sample steady states of a metabolic network with dingo,

```python
from dingo import metabolic_network, polytope_sampler

model = metabolic_network('path/to/model_file')
sampler = polytope_sampler(model)
steady_states = sampler.generate_steady_states()
```

The output variable `steady_states` is a `numpy` array that contains the steady states of the model column-wise. You could ask from the `sampler` for more statistical guarantees on sampling,  

```python
steady_states = sampler.generate_steady_states(ess=2000, psrf = True)
```

The `ess` stands for the effective sample size (ESS) (default value is `1000`) and the `psrf` is a flag to request an upper bound equal to 1.1 for the value of the  *potential scale reduction factor* of each marginal flux (default option is `False`).  

You could also ask for parallel MMMCS algorithm,

```python
steady_states = sampler.generate_steady_states(ess=2000, psrf = True, parallel_mmcs = True, num_threads = 2)
```

The default option is to run the sequential MMCS algorithm.  

 dingo provides two classes:  

- `metabolic_network` represents a metabolic network
- `polytope_sampler` can be used to sample from the flux space of a metabolic network or from a general convex polytope.

### Perform FBA and FVA methods

To apply FVA and FBA methods you have to use the class `metabolic_network`,

```python
from dingo import metabolic_network

model = metabolic_network('path/to/model_file')
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

max_biomass_flux_vector = fva_output[0]
max_biomass_objective = fva_output[1]
```

while the output vectors are the same with the previous example.  