In the `community_models` directory, you will find models in `.json` format (`.mat` are also supported) 
that will be used as input for the `CommunityMetabolicNetwork` class of `dingo`. 

The class will get all the `.json` files of a directory (or `.mat` accordingly) and build a united model
from all of them. 

In our example case, we use two BIGG models `e_coli_core` and `` as they are among those with the lowest number of reactions, 
so not to have long computation time for it. 
Even in that case, the computational time needed is not negligible, especially if `gurobi` is not available. 
(more than 10 minutes)

Here is how to use `dingo` for a community. 

```python
model         = dingo.CommunityMetabolicNetwork.buildModelList(dir, "json")
sampler       = dingo.CommunityPolytopeSampler(model)
steady_states = sampler.generate_steady_states()
```
