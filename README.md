<p align="center"><img src="doc/logo/dingo.jpg" width="260" height="260"></p>

**dingo** is a Python package that analyzes metabolic networks.
It relies on high dimensional sampling with Markov Chain Monte Carlo (MCMC)
methods and fast optimization methods to analyze the possible states of a
metabolic network. To perform MCMC sampling, `dingo` relies on the `C++` library
[volesti](https://github.com/GeomScale/volume_approximation), which provides
several algorithms for sampling convex polytopes.
`dingo` also performs two standard methods to analyze the flux space of a
metabolic network, namely Flux Balance Analysis and Flux Variability Analysis.

`dingo` is part of [GeomScale](https://geomscale.github.io/) project.

[![unit-tests](https://github.com/GeomScale/dingo/workflows/dingo-ubuntu/badge.svg)](https://github.com/GeomScale/dingo/actions?query=workflow%3Adingo-ubuntu)
[![Tutorial In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeomScale/dingo/blob/develop/tutorials/dingo_tutorial.ipynb)
[![Chat](https://badges.gitter.im/geomscale.png)](https://gitter.im/GeomScale/community?utm_source=share-link&utm_medium=link&utm_campaign=share-link)


## Installation 

### LP solver (optional, probably better performance)

`dingo` makes use of [`pyoptinterface`](https://metab0t.github.io/PyOptInterface/) to interface with a range of linear programming solvers. 

The default solver is [`highs`](https://highs.dev/#get-started). 

However, one may switch to other solvers that `pyoptinterface` supports, for example the commonly used [`gurobi`](https://www.gurobi.com/). 
Yet, in that case a Gurobi license is required. 

> **Get a Gurobi license**
> 
> If you are affiliated in an academic insitute, you can generate a **free academic license**.
> 
> First, register and/or login to your [Gurobi account](https://portal.gurobi.com/iam/login/), and 
> 
> * if you are about to use `dingo` as a container, get a [**Web License Service (WLS) academic license**](https://support.gurobi.com/hc/en-us/articles/13210193318033-What-is-an-Academic-WLS-license)
> * otherwise, you should go for the typical [**free academic license**](https://www.gurobi.com/academics)
>  
> 🔴 In both cases, make sure you are connected to the internet of an academic institution.


### Installation (on Linux)

**Note:** Python version should be 3.8.x. You can check this by running the following command in your terminal:
```bash
python --version
```

If you have a different version of Python installed, you'll need to install it ([start here](https://linuxize.com/post/how-to-install-python-3-8-on-ubuntu-18-04/))
and update-alternatives ([start here](https://linuxhint.com/update_alternatives_ubuntu/)).

Clone the `dingo` repo by 

```
git clone https://github.com/GeomScale/dingo.git
```


and load the submodules that `dingo` uses:

````bash
cd dingo
git submodule update --init
````

You will then need to download and unzip the [Boost C++](https://www.boost.org/) library:
```
wget -O boost_1_76_0.tar.bz2 https://archives.boost.io/release/1.76.0/source/boost_1_76_0.tar.bz2
tar xjf boost_1_76_0.tar.bz2
rm boost_1_76_0.tar.bz2
```

You will also need to download and unzip the [`lpsolve`](https://lpsolve.sourceforge.net/5.5/) library:
```
wget https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz
tar xzvf lp_solve_5.5.2.11_source.tar.gz
rm lp_solve_5.5.2.11_source.tar.gz
```

Then, you need to install the dependencies for the [PySPQR](https://github.com/yig/PySPQR) library; 
for this, you will most likely need `sudo` rights:

```bash
sudo apt-get update -y
sudo apt-get install -y libsuitesparse-dev
```

To install the Python dependencies, `dingo` is using [Poetry](https://python-poetry.org/),
```
curl -sSL https://install.python-poetry.org | python - --version 1.3.2
poetry shell
poetry install
```

otherwise, you may try:

```
python setup.py install --user
```


Last, in case you are about to use Gurobi, remember to install the Python interface of Gurobi, [`gurobipy`](https://www.gurobi.com/resources/faq/gurobipy):

```
pip install -i https://pypi.gurobi.com gurobipy
```



## Using `dingo` as a Docker container

To use `dingo` as a container, you need to [install Docker](https://docs.docker.com/engine/install/), 
or [Docker desktop](https://docs.docker.com/desktop/), first.

Then you can clone the `dingo` repo and build its Docker image:

```
git clone https://github.com/GeomScale/dingo.git
cd dingo 
docker build -f Dockerfile -t dingo .
```

Once the image is built, you may run:

```
docker run --rm -it -v <path_to_your_model>:/data dingo
```

or, if you are using Gurobi, you may run:

```
docker run --rm -it -v <path_to_WLS_license>:/opt/gurobi/gurobi.lic -v <path_to_your_model>:/data dingo
```

> **Remember** for this use need the WLS Gurobi license. 
> 
> This would look something like this:
>
> ```
> # Gurobi WLS license file
> # Your credentials are private and should not be shared or copied to public repositories.
> # Visit https://license.gurobi.com/manager/doc/overview for more information.
> WLSACCESSID=d5419c87-0d36-4a93-9385-773f5483b3c1
> WLSSECRET=afa5d95f-ad0b-4a38-9550-a8913aacb7c0
> LICENSEID=000000
>```



## Unit tests

Now, you can run the unit tests by the following commands (with the default solver `highs`):
```
python tests/fba.py
python tests/full_dimensional.py
python tests/max_ball.py
python tests/scaling.py
python tests/rounding.py
python tests/sampling.py
```

Or, assuming you have installed Gurobi successfully, or an other `pyoptinterface`-supported solver, you may run:
```
python tests/fba.py gurobi
python tests/full_dimensional.py gurobi
python tests/max_ball.py gurobi
python tests/scaling.py gurobi
python tests/rounding.py gurobi
python tests/sampling.py gurobi
```

## Tutorial

You may check out `dingo`'s main features through a GitHub codespace.
To do this, you may click [here](https://github.com/codespaces/new?repo=GeomScale/dingo&ref=main) and fire a new codespace
by clicking on the "Create codespace" button. 

This will take a few minutes (~5'). 

Once the codespace is ready, you may try to follow the [`dingo_tutorial`](./tutorials/dingo_tutorial.ipynb) Jupyter notebook. 


## Documentation


It quite simple to use dingo in your code. In general, dingo provides two classes:

- `metabolic_network` represents a metabolic network
- `polytope_sampler` can be used to sample from the flux space of a metabolic network or from a general convex polytope.

 The following script shows how you could sample steady states of a metabolic network with dingo. To initialize a metabolic network object you have to provide the path to the `json` file as those in [BiGG](http://bigg.ucsd.edu/models) dataset or the `mat` file (using the `matlab` wrapper in folder `/ext_data` to modify a standard `mat` file of a model as those in BiGG dataset):

```python
from dingo import MetabolicNetwork, PolytopeSampler

model = MetabolicNetwork.from_json('path/to/model_file.json')
sampler = PolytopeSampler(model)
steady_states = sampler.generate_steady_states()
```

`dingo` can also load a model given in `.sbml` format using the following command,

```python
model = MetabolicNetwork.from_sbml('path/to/model_file.sbml')
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


#### Rounding the polytope

`dingo` provides three methods to round a polytope: (i) Bring the polytope to John position by apllying to it the transformation that maps the largest inscribed ellipsoid of the polytope to the unit ball, (ii) Bring the polytope to near-isotropic position by using uniform sampling with Billiard Walk, (iii) Apply to the polytope the transformation that maps the smallest enclosing ellipsoid of a uniform sample from the interior of the polytope to the unit ball.

```python
from dingo import MetabolicNetwork, PolytopeSampler

model = MetabolicNetwork.from_json('path/to/model_file.json')
sampler = PolytopeSampler(model)
A, b, N, N_shift = sampler.get_polytope()

A_rounded, b_rounded, Tr, Tr_shift = sampler.round_polytope(A, b, method="john_position")
A_rounded, b_rounded, Tr, Tr_shift = sampler.round_polytope(A, b, method="isotropic_position")
A_rounded, b_rounded, Tr, Tr_shift = sampler.round_polytope(A, b, method="min_ellipsoid")
```

Then, to sample from the rounded polytope, the user has to call the following static method of PolytopeSampler class,

```python
samples = sample_from_polytope(A_rounded, b_rounded)
```

Last you can map the samples back to steady states,

```python
from dingo import map_samples_to_steady_states

steady_states = map_samples_to_steady_states(samples, N, N_shift, Tr, Tr_shift)
```

#### Other MCMC sampling methods

To use any other MCMC sampling method that `dingo` provides you can use the following piece of code:

```python
sampler = polytope_sampler(model)
steady_states = sampler.generate_steady_states_no_multiphase() #default parameters (method = 'billiard_walk', n=1000, burn_in=0, thinning=1)
```

The MCMC methods that dingo (through `volesti` library) provides are the following: (i) 'cdhr': Coordinate Directions Hit-and-Run, (ii) 'rdhr': Random Directions Hit-and-Run,
(iii) 'billiard_walk', (iv) 'ball_walk', (v) 'dikin_walk', (vi) 'john_walk', (vii) 'vaidya_walk'.



#### Switch the linear programming solver

We use `pyoptinterface` to interface with the linear programming solvers. 
To switch the solver that `dingo` uses, you can use the `set_default_solver` function. 
The default solver is `highs` and you can switch to `gurobi` by running:

```python
from dingo import set_default_solver
set_default_solver("gurobi")
```

You can also switch to other solvers that `pyoptinterface` supports, but we recommend using `highs` or `gurobi`. 
If you have issues with the solver, you can check the `pyoptinterface` [documentation](https://metab0t.github.io/PyOptInterface/getting_started.html).

### Apply FBA and FVA methods

To apply FVA and FBA methods you have to use the class `metabolic_network`,

```python
from dingo import MetabolicNetwork

model = MetabolicNetwork.from_json('path/to/model_file.json')
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
model.objective_function(obj_fun)

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

model = MetabolicNetwork.from_json('path/to/e_coli_core.json')
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

![histogram](./doc/e_coli_aconta.png)

### Plot a copula between two fluxes

The generated steady states can be used to estimate and plot the copula between two fluxes. You can plot the copula using the samples,

```python
from dingo import plot_copula

model = MetabolicNetwork.from_json('path/to/e_coli_core.json')
sampler = PolytopeSampler(model)
steady_states = sampler.generate_steady_states(ess = 3000)

# plot the copula between the 13th (PPC) and the 14th (ACONTa) reaction in e-coli
reactions = model.reactions

data_flux2=[steady_states[12],reactions[12]]
data_flux1=[steady_states[13],reactions[13]]

plot_copula(data_flux1, data_flux2, n=10)
```

The default number of cells is 5x5=25. dingo uses the package `plotly` for plotting.

![histogram](./doc/aconta_ppc_copula.png)


## Citation

Apostolos Chalkis, Vissarion Fisikopoulos, Elias Tsigaridas, Haris Zafeiropoulos, dingo: a Python package for metabolic flux sampling, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae037, [https://doi.org/10.1093/bioadv/vbae037](https://doi.org/10.1093/bioadv/vbae037)


