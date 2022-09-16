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
[![Chat](https://badges.gitter.im/boostorg/geometry.png)](https://gitter.im/GeomScale/community?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

## Installation

To load the submodules that dingo uses, run

````unix
git submodule update --init
````

You will need to download and unzip the Boost library:
```
wget -O boost_1_76_0.tar.bz2 https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2
tar xjf boost_1_76_0.tar.bz2
rm boost_1_76_0.tar.bz2
```

Then, you need to install the dependencies for the PySPQR library; for Debian/Ubuntu Linux, run

```
apt-get install libsuitesparse-dev
```

To install the Python dependencies, install [Poetry](https://python-poetry.org/). Then, run
```
poetry shell
poetry install
```

To exploit the fast implementations of dingo, you have to install the [Gurobi solver](https://www.gurobi.com/). Run

```
pip3 install -i https://pypi.gurobi.com gurobipy
```

Then, you will need a [license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). For more information, we refer to the Gurobi [download center](https://www.gurobi.com/downloads/).  




## Unit tests

Now, you can run the unit tests by the following commands:  
```
python3 tests/fba.py
python3 tests/full_dimensional.py
python3 tests/max_ball.py
python3 tests/scaling.py
```

If you have installed Gurobi successfully, then run
```
python3 tests/fast_implementation_test.py
```



## Documentation

Read the [dingo documentation](https://github.com/GeomScale/dingo/tree/develop/doc)

You can also have a look at our [Google Colab notebook](https://colab.research.google.com/github/GeomScale/dingo/blob/develop/tutorials/dingo_tutorial.ipynb) 
on how to use `dingo`.
