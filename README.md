<p align="center"><img src="doc/logo/dingo.jpg" width="260" height="260"></p>

**dingo** is a python package that analyzes metabolic networks.
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

To load the submodules that dingo uses run,

````unix
git submodule update --init
````

You will need to download and unzip the boost library,
```
wget -O boost_1_76_0.tar.bz2 https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2
tar xjf boost_1_76_0.tar.bz2
rm boost_1_76_0.tar.bz2
```

Then, you need to install the dependencies for PySPQR library; for debian/ubuntu linux run:

```
apt-get install libsuitesparse-dev
```

To install the python dependencies install [poetry](https://python-poetry.org/). Then, run:  
```
poetry shell
poetry install
```

To exploit the fast implementations of dingo you have to install the [Gurobi solver](https://www.gurobi.com/). Run:  

```
pip3 install -i https://pypi.gurobi.com gurobipy
```

Then, you will need a [license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). For more information we refer to gurobi [download center](https://www.gurobi.com/downloads/).  




## Unit tests

Now, you can run the unit tests by the following commands:  
```
python3 tests/unit_tests.py
```

If you have installed gurobi sucesfully then run:  
```
python3 tests/fast_implementation_test.py
```



## Documentation

Read [dingo documentation](https://github.com/GeomScale/dingo/tree/develop/doc)

You can also have a look at our [Google Collab notebook](https://colab.research.google.com/drive/1AdXCo6tMBV4lPDYWWOXhzcM0wg30OOIx?usp=sharing) 
on how to use `dingo`.
