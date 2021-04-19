# dingo

A python library for metabolic networks sampling and analysis.  

## Installation

To load the submodules that dingo uses run: `git submodule update --init`.  

You will need to download and unzip the boost library:
```
wget -O boost_1_75_0.tar.bz2 https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.bz2
tar xjf boost_1_75_0.tar.bz2
rm boost_1_75_0.tar.bz2
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

