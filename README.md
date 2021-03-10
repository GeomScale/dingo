# dingo

A python library for metabolic networks sampling and analysis.  

## Installation

### Dependencies

To load the submodules that dingo uses run: `git submodule update --init`.  

For debian/ubuntu linux run: `apt-get install libsuitesparse-dev`.  

To install the python dependencies install [poetry](https://python-poetry.org/). Then, run:  
```
poetry install
```

To exploit the fast implementations of dingo you have to install the [Gurobi solver](https://www.gurobi.com/). Run:  

```
pip3 install -i https://pypi.gurobi.com gurobipy
```

Then, you will need a [license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). For more information we refer to gurobi [download center](https://www.gurobi.com/downloads/).  


### Install *dingo*

After getting the dependencies run:

```
python3 setup.py install --user
```

## Unit tests

You can run the unit tests by the following commands:  
```
python3 tests/unit_tests.py
```

If you have installed gurobi sucesfully then run:  
```
python3 tests/fast_implementation_test.py
```

