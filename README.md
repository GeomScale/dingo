# dingo

A python library for metabolic networks sampling and analysis.  

## Installation

### Dependencies

You need `cython`, `numpy` and `setuptools`. In `debian` systems you can get then by
```
sudo apt-get install cython libsuitesparse-dev
```

Finally, the [Gurobi solver](https://www.gurobi.com/) needs to be installed as well to exploit fast implementations.
Through the [Download center](https://www.gurobi.com/downloads/) you may download its Python interface along with a license; the latter is a text file called `gurobi.lic`.  

You can follow the inscriptions described [here](https://support.gurobi.com/hc/en-us/articles/360044290292-Installing-Gurobi-for-Python) to get ```gurobipy```. 

Once you complete these steps, make sure that `gurobipy` is now available for your Python environment. 

```
$ python3  
>>> import gurobipy
```

### Load `volesti` submodule
To load `volesti` submodule run:  
 
```
git submodule update --init
```

### Dependencies for PySPQR

```
sudo apt-get install libsuitesparse-dev
pip3 install cffi
```

### Install *dingo*

After getting the dependencies run:

```
python3 setup.py install --user
```

## Run an example


