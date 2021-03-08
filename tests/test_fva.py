#!/usr/bin/python3.6
import os
import json
import numpy as np
import scipy.sparse as sp
from dingo.loading_models import read_json_file
from dingo.fva import slow_fva, fast_fva


current_directory = os.getcwd()
input_file_json = current_directory +  '/ext_data/e_coli_core.json'

print("\n importing e_coli model... \n")
e_coli_network = read_json_file(input_file_json)

lb = e_coli_network[0]
ub = e_coli_network[1]
S = e_coli_network[2]

res = fast_fva(lb, ub, S)
print("shape of A is:")
print(res[0].shape[0], res[0].shape[1])
print("new b size is:")
print(res[1].size)
print("new Aeq size is:")
print(res[2].shape[0], res[2].shape[1])
print("new beq size is:")
print(res[3].size)
print("minimum values of fluxes:")
print(res[4])
print("maximum values of fluxes:")
print(res[5])


