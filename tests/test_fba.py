#!/usr/bin/python3.6

import numpy as np
from dingo.loading_models import read_json_file
from dingo.fba import slow_fba 


current_directory = os.getcwd()
input_file_json = current_directory +  '/ext_data/e_coli_core.json'

print("\n importing e_coli model... \n")
e_coli_network = read_json_file(input_file_json)

lb = e_coli_network[0]
ub = e_coli_network[1]
S = e_coli_network[2]

n = S.shape[1]

obj_fun =  np.ones(n, dtype=np.float)
print("\n this is the objective function: \n")
print(obj_fun)

res = slow_fba(lb, ub, S, obj_fun)
print("FBA optimal solution:")
print(res[0])
print("FBA optimal value:")
print(res[1])

