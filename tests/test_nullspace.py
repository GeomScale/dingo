#!/usr/bin/env python3
import os
import json
import numpy as np
from src.fva import slow_fva
from src.loading_models import read_json_file, read_mat_file
from src.nullspace import nullspace_dense


current_directory = os.getcwd()
input_file_json = current_directory +  '/ext_data/e_coli_core.json'


e_coli_network = volestipy.read_json_file(input_file_json)

A = e_coli_network[0]
b = e_coli_network[1]
Aeq = e_coli_network[2]
beq = e_coli_network[3]

fva_res = slow_fva(A, b, Aeq, beq)

A = fva_res[0]
b = fva_res[1]
Aeq = fva_res[2]
beq = fva_res[3]

nullspace_res_dense = nullspace_dense(Aeq, beq)
#nullspace_res_sparse = nullspace_sparse(Aeq, beq)

print("matrix of nullspace from dense N shape: \n", nullspace_res_dense[0].shape)
#print("matrix of nullspace from sparse N shape: \n", nullspace_res_sparse[0].shape)

