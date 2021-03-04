#!/usr/bin/env python3
import os
import json
import numpy as np
from dingo.loading_models import read_json_file, read_mat_file
from dingo.nullspace import nullspace_sparse
from dingo.scaling import gmscale, apply_scaling, remove_almost_redundant_facets
from dingo.fva import slow_fva 
from dingo.inner_ball import slow_inner_ball


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

nullspace_res_sparse = nullspace_sparse(Aeq, beq)

N = nullspace_res_sparse[0]
N_shift = nullspace_res_sparse[1]

product = np.dot(A, N_shift)
b = np.subtract(b, product)
A = np.dot(A, N)

res = remove_almost_redundant_facets(A, b)
A = res[0]
b = res[1]

res = gmscale(A, 0, 0.99):
res = apply_scaling(A, b, res[0], res[1])
A = res[0]
b = res[1]

res = remove_almost_redundant_facets(A, b)
A = res[0]
b = res[1]

print("shape of the matrix of the full dimensional polytope, A: \n", A.shape)
print("shape of the vector of the full dimensional polytope, b: \n", b.shape)

max_ball = slow_inner_ball(A, b)
print("The center point of the maximum inscribed ball:")
print(max_ball[0])
print("The radius of the maximum inscribed ball:")
print(max_ball[1])
