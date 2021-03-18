# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import scipy
import numpy as np
from dingo.fva import slow_fva
from dingo.fba import slow_fba
from dingo.loading_models import read_json_file
from dingo.inner_ball import slow_inner_ball
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import gmscale, apply_scaling, remove_almost_redundant_facets

class TestStringMethods(unittest.TestCase):

    def test_max_bal_computation(self):

        m = 2
        n = 5

        A = np.zeros((2*n, n), dtype='float')
        A[0:n] = np.eye(n)
        A[n:] -=  np.eye(n,n, dtype='float')
        b = np.ones(2*n, dtype='float')
        
        max_ball = slow_inner_ball(A,b)

        self.assertTrue(abs(max_ball[1] - 1) < 1e-10)


    def test_ecoli_to_full_dimensional_polytope(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + '/ext_data/e_coli_core.json'

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]

        fva_res = slow_fva(lb, ub, S)

        A = fva_res[0]
        b = fva_res[1]
        Aeq = fva_res[2]
        beq = fva_res[3]

        nullspace_res_dense = nullspace_dense(Aeq, beq)
        nullspace_res_sparse = nullspace_sparse(Aeq, beq)

        self.assertEqual(nullspace_res_dense[0].shape, (95, 24))
        self.assertEqual(nullspace_res_sparse[0].shape, (95, 24))

        N = nullspace_res_sparse[0]
        N_shift = nullspace_res_sparse[1]

        product = np.dot(A, N_shift)
        b = np.subtract(b, product)
        A = np.dot(A, N)

        res = remove_almost_redundant_facets(A, b)
        A = res[0]
        b = res[1]

        res = gmscale(A, 0.99)
        res = apply_scaling(A, b, res[0], res[1])
        A = res[0]
        b = res[1]

        res = remove_almost_redundant_facets(A, b)
        A = res[0]
        b = res[1]

        self.assertEqual(A.shape, (174, 24))
        self.assertEqual(b.size, 174)

        max_ball = slow_inner_ball(A, b)

        self.assertTrue(abs(max_ball[1] - 187.2556753095123) < 1e-04)
        self.assertTrue(abs(scipy.linalg.norm(max_ball[0]) - 37645.55778670023) < 1e-10)


    def test_fba(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + '/ext_data/e_coli_core.json'

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]

        n = S.shape[1]

        obj_fun =  np.ones(n, dtype='float')
        res = slow_fba(lb, ub, S, obj_fun)

        self.assertTrue(abs(res[1] - 3103.219983067629) < 1e-04)
    

    def test_scaling(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + '/ext_data/e_coli_core.json'

        e_coli_network = read_json_file(input_file_json)
        S = e_coli_network[2]

        res = gmscale(S, 0.99)

        self.assertTrue(abs(scipy.linalg.norm(res[0]) - 15.285577732002883) < 1e-10)
        self.assertTrue(abs(scipy.linalg.norm(res[1]) - 23.138373030721855) < 1e-10)


if __name__ == '__main__':
    unittest.main()
