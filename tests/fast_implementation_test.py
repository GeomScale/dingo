# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import scipy
import numpy as np
from dingo.loading_models import read_json_file
from dingo.gurobi_based_implementations import fast_inner_ball, fast_fba, fast_fva

class TestStringMethods(unittest.TestCase):

    def test_fast_max_bal_computation(self):

        m = 2
        n = 5

        A = np.zeros((2*n, n), dtype='float')
        A[0:n] = np.eye(n)
        A[n:] -=  np.eye(n,n, dtype='float')
        b = np.ones(2*n, dtype='float')
        
        max_ball = fast_inner_ball(A,b)

        self.assertTrue(abs(max_ball[1] - 1) < 1e-10)


    def test_fast_fva(self):

        current_directory = os.getcwd()
        input_file_json = current_directory +  '/ext_data/e_coli_core.json'

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]

        fva_res = fast_fva(lb, ub, S)

        A = fva_res[0]
        b = fva_res[1]
        Aeq = fva_res[2]
        beq = fva_res[3]

        self.assertEqual(Aeq.shape, (80, 95))
        self.assertEqual(A.shape, (190, 95))
        self.assertEqual(beq.size, 80)
        self.assertEqual(b.size, 190)


    def test_fast_fba(self):

        current_directory = os.getcwd()
        input_file_json = current_directory +  '/ext_data/e_coli_core.json'

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]

        n = S.shape[1]

        obj_fun =  np.ones(n, dtype='float')
        res = fast_fba(lb, ub, S, obj_fun)

        self.assertTrue(abs(res[1] - 3103.2200000000003) < 1e-10)


if __name__ == '__main__':
    unittest.main()
