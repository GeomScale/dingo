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
from dingo.nullspace import nullspace_sparse
from dingo.scaling import gmscale, apply_scaling, remove_almost_redundant_facets
from dingo import HPolytope


class TestStringMethods(unittest.TestCase):
    def test_fast_max_bal_computation(self):

        m = 2
        n = 5

        A = np.zeros((2 * n, n), dtype="float")
        A[0:n] = np.eye(n)
        A[n:] -= np.eye(n, n, dtype="float")
        b = np.ones(2 * n, dtype="float")

        max_ball = fast_inner_ball(A, b)

        self.assertTrue(abs(max_ball[1] - 1) < 1e-08)

    def test_fast_fva(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]
        biomass_index = e_coli_network[5]
        biomass_function = e_coli_network[6]

        fva_res = fast_fva(lb, ub, S, biomass_function)

        A = fva_res[0]
        b = fva_res[1]
        Aeq = fva_res[2]
        beq = fva_res[3]

        self.assertEqual(Aeq.shape, (80, 95))
        self.assertEqual(A.shape, (191, 95))
        self.assertEqual(beq.size, 80)
        self.assertEqual(b.size, 191)
    
    def test_ecoli_to_full_dimensional_polytope(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]
        biomass_index = e_coli_network[5]
        biomass_function = e_coli_network[6]

        fva_res = fast_fva(lb, ub, S, biomass_function)

        A = fva_res[0]
        b = fva_res[1]
        Aeq = fva_res[2]
        beq = fva_res[3]

        nullspace_res_sparse = nullspace_sparse(Aeq, beq)

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

        self.assertEqual(A.shape, (175, 24))
        self.assertEqual(b.size, 175)

        max_ball = fast_inner_ball(A, b)

        self.assertTrue(abs(max_ball[1] - 0.00011558487362508822) < 1e-05)
        self.assertTrue(abs(scipy.linalg.norm(max_ball[0]) - 75962.67878547237) < 1e-05)

        p = HPolytope(A,b)

        self.assertEqual(p.dimension(), 24)

        p.fast_mmcs()

    def test_fast_fba(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        e_coli_network = read_json_file(input_file_json)

        lb = e_coli_network[0]
        ub = e_coli_network[1]
        S = e_coli_network[2]
        biomass_index = e_coli_network[5]
        biomass_function = e_coli_network[6]

        res = fast_fba(lb, ub, S, biomass_function)

        self.assertTrue(abs(res[1] - 0.8739215069684305) < 1e-08)


if __name__ == "__main__":
    unittest.main()
