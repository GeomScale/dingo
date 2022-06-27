# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import scipy
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.inner_ball import slow_inner_ball
from dingo.scaling import gmscale


class TestDingoMethods(unittest.TestCase):
    def test_max_bal_computation(self):

        m = 2
        n = 5

        A = np.zeros((2 * n, n), dtype="float")
        A[0:n] = np.eye(n)
        A[n:] -= np.eye(n, n, dtype="float")
        b = np.ones(2 * n, dtype="float")

        max_ball = slow_inner_ball(A, b)

        self.assertTrue(abs(max_ball[1] - 1) < 1e-04)

    def test_ecoli_to_full_dimensional_polytope(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)
        model.set_slow_mode()

        sampler = PolytopeSampler(model)
        sampler.set_slow_mode()
        sampler.get_polytope()

        self.assertEqual(sampler.A.shape[0], 175)
        self.assertEqual(sampler.A.shape[1], 24)

    def test_fba(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)
        model.set_slow_mode()

        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_scaling(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)

        res = gmscale(model.S, 0.99)

        self.assertTrue(abs(scipy.linalg.norm(res[0]) - 15.285577732002883) < 1e-03)
        self.assertTrue(abs(scipy.linalg.norm(res[1]) - 23.138373030721855) < 1e-03)


if __name__ == "__main__":
    unittest.main()
