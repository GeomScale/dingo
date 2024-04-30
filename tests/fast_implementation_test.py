# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.gurobi_based_implementations import fast_inner_ball

class TestFastMethods(unittest.TestCase):

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

        model = MetabolicNetwork.from_json(input_file_json)
        model.set_fast_mode()

        res = model._fva()

        self.assertTrue(abs(res[3] - 0.8739215069684305) < 1e-08)
        self.assertEqual(res[0].size, 95)
        self.assertEqual(res[1].size, 95)

        fva_df = model.fva()
        biomass_function = model.reactions[model.biomass_index]
        self.assertTrue(fva_df.loc[biomass_function]["maximum"] - fva_df.loc[biomass_function]["minimum"] <= 1e-03)
        self.assertEqual(fva_df.shape, (95,2))


    def test_ecoli_to_full_dimensional_polytope(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)
        model.set_fast_mode()

        sampler = PolytopeSampler(model)
        sampler.set_fast_mode()

        steady_states = sampler.generate_steady_states()

        self.assertEqual(sampler.A.shape[0], 26)
        self.assertEqual(sampler.A.shape[1], 24)
        self.assertEqual(steady_states.shape[0], 95)


    def test_fast_fba(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)
        model.set_fast_mode()
        res = model._fba()
        self.assertTrue(abs(res[1] - 0.8739215069684305) < 1e-08)

        fba_df = model.fba()
        biomass_function = model.reactions[model.biomass_index]
        self.assertTrue(abs(fba_df.loc[biomass_function]["fluxes"] - 0.8739215069684305) < 1e-08)


if __name__ == "__main__":
    unittest.main()
