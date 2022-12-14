# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Vissarion Fisikopoulos

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
from dingo import MetabolicNetwork, PolytopeSampler

class TestFullDim(unittest.TestCase):

    def test_get_full_dim_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json( input_file_json )
        sampler = self.get_polytope_from_model(model)

        self.assertEqual(sampler.A.shape[0], 175)
        self.assertEqual(sampler.A.shape[1], 24)

    def test_get_full_dim_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml( input_file_sbml )
        sampler = self.get_polytope_from_model( model )

        self.assertEqual(sampler.A.shape[0], 175)
        self.assertEqual(sampler.A.shape[1], 24)


    def test_get_full_dim_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat( input_file_mat )
        sampler = self.get_polytope_from_model( model )

        self.assertEqual(sampler.A.shape[0], 175)
        self.assertEqual(sampler.A.shape[1], 24)

    @staticmethod
    def get_polytope_from_model (met_model):

        met_model.set_slow_mode()
        sampler = PolytopeSampler(met_model)
        sampler.set_slow_mode()
        sampler.get_polytope()

        return sampler

if __name__ == "__main__":
    unittest.main()
