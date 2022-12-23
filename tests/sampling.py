# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis
# Copyright (c) 2022 Vissarion Fisikopoulos
# Copyright (c) 2022 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
from dingo import MetabolicNetwork, PolytopeSampler


class TestSampling(unittest.TestCase):

    def test_sample_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json( input_file_json )
        sampler = PolytopeSampler(model)

        steady_states = sampler.generate_steady_states(ess = 1000, psrf = True) 

        self.assertTrue( steady_states.shape[0] == 95 )
        self.assertTrue( abs( steady_states[12].mean()  - 2.504 ) < 1e-03 )


    def test_sample_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat)
        sampler = PolytopeSampler(model)

        steady_states = sampler.generate_steady_states(ess = 1000, psrf = True) 

        self.assertTrue( steady_states.shape[0] == 95 )
        self.assertTrue( abs( steady_states[12].mean()  - 2.504 ) < 1e-03 )


    def test_sample_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml( input_file_sbml )
        sampler = PolytopeSampler(model)

        steady_states = sampler.generate_steady_states(ess = 1000, psrf = True) 

        self.assertTrue( steady_states.shape[0] == 95 )
        self.assertTrue( abs( steady_states[12].mean()  - 2.504 ) < 1e-03 )



if __name__ == "__main__":
    unittest.main()
