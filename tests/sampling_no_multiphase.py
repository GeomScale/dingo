# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis
# Copyright (c) 2022-2025 Vissarion Fisikopoulos
# Copyright (c) 2022 Haris Zafeiropoulos
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.pyoptinterface_based_impl import set_default_solver

def test_shape_and_non_zero(steady_states, shape0, shape1, testing_class):
    testing_class.assertTrue( steady_states.shape[0] == shape0 )
    testing_class.assertTrue( steady_states.shape[1] == shape1 )
    testing_class.assertTrue( np.any(np.abs(steady_states) > 1e-12) )

def sampling(model, testing_class):
    sampler = PolytopeSampler(model)

    #Gaussian hmc sampling
    steady_states = sampler.generate_steady_states_no_multiphase(method = 'gaussian_hmc_walk', n=1000)
    test_shape_and_non_zero(steady_states, 95, 1000, testing_class)

    #exponential hmc sampling
    steady_states = sampler.generate_steady_states_no_multiphase(method = 'exponential_hmc_walk', n=1000, variance=50)
    test_shape_and_non_zero(steady_states, 95, 1000, testing_class)
    
    #hmc sampling with Gaussian distribution
    steady_states = sampler.generate_steady_states_no_multiphase(method = 'hmc_leapfrog_gaussian', n=1000)
    test_shape_and_non_zero(steady_states, 95, 1000, testing_class)

    #hmc sampling with exponential distribution
    steady_states = sampler.generate_steady_states_no_multiphase(method = 'hmc_leapfrog_exponential', n=1000, variance=50)
    test_shape_and_non_zero(steady_states, 95, 1000, testing_class)

    #steady_states[12].mean() seems to have a lot of discrepancy between experiments, so we won't check the mean for now
    #self.assertTrue( abs( steady_states[12].mean()  - 2.504 ) < 1e-03 )

class TestSampling(unittest.TestCase):

    def test_sample_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json( input_file_json )
        sampling(model, self)

    def test_sample_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat)
        sampling(model, self)

    def test_sample_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml( input_file_sbml )
        sampling(model, self)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
