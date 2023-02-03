# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis
# Copyright (c) 2022 Vissarion Fisikopoulos
# Copyright (c) 2022 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.gurobi_based_implementations import fast_inner_ball

class TestSampling(unittest.TestCase):

    def test_rounding_min_ellipsoid(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json( input_file_json )
        sampler = PolytopeSampler(model)

        A, b, N, N_shift = sampler.get_polytope()

        A_rounded, b_rounded, Tr, Tr_shift = sampler.round_polytope(A, b, method="min_ellipsoid")

        self.assertTrue( A_rounded.shape[0] == 26 )
        self.assertTrue( A_rounded.shape[1] == 24 )

        self.assertTrue( b.size == 26 )
        self.assertTrue( N_shift.size == 95 )
        self.assertTrue( b_rounded.size == 26 )
        self.assertTrue( Tr_shift.size == 24 )
        

        self.assertTrue( N.shape[0] == 95 )
        self.assertTrue( N.shape[1] == 24 )

        self.assertTrue( Tr.shape[0] == 24 )
        self.assertTrue( Tr.shape[1] == 24 )

        samples = sampler.sample_from_polytope_no_multiphase(
            A_rounded, b_rounded, method = 'billiard_walk', n=1000, burn_in=10, thinning=1
        )

        temp = np.full((samples.shape[0], samples.shape[1]), Tr_shift)
        Tr_samples = Tr.dot(samples.T) + temp.T
        all_points_in = True
        for i in range(Tr_samples.shape[1]):
            if np.any(A.dot(Tr_samples[:,i]) - b > 1e-05):
                all_points_in = False
                break
        
        self.assertTrue( all_points_in )

    def test_rounding_john_position(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json( input_file_json )
        sampler = PolytopeSampler(model)

        A, b, N, N_shift = sampler.get_polytope()

        A_rounded, b_rounded, Tr, Tr_shift = sampler.round_polytope(A, b, method="john_position")

        self.assertTrue( A_rounded.shape[0] == 26 )
        self.assertTrue( A_rounded.shape[1] == 24 )

        self.assertTrue( b.size == 26 )
        self.assertTrue( N_shift.size == 95 )
        self.assertTrue( b_rounded.size == 26 )
        self.assertTrue( Tr_shift.size == 24 )
        

        self.assertTrue( N.shape[0] == 95 )
        self.assertTrue( N.shape[1] == 24 )

        self.assertTrue( Tr.shape[0] == 24 )
        self.assertTrue( Tr.shape[1] == 24 )

        samples = sampler.sample_from_polytope_no_multiphase(
            A_rounded, b_rounded, method = 'billiard_walk', n=1000, burn_in=10, thinning=1
        )

        temp = np.full((samples.shape[0], samples.shape[1]), Tr_shift)
        Tr_samples = Tr.dot(samples.T) + temp.T
        all_points_in = True
        for i in range(Tr_samples.shape[1]):
            if np.any(A.dot(Tr_samples[:,i]) - b > 1e-05):
                all_points_in = False
                break
        
        self.assertTrue( all_points_in )




if __name__ == "__main__":
    unittest.main()