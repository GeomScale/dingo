# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Vissarion Fisikopoulos

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import sys
import scipy
import numpy as np
from dingo import MetabolicNetwork
from dingo.scaling import gmscale
from dingo.pyoptinterface_based_impl import set_default_solver


class TestScaling(unittest.TestCase):
    
    def test_scale_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)

        json_res = gmscale(model.S, 0.99)

        self.assertTrue(abs(scipy.linalg.norm(json_res[0]) - 15.285577732002883) < 1e-03)
        self.assertTrue(abs(scipy.linalg.norm(json_res[1]) - 23.138373030721855) < 1e-03)

    def test_scale_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"

        model = MetabolicNetwork.from_mat(input_file_mat)

        mat_res = gmscale(model.S, 0.99)

        self.assertTrue(abs(scipy.linalg.norm(mat_res[0]) - 15.285577732002883) < 1e-03)
        self.assertTrue(abs(scipy.linalg.norm(mat_res[1]) - 23.138373030721855) < 1e-03)

    def test_scale_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"

        model = MetabolicNetwork.from_sbml(input_file_sbml)

        sbml_res = gmscale(model.S, 0.99)

        self.assertTrue(abs(scipy.linalg.norm(sbml_res[0]) - 15.285577732002883) < 1e-03)
        self.assertTrue(abs(scipy.linalg.norm(sbml_res[1]) - 23.138373030721855) < 1e-03)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
