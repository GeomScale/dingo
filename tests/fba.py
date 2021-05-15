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


class TestFba(unittest.TestCase):
    
    def test_ecoli(self):

        current_directory = os.getcwd()
        input_file_json = current_directory + "/ext_data/e_coli_core.json"

        model = MetabolicNetwork.from_json(input_file_json)
        model.set_slow_mode()

        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

if __name__ == "__main__":
    unittest.main()
