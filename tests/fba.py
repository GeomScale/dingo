# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
from dingo import MetabolicNetwork

class TestFba(unittest.TestCase):
    
    def test_fba_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json(input_file_json)
        model.set_slow_mode()
        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_fba_mat(self):
        
        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat)
        model.set_slow_mode()

        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_fba_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)
        model.set_slow_mode()

        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_modify_medium(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)
        model.set_slow_mode()

        initial_medium = model.medium
        initial_fba = model.fba()[-1]
        glc_index = model.reactions.index("EX_glc__D_e")

        new_media = initial_medium.copy()
        new_media["EX_glc__D_e"] = 1.5

        model.medium = new_media

        if model.lb[glc_index] != -1.5:
            self.assertTrue(model.lb[glc_index] == -1.5)

        self.assertTrue(initial_fba - model.fba()[-1] > 0)


if __name__ == "__main__":
    unittest.main()
