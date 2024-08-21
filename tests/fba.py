# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import sys
from dingo import MetabolicNetwork
from dingo.pyoptinterface_based_impl import set_default_solver

class TestFba(unittest.TestCase):

    def test_fba_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json(input_file_json)
        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_fba_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat)

        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_fba_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)

        res = model.fba()

        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

    def test_modify_medium(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)

        initial_medium = model.medium
        initial_fba = model.fba()[-1]

        e_coli_core_medium_compound_indices = {
            "EX_co2_e" : 46,
            "EX_glc__D_e" : 51,
            "EX_h_e" : 54,
            "EX_h2o_e" : 55,
            "EX_nh4_e" : 58,
            "EX_o2_e" : 59,
            "EX_pi_e" : 60
        }

        glc_index = model.reactions.index("EX_glc__D_e")
        o2_index = model.reactions.index("EX_o2_e")

        new_media = initial_medium.copy()
        new_media["EX_glc__D_e"] = 1.5
        new_media["EX_o2_e"] = -0.5

        model.medium = new_media

        updated_media = model.medium
        updated_medium_indices = {}
        for reac in updated_media:
            updated_medium_indices[reac] = model.reactions.index(reac)

        self.assertTrue(updated_medium_indices == e_coli_core_medium_compound_indices)

        self.assertTrue(model.lb[glc_index] == -1.5 and model.lb[o2_index] == 0.5)

        self.assertTrue(initial_fba - model.fba()[-1] > 0)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
