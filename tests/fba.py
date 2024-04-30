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
        res = model._fba()
        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

        biomass_function = model.reactions[model.biomass_index]
        fba_df = model.fba()
        self.assertTrue(abs(fba_df.loc[biomass_function]["fluxes"] - 0.8739215067486387) < 1e-03)

    def test_fba_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat)
        model.set_slow_mode()
        res = model._fba()
        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

        biomass_function = model.reactions[model.biomass_index]
        fba_df = model.fba()
        self.assertTrue(abs(fba_df.loc[biomass_function]["fluxes"] - 0.8739215067486387) < 1e-03)

    def test_fba_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)
        model.set_slow_mode()
        res = model._fba()
        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)

        biomass_function = model.reactions[model.biomass_index]
        fba_df = model.fba()
        self.assertTrue(abs(fba_df.loc[biomass_function]["fluxes"] - 0.8739215067486387) < 1e-03)

    def test_modify_medium(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)

        # Check if script is running in GitHub action
        if os.getenv('CI', 'false').lower() == 'true':
            model.set_slow_mode()

        initial_medium = model.medium
        initial_fba = model._fba()[-1]
        initial_fba_df = model.fba()

        # Original indices of the exchange reactions
        e_coli_core_medium_compound_indices = {
            "EX_co2_e" : 46,
            "EX_glc__D_e" : 51,
            "EX_h_e" : 54,
            "EX_h2o_e" : 55,
            "EX_nh4_e" : 58,
            "EX_o2_e" : 59,
            "EX_pi_e" : 60
        }

        # Edit and assign new medium to model
        new_media = initial_medium.copy()
        new_media["EX_glc__D_e"] = 35
        new_media["EX_o2_e"] = 0.5
        model.medium = new_media

        # Check if indices are affected
        updated_media = model.medium
        updated_medium_indices = {}
        for reac in updated_media:
            updated_medium_indices[reac] = model.reactions.index(reac)

        glc_index = model.reactions.index("EX_glc__D_e")
        o2_index = model.reactions.index("EX_o2_e")

        self.assertTrue(updated_medium_indices == e_coli_core_medium_compound_indices)
        self.assertTrue(model.lb[glc_index] == -35 and model.lb[o2_index] == -0.5)

        # Check if optimal value is affected
        self.assertTrue(
            abs((model._fba()[-1] - initial_fba) - 0.1172) <= 1e-03
        )


        biomass_function = model.reactions[model.biomass_index]
        fba_df = model.fba()
        self.assertTrue(
            abs((fba_df.loc[biomass_function]["fluxes"] - initial_fba_df.loc[biomass_function]["fluxes"]) - 0.1172) <= 1e-03
        )


if __name__ == "__main__":
    unittest.main()
