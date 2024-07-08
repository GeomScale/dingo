
from cobra.io import load_json_model
from dingo.preprocess import PreProcess
import unittest
import numpy as np


class TestPreprocess(unittest.TestCase):

    def test_preprocess(self):

        model = load_json_model("ext_data/e_coli_core.json")
        obj = PreProcess(model)
        
        
        # find reaction ids from the loaded model
        initial_reactions_ids = []    
        for reaction in model.reactions:
            reaction_id = reaction.id
            initial_reactions_ids.append(reaction_id)        


        # calculate the count of removed reactions with extend set to 0
        removed_reactions, dingo_model = obj.reduce(extend=0)
        removed_reactions_count = len(removed_reactions)
        self.assertTrue( 46 - removed_reactions_count == 0 )        

        # calculate the count of removed reactions with extend set to 0 from the dingo model       
        dingo_removed_reactions = np.sum((dingo_model[0] == 0) & (dingo_model[1] == 0))
        self.assertTrue( 46 - dingo_removed_reactions == 0 )
        
        # calculate the count of reactions with bounds equal to 0 with extend set to 0
        zero_flux_count = 0
        for reaction_id in initial_reactions_ids:
            bounds = model.reactions.get_by_id(reaction_id).bounds
            if ((bounds[0] == 0) and (bounds[1] == 0)):
                zero_flux_count += 1
            
        self.assertTrue( 46 - zero_flux_count == 0 )
        
        fba_solution = model.optimize()   
        self.assertTrue(abs(fba_solution.objective_value - 0.8739215067486387) < 1e-03)



        # calculate the count of removed reactions with extend set to 1
        removed_reactions, dingo_model = obj.reduce(extend=1)
        removed_reactions_count = len(removed_reactions)
        self.assertTrue( 47 - removed_reactions_count == 0 )
        
        # calculate the count of removed reactions with extend set to 1 from the dingo model       
        dingo_removed_reactions = np.sum((dingo_model[0] == 0) & (dingo_model[1] == 0))
        self.assertTrue( 47 - dingo_removed_reactions == 0 )

                   
        # calculate the count of reactions with bounds equal to 0 with extend set to 1        
        zero_flux_count = 0
        for reaction_id in initial_reactions_ids:
            bounds = model.reactions.get_by_id(reaction_id).bounds
            if ((bounds[0] == 0) and (bounds[1] == 0)):
                zero_flux_count += 1
            
        self.assertTrue( 47 - zero_flux_count == 0 )
      
        fba_solution = model.optimize()   
        self.assertTrue(abs(fba_solution.objective_value - 0.8739215067486387) < 1e-03)
        

if __name__ == "__main__":
    unittest.main()

 