
from cobra.io import load_json_model
from dingo.preprocess import PreProcess
import unittest
import numpy as np


class TestPreprocess(unittest.TestCase):

    def test_preprocess(self):

        # load cobra model
        cobra_model = load_json_model("ext_data/e_coli_core.json")

        # find reaction ids from the loaded model
        initial_reactions_ids = []    
        for reaction in cobra_model.reactions:
            reaction_id = reaction.id
            initial_reactions_ids.append(reaction_id)


        # call the reduce function from the PreProcess class 
        # with extend=0 to remove reactions from the model        
        obj = PreProcess(cobra_model)  
        removed_reactions, dingo_model = obj.reduce(extend=0)      
        
        # calculate the count of removed reactions with extend set to 0        
        removed_reactions_count = len(removed_reactions)
        self.assertTrue( 46 - removed_reactions_count == 0 )

        # calculate the count of reactions with bounds equal to 0 
        # with extend set to 0 from the dingo model
        dingo_removed_reactions = np.sum((dingo_model.lb == 0) & (dingo_model.ub == 0))
        self.assertTrue( 46 - dingo_removed_reactions == 0 )
        
        # perform an FBA to check the result after reactions removal
        res = dingo_model.fba()
        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)       



        # load models in cobra and dingo format again to restore bounds
        cobra_model = load_json_model("ext_data/e_coli_core.json")        

        # call the reduce function from the PreProcess class 
        # with extend=1 to remove additional reactions from the model        
        obj = PreProcess(cobra_model)        
        removed_reactions, dingo_model = obj.reduce(extend=1)        
    
        # calculate the count of removed reactions with extend set to 1        
        removed_reactions_count = len(removed_reactions)
        self.assertTrue( 47 - removed_reactions_count == 0 )

        # calculate the count of reactions with bounds equal to 0 
        # with extend set to 1 from the dingo model
        dingo_removed_reactions = np.sum((dingo_model.lb == 0) & (dingo_model.ub == 0))
        self.assertTrue( 47 - dingo_removed_reactions == 0 )
        
        # perform an FBA to check the result after reactions removal
        res = dingo_model.fba()
        self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)       


if __name__ == "__main__":
    unittest.main()

 