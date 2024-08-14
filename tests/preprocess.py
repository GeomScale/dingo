
from cobra.io import load_json_model
from dingo import MetabolicNetwork
from dingo.preprocess import PreProcess
import unittest
import numpy as np

class TestPreprocess(unittest.TestCase):

    def test_preprocess(self):

        # load cobra model
        cobra_model = load_json_model("ext_data/e_coli_core.json")
        
        # convert cobra to dingo model
        initial_dingo_model = MetabolicNetwork.from_cobra_model(cobra_model)
        
        # perform an FBA to find the initial FBA solution
        initial_fba_solution = initial_dingo_model.fba()[1]


        # call the reduce function from the PreProcess class 
        # with extend=False to remove reactions from the model        
        obj = PreProcess(cobra_model, tol=1e-5, open_exchanges=False, verbose=False)  
        removed_reactions, final_dingo_model = obj.reduce(extend=False)      
        
        # calculate the count of removed reactions with extend set to False        
        removed_reactions_count = len(removed_reactions)
        self.assertTrue( 46 - removed_reactions_count == 0 )

        # calculate the count of reactions with bounds equal to 0 
        # with extend set to False from the dingo model
        dingo_removed_reactions = np.sum((final_dingo_model.lb == 0) & (final_dingo_model.ub == 0))
        self.assertTrue( 46 - dingo_removed_reactions == 0 )
        
        # perform an FBA to check the solution after reactions removal
        final_fba_solution = final_dingo_model.fba()[1]
        self.assertTrue(abs(final_fba_solution - initial_fba_solution) < 1e-03)       


        # load models in cobra and dingo format again to restore bounds
        cobra_model = load_json_model("ext_data/e_coli_core.json")
        
        # convert cobra to dingo model
        initial_dingo_model = MetabolicNetwork.from_cobra_model(cobra_model)
     
        # call the reduce function from the PreProcess class 
        # with extend=True to remove additional reactions from the model        
        obj = PreProcess(cobra_model, tol=1e-6, open_exchanges=False, verbose=False)        
        removed_reactions, final_dingo_model = obj.reduce(extend=True)        
    
        # calculate the count of removed reactions with extend set to True        
        removed_reactions_count = len(removed_reactions)
        self.assertTrue( 46 - removed_reactions_count <= 0 )

        # calculate the count of reactions with bounds equal to 0 
        # with extend set to True from the dingo model
        dingo_removed_reactions = np.sum((final_dingo_model.lb == 0) & (final_dingo_model.ub == 0))
        self.assertTrue( 46 - dingo_removed_reactions <= 0 )
        
        # perform an FBA to check the result after reactions removal
        final_fba_solution = final_dingo_model.fba()[1]
        self.assertTrue(abs(final_fba_solution - initial_fba_solution) < 1e-03)   


if __name__ == "__main__":
    unittest.main()

 