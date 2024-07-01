
from cobra.io import load_json_model
from dingo.preprocess import PreProcess
import unittest
import os


class TestPreprocess(unittest.TestCase):

    def test_preprocess(self):

        model = load_json_model("ext_data/e_coli_core.json")
        obj = PreProcess(model)
        
        essentials = len(obj.essential_reactions)

        self.assertTrue( 27-essentials < 0.01)
        
        

if __name__ == "__main__":
    unittest.main()

 