
from cobra.io import load_json_model
from dingo.preprocess import PreProcess
import unittest
import os


class TestPreprocess(unittest.TestCase):

    def test_preprocess(self):

        model = load_json_model("ext_data/e_coli_core.json")
        obj = PreProcess(model)

        removed_reactions_count = len(obj.removed_reactions_ids())
        
        self.assertTrue( 55 - removed_reactions_count == 0 )
        

if __name__ == "__main__":
    unittest.main()

 