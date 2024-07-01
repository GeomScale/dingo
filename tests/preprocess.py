
import unittest
import os
from dingo.preprocess import PreProcess

class TestPreprocess(unittest.TestCase):

    def test_preprocess(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        PreProcess.blocked(input_file_json)
        #model = MetabolicNetwork.from_json(input_file_json)

        #self.assertTrue()
        
        

if __name__ == "__main__":
    unittest.main()

 