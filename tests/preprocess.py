
import unittest
import os
from dingo.preprocess import PreProcess

class TestFba(unittest.TestCase):

    def test_fba_json(self):

        #input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        PreProcess.blocked(input_file_json)
        #model = MetabolicNetwork.from_json(input_file_json)
        #model.set_slow_mode()
        #res = model.fba()

        #self.assertTrue(abs(res[1] - 0.8739215067486387) < 1e-03)
        
        

if __name__ == "__main__":
    unittest.main()

 