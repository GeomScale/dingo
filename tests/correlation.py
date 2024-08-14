
from dingo.utils import correlated_reactions
from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np
import unittest

class TestCorrelation(unittest.TestCase):

    def test_correlation(self):
        
        dingo_model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
        reactions = dingo_model.reactions

        sampler = PolytopeSampler(dingo_model)
        steady_states = sampler.generate_steady_states()
        
        # calculate correlation matrix with filtering from copula indicator
        corr_matrix, indicator_dict = correlated_reactions(steady_states,
                                        reactions = reactions,
                                        indicator_cutoff = 5,
                                        pearson_cutoff = 0.999999,
                                        lower_triangle = False,
                                        verbose = False)

        # sum values in the diagonal of the correlation matrix ==> 95*pearson ==> 95*1
        self.assertTrue(np.trace(corr_matrix) == len(reactions))
        # rows and columns must be equal to model reactions
        self.assertTrue(corr_matrix.shape[0] == len(reactions))
        self.assertTrue(corr_matrix.shape[1] == len(reactions))
                
        
        dingo_model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
        reactions = dingo_model.reactions

        sampler = PolytopeSampler(dingo_model)
        steady_states = sampler.generate_steady_states()
        
        # calculate correlation matrix without filtering from copula indicator
        corr_matrix = correlated_reactions(steady_states, 
                                           indicator_cutoff = 0,
                                           pearson_cutoff = 0,
                                           lower_triangle = True,
                                           verbose = False)

        # sum values in the diagonal of the correlation matrix ==> 95*pearson ==> 95*1
        self.assertTrue(np.trace(corr_matrix) == len(reactions))
        # rows and columns must be equal to model reactions
        self.assertTrue(corr_matrix.shape[0] == len(reactions))
        self.assertTrue(corr_matrix.shape[1] == len(reactions))
   

if __name__ == "__main__":
    unittest.main()
