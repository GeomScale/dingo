
from dingo.illustrations import plot_corr_matrix
from dingo.utils import correlated_reactions
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.preprocess import PreProcess
from cobra.io import load_json_model
import numpy as np
import unittest

class TestCorrelation(unittest.TestCase):

    def test_correlation(self):
        dingo_model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
        reactions = dingo_model.reactions

        sampler = PolytopeSampler(dingo_model)
        steady_states = sampler.generate_steady_states()

        corr_matrix = correlated_reactions(steady_states, indicator_cutoff=0)

        # sum values in the diagonal of the correlation matrix ==> 95*pearson ==> 95*1
        self.assertTrue(np.trace(corr_matrix) == len(reactions))
        # rows and columns must be equal to model reactions
        self.assertTrue(corr_matrix.shape[0] == len(reactions))
        self.assertTrue(corr_matrix.shape[1] == len(reactions))

        plot_corr_matrix(corr_matrix, reactions)



        cobra_model = load_json_model('ext_data/e_coli_core.json')          
        obj = PreProcess(cobra_model, tol=1e-5, open_exchanges=False)  
        removed_reactions, reduced_dingo_model = obj.reduce(extend=True)      
        reactions = reduced_dingo_model.reactions

        sampler = PolytopeSampler(reduced_dingo_model)
        steady_states = sampler.generate_steady_states()

        corr_matrix, filtered_corr_matrix = correlated_reactions(
                                                steady_states,  
                                                pearson_cutoff = 0.6,
                                                indicator_cutoff = 2, 
                                                cells = 10,
                                                cop_coeff = 0.3)
        
        # sum values in the diagonal of the correlation matrix ==> 95*pearson ==> 95*1                                        
        self.assertTrue(np.trace(corr_matrix) == len(reactions))
        # rows and columns must be equal to model reactions        
        self.assertTrue(filtered_corr_matrix.shape[0] == len(reactions))
        self.assertTrue(filtered_corr_matrix.shape[1] == len(reactions))

        plot_corr_matrix(corr_matrix, reactions)
        plot_corr_matrix(filtered_corr_matrix, reactions)


if __name__ == "__main__":
    unittest.main()
