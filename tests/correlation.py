
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
        
        
        # calculate correlation matrix without filtering from copula indicator
        corr_matrix = correlated_reactions(steady_states, 
                                           indicator_cutoff=0,
                                           pearson_cutoff = 0,
                                           lower_triangle = False)

        plot_corr_matrix(corr_matrix, reactions, format="svg")
        
        
        
        # calculate correlation matrix without filtering from copula indicator
        corr_matrix = correlated_reactions(steady_states, 
                                           indicator_cutoff=2,
                                           pearson_cutoff = 0.5,
                                           lower_triangle = False)

        plot_corr_matrix(corr_matrix, reactions, format="svg")

        # sum values in the diagonal of the correlation matrix ==> 95*pearson ==> 95*1
        self.assertTrue(np.trace(corr_matrix) == len(reactions))
        # rows and columns must be equal to model reactions
        self.assertTrue(corr_matrix.shape[0] == len(reactions))
        self.assertTrue(corr_matrix.shape[1] == len(reactions))



        cobra_model = load_json_model('ext_data/e_coli_core.json')          
          
        # reduce the metabolic model
        obj = PreProcess(cobra_model, tol=1e-5, open_exchanges=False)  
        removed_reactions, reduced_dingo_model = obj.reduce(extend=False)      
        reactions = reduced_dingo_model.reactions

        for reaction in reactions:
            index = reactions.index(reaction)
            if reaction in removed_reactions:
                reactions[index] = None

        sampler = PolytopeSampler(reduced_dingo_model)
        steady_states = sampler.generate_steady_states()

        # calculate correlation matrix with additional filtering from copula indicator       
        filtered_corr_matrix = correlated_reactions(
                                steady_states,  
                                pearson_cutoff = 0,
                                indicator_cutoff = 0, 
                                cells = 10,
                                cop_coeff = 0.3,
                                lower_triangle = False)
                
        plot_corr_matrix(filtered_corr_matrix, reactions, format="svg")
        

if __name__ == "__main__":
    unittest.main()
