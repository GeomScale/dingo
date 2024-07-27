
from dingo.illustrations import plot_corr_matrix
from dingo.utils import correlated_reactions
from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np
from dingo import plot_copula


model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
reactions = model.reactions

sampler = PolytopeSampler(model)
steady_states = sampler.generate_steady_states()

corr_matrix, updated_corr_matrix = correlated_reactions(steady_states,  
                                                        pearson_cutoff = 0.80,
                                                        indicator_cutoff = 9e1000, 
                                                        n = 11)


arr = corr_matrix == updated_corr_matrix
rows, cols = np.where(arr == True)
for i in range(len(rows)):
    row = rows[i] 
    col = cols[i]
    print(reactions[row], reactions[col])
    #data_flux2=[steady_states[row],reactions[row]]
    #data_flux1=[steady_states[col],reactions[col]]

    #plot_copula(data_flux1, data_flux2, n=10)


plot_corr_matrix(corr_matrix, reactions, color="RdYlBu")
plot_corr_matrix(updated_corr_matrix, reactions, color="RdYlBu")