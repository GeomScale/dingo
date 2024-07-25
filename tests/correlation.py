
from dingo.illustrations import plot_corr_matrix
from dingo.utils import correlated_reactions
from dingo import MetabolicNetwork, PolytopeSampler


model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
reactions = model.reactions

sampler = PolytopeSampler(model)
steady_states = sampler.generate_steady_states()

plot_corr_matrix(steady_states, reactions, color="RdYlBu")
correlated_reactions(steady_states, reactions, pearson_cutoff = 0.7, n = 10)


