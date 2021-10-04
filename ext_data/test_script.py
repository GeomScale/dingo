import unittest
import os
import scipy
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo import plot_histogram
#from dingo.gurobi_based_implementations import fast_inner_ball


current_directory = os.getcwd()
input_file_json = current_directory + "/ext_data/iSB619.json"

model = MetabolicNetwork.from_json(input_file_json)
model.set_fast_mode()
#model.biomass_function = np.zeros(model.S.shape[1])

sampler = PolytopeSampler(model)
sampler.set_fast_mode()

steady_states = sampler.generate_steady_states(ess = 3000, psrf=True, parallel_mmcs = False, num_threads = 1)

reactions = model.reactions
plot_histogram(
        steady_states[13],
        reactions[13],
        n_bins = 60,
        )
