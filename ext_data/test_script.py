import unittest
import os
import scipy
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
#from dingo.gurobi_based_implementations import fast_inner_ball


current_directory = os.getcwd()
input_file_json = current_directory + "/ext_data/e_coli_core.json"

model = MetabolicNetwork.from_json(input_file_json)
model.set_fast_mode()

sampler = PolytopeSampler(model)
sampler.set_fast_mode()

steady_states = sampler.generate_steady_states(ess = 1000, psrf=True, parallel_mmcs = False, num_threads = 1)
