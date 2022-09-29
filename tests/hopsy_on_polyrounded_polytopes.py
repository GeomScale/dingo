#!/usr/bin/python3

import os, sys

import scipy
import numpy as np
import time

import hopsy
import matplotlib.pyplot as plt

import pickle
import copy 


# General
num_of_samples = 10

## Import polyrounded polytope
name = sys.argv[1].split("/")[-1].replace(".xml.pckl", "")
polyrounded_polytope_file = os.getcwd() + "/" + sys.argv[1]
with open(polyrounded_polytope_file, "rb") as f:
    obj = pickle.load(f)
polytope = obj[0]


## Sampling with Hopsy
polyround_A = polytope.A.to_numpy()
polyround_b = polytope.b.to_numpy()

problem        = hopsy.Problem(polyround_A, polyround_b)
starting_point = hopsy.compute_chebyshev_center(problem)
proposal       = hopsy.UniformCoordinateHitAndRunProposal
thinning_value = polyround_A.shape[1]*100

counter = 0
total_time = 0
ess_check = True
n_samples_per_step = 10



# for i in range(3):

#     print("\n\n ~~~ \n\nITERATION #", str(i))

#     if counter == 0:
#         markov_chain = hopsy.MarkovChain(problem, proposal, starting_point = starting_point)
#     else:
#         markov_chain = hopsy.MarkovChain(problem, proposal, starting_point = last_point_of_previous_chain)
    
#     rng = hopsy.RandomNumberGenerator(seed = 42) 
#     t_0 = time.process_time()

#     accrate, states = hopsy.sample(markov_chain, rng, n_samples = n_samples_per_step, thinning = thinning_value)

#     t_1         = time.process_time()
#     total_time += t_1 - t_0

#     if counter == 0:
        
#         total_samples     = states        # each sample is added as an extra chain, e.g (1, n_samples_per_step, d) --> (2, n_samples_per_step, d) 
#         unified_chain     = states[0]     # put all samples in a 2d array ((i+1)*n_samples_per_step , d)

#     else:
#         total_samples     = np.concatenate((total_samples, states), axis = 0)
#         unified_chain     = np.concatenate((unified_chain, states[0]), axis = 0)

#     chains_on_the_run = np.array(np.split(unified_chain, 5))  # split the unified_chain in 5 chains 

#     j = i*n_samples_per_step

#     print("\nSamples as coming from hopsy iterations: \n")
#     print(total_samples[:,:,:3])
    
    
#     print("\nThe new points added in the unified chain: \n")
#     print(unified_chain[-num_of_samples:,:3])

#     # print("\nThe first 2 points in the 5 chains: \n", chains_on_the_run[:,:2,:3])
#     print("\n\npoints in the 5 chains: \n", chains_on_the_run[:,:,:3])


#     print(unified_chain.shape, chains_on_the_run.shape)

#     last_point_of_previous_chain = states[0][-1,:]
#     counter += 1




sys.exit(0)



while ess_check: 

    if counter == 0:
        markov_chain = hopsy.MarkovChain(problem, proposal, starting_point = starting_point)
    else:
        markov_chain = hopsy.MarkovChain(problem, proposal, starting_point = last_point_of_previous_chain)
    
    rng = hopsy.RandomNumberGenerator(seed = 42) 
    t_0 = time.process_time()

    accrate, states = hopsy.sample(markov_chain, rng, n_samples = n_samples_per_step, thinning = thinning_value)

    t_1         = time.process_time()
    total_time += t_1 - t_0

    if counter == 0:
        
        total_samples     = states        # each sample is added as an extra chain, e.g (1, n_samples_per_step, d) --> (2, n_samples_per_step, d) 
        unified_chain     = states[0]     # put all samples in a 2d array ((i+1)*n_samples_per_step , d)

    else:
        total_samples     = np.concatenate((total_samples, states), axis = 0)
        unified_chain     = np.concatenate((unified_chain, states[0]), axis = 0)

    chains_on_the_run = np.array(np.split(unified_chain, 5))  # split the unified_chain in 5 chains 

    print("hello friend")

    ess_local         = hopsy.ess(chains_on_the_run)

    if ess_local.min() > 1000: 

        final_chains = chains_on_the_run
        ess_check    = False

    last_point_of_previous_chain = states[0][-1,:]
    
    counter += 1


rhat = hopsy.rhat(final_chains)
ess = hopsy.ess(final_chains)

with open("hopsy_samples/total_samples_" + name + ".pckl", "wb") as hopsy_samples_file: 
        pickle.dump(total_samples, hopsy_samples_file)


with open("hopsy_samples/total_samples_" + name + ".txt", "w") as stats:
    stats.write("model: " + polyrounded_polytope_file + "\n")
    stats.write("polyround_A.shape[0] : " + str(polyround_A.shape[0])+ "\n")
    stats.write("d (polyround_A.shape[1]): " + str(polyround_A.shape[1])+ "\n")
    stats.write("\n~~~~~\n")
    stats.write("hopsy sampling time: " + str(total_time))
    stats.write("\n~~~~~\n")
    stats.write("rhat.max : " + str(rhat.max())+ "\n")
    stats.write("ess.min: " + str(ess.min())+ "\n")




