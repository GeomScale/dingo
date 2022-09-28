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
name = sys.argv[1].split("/")[-1]
polyrounded_polytope_file = os.getcwd() + "/" + sys.argv[1]
with open(polyrounded_polytope_file, "rb") as f:
    obj = pickle.load(f)
polytope = obj[0]


## Sampling with Hopsy
polyround_A = polytope.A.to_numpy()
polyround_b = polytope.b.to_numpy()

sampling_hopsy_t1 = time.process_time()

problem = hopsy.Problem(polyround_A, polyround_b)

number_of_chains = 4
# rng = [hopsy.RandomNumberGenerator(seed=i*42) for i in range(number_of_chains)]
starting_point = hopsy.compute_chebyshev_center(problem)

# ------- Uniform distribution
proposal = hopsy.UniformCoordinateHitAndRunProposal
# markov_chains = [hopsy.MarkovChain(problem, proposal, starting_point = starting_point) \
#                 for i in range(number_of_chains)
#     ]


# # ----- Gaussian 
# proposal = hopsy.GaussianCoordinateHitAndRunProposal
# markov_chains = [hopsy.MarkovChain(problem, proposal, starting_point=starting_point) for i in range(number_of_chains)]

# tuning_target = hopsy.AcceptanceRateTarget(markov_chains)
# TS = hopsy.ThompsonSamplingTuning()
# stepsize, tuning_posterior = hopsy.tune(method=TS, target=tuning_target, rngs=rng)


# # sets stepsize for every markov chain
# for mc in markov_chains:
#     mc.proposal.stepsize = stepsize

# # ------  end of Gaussian


thinning_value = polyround_A.shape[1]*100
# accrate, states = hopsy.sample(markov_chains, rng, n_samples=num_of_samples, thinning=thinning_value)

counter = 0
total_time = 0
ess_check = True
n_samples_per_step = 10000

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
            
            timecost_per_point = float(total_time / n_samples_per_step)

            total_samples = states.copy()
            chains_on_the_run = np.array(np.split(states[0], 5))
            ess_local         = hopsy.ess(chains_on_the_run) # needs to be a 3-d array

        else:
            total_samples = np.concatenate((total_samples, states), axis = 0) 

            tmp               = np.concatenate((total_samples[0], states[0]), axis = 0)
            chains_on_the_run = np.array(np.split(tmp, 5))      
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



print("model: ", polyrounded_polytope_file)
print("polyround_A.shape[0] : ", polyround_A.shape[0])
print("d (polyround_A.shape[1]): ", polyround_A.shape[1])
print("time cost per sampling point: ", str(timecost_per_point))
print("\n~~~~~\n")
print("hopsy sampling time: ", str(total_time))
print("\n~~~~~\n")
print("rhat.max : ", rhat.max())
print("ess.min: ", ess.min())




