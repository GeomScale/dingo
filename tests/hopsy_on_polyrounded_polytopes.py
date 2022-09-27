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
for i in range(number_of_chains):

    if counter == 0:
        markov_chain = hopsy.MarkovChain(problem, proposal, starting_point = starting_point)
    else:
        markov_chain = hopsy.MarkovChain(problem, proposal, starting_point = last_point_of_previous_chain)
    
    rng = hopsy.RandomNumberGenerator(seed = i*42) 
    accrate, states = hopsy.sample(markov_chain, rng, n_samples = num_of_samples, thinning = thinning_value)

    if counter == 0:
        total_samples = states.copy()
    else: 
        total_samples = np.vstack([total_samples, states])

    last_point_of_previous_chain = states[0][-1,:]
    counter += 1


sampling_hopsy_t2 = time.process_time()
rhat = hopsy.rhat(total_samples)
ess = hopsy.ess(total_samples)

print("model: ", polyrounded_polytope_file, "\n")
print("polyround_A.shape[0] : ", polyround_A.shape[0], "\n")
print("d (polyround_A.shape[1]): ", polyround_A.shape[1], "\n")
print("hopsy sampling time: ", str(sampling_hopsy_t2 - sampling_hopsy_t1), "\n")
print("rhat.max : ", rhat.max(), "\n")
print("ess.min: ", ess.min(), "\n")

