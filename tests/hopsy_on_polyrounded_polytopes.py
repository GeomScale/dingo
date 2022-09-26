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

number_of_chains = 1
rng = [hopsy.RandomNumberGenerator(seed=i*42) for i in range(number_of_chains)]
starting_point = hopsy.compute_chebyshev_center(problem)

# ------- Uniform distribution
proposal = hopsy.UniformCoordinateHitAndRunProposal
markov_chains = [hopsy.MarkovChain(problem, proposal, starting_point=starting_point) for i in range(number_of_chains)]


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

accrate, states = hopsy.sample(markov_chains, rng, n_samples=num_of_samples, thinning=thinning_value)

sampling_hopsy_t2 = time.process_time()

rhat = hopsy.rhat(states)
ess = hopsy.ess(states)

print("model: ", polyrounded_polytope_file, "\n")
print("polyround_A.shape[0] : ", polyround_A.shape[0], "\n")
print("d (polyround_A.shape[1]): ", polyround_A.shape[1], "\n")
print("hopsy sampling time: ", str(sampling_hopsy_t2 - sampling_hopsy_t1), "\n")
print("rhat.max : ", rhat.max(), "\n")
print("ess.min: ", ess.min(), "\n")



# sampling_hopsy_tB1 = time.process_time()
# accrate2, states2 = hopsy.sample(markov_chains, rng, n_samples=num_of_samples, thinning=8*d*d)
# sampling_hopsy_tB2 = time.process_time()

# # when doing mcmc, assessing convergence is very important!
# #rhat = hopsy.rhat(states, series=10)
# #print('Acceptance rates (per chain):', *accrate)

# # checks convergence and ESS
# rhat2 = hopsy.rhat(states2)
# #print(rhat.max())
# #assert((rhat<=1.1).all()) # asserts that convergence has been reached. Here we use a strict limit of 1.01

# ess2 = hopsy.ess(states2)
# #print(ess.min())
# #assert((ess >= 400).all()) # asserts that we have reached the minimum number of effective samples we want, in this case 400.

# print(sampling_hopsy_tB2 - sampling_hopsy_tB1,
#     rhat2.max(), ess2.min())
