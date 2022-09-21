import os, sys

import scipy
import numpy as np
import time

from dingo.gurobi_based_implementations import  fast_fba
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.scaling import gmscale
from dingo.utils import apply_scaling
from dingo.utils import remove_almost_redundant_facets

from volestipy import HPolytope

from PolyRound.api import PolyRoundApi
from PolyRound.static_classes.lp_utils import ChebyshevFinder
from PolyRound.settings import PolyRoundSettings

import hopsy
import matplotlib.pyplot as plt

import pickle
import copy 

def is_in(point: np.ndarray, polytope: HPolytope, loop: bool = True):
    if loop:
        for i in range(polytope.A().shape[0]):
            h = polytope.A()[i, :]
            if h.dot(point) > polytope.b()[i]:
                return False
        return True
    else:
        return (polytope.A().dot(point) <= polytope.b()).all()

# Create a settings object with the default settings.
settings = PolyRoundSettings()

#models = ["e_coli_core","iLJ478","iSB619","iHN637","iJN678","iNF517","iJN746","iAB_RBC_283","iJR904","iAT_PLT_636","iSDY_1059","iAF1260","iEC1344_C","iJO1366","iBWG_1329","iML1515","RECON1","Recon2","Recon3D_latest"]
# iNF517 segfaults
#models = ["e_coli_core","iLJ478","iSB619","iHN637","iJN678","iJN746","iAB_RBC_283"]
model = str(sys.argv[1])
with_polyround = str(sys.argv[2])

#for model in models:

## Import model and create Polytope object

current_directory = os.getcwd()
input_file = current_directory + "/ext_data/" + model  + ".xml"

polytope = PolyRoundApi.sbml_to_polytope(input_file)

## Simplify and tranform with PolyRound

print('Simplify')

simpl_trnsf_t1 = time.process_time()

# Remove redundant constraints and refunction inequality constraints that are de-facto equalities.
# Due to these inequalities, the polytope is empty (distance from chebyshev center to boundary is zero)
#x, dist = ChebyshevFinder.chebyshev_center(polytope, settings)
#print(dist)
simplified_polytope = PolyRoundApi.simplify_polytope(polytope)
# The simplified polytope has non-zero border distance
#x, dist = ChebyshevFinder.chebyshev_center(simplified_polytope, settings)
#print(dist)
transformed_polytope = PolyRoundApi.transform_polytope(simplified_polytope)
# The distance from the chebyshev center to the boundary changes in the new coordinate system
#x, dist = ChebyshevFinder.chebyshev_center(transformed_polytope, settings)
#print(dist)
#print(transformed_polytope.A.shape)
#print("Simplify + transform time= ", time.process_time() - simpl_trnsf_t1)

simpl_trnsf_t2 = time.process_time()

with open("simplified_transformed_polytopes/" + model + ".pkl",'wb') as f:
    pickle.dump(transformed_polytope, f)

## Sampling with dingo

sampling_dingo_t1 = time.process_time()

sampling_dingo_t1=0
sampling_dingo_t1=0
transforming_dingo_t1=0
transforming_dingo_t2=0
num_of_samples=10

print('Sampling')

if (with_polyround == "1"):
    sampling_dingo_t1 = time.process_time()
    
    A = copy.deepcopy(transformed_polytope.A.to_numpy())
    b = copy.deepcopy(transformed_polytope.b.to_numpy())
    
    res = gmscale(A, 0.99)
    res = apply_scaling(A, b, res[0], res[1])
    As = res[0]
    bs = res[1]
    
    res = remove_almost_redundant_facets(As, bs)
    Ass = res[0]
    bss = res[1]
    samples_transformed, dingo_psrf = PolytopeSampler.sample_from_polytope(Ass, bss,
	                                                                       num_of_samples, True, False, 1)
    with open("dingo_samples_polyround_transformation/" + model + ".pkl",'wb') as f:
        pickle.dump(samples_transformed, f)
    sampling_dingo_t2 = time.process_time()
else:
    transforming_dingo_t1 = time.process_time()

    input_file_json = current_directory + "/ext_data/" + model  + ".json"

    mn_model = MetabolicNetwork.from_json(input_file_json)
    mn_model.set_fast_mode()

    sampler = PolytopeSampler(mn_model)
    sampler.set_fast_mode()

    dingo_transformed_polytope = sampler.get_polytope()

    transforming_dingo_t2 = time.process_time()

    with open("dingo_simplified_transformed_polytopes/" + model + ".pkl",'wb') as f:
        pickle.dump(dingo_transformed_polytope, f)

    sampling_dingo_t1 = time.process_time()
    
    A = copy.deepcopy(sampler._A)
    b = copy.deepcopy(sampler._b)
    
    samples_transformed, dingo_psrf = PolytopeSampler.sample_from_polytope(A,
                                                                           b,
	                                                                       num_of_samples, True, False, 1)
    sampling_dingo_t2 = time.process_time()

    with open("dingo_samples/" + model + ".pkl",'wb') as f:
        pickle.dump(samples_transformed, f)

#samples_transformed = np.transpose(samples_transformed)
#P = HPolytope(transformed_polytope.A.to_numpy(),transformed_polytope.b.to_numpy())
#print(is_in(samples_transformed[2],P))
#print("Sampling time= ", time.process_time() - sample_t)

sampling_dingo_t2 = time.process_time()


## Rounding with PolyRound

rounding_t1 = time.process_time()
rounded_polytope = PolyRoundApi.round_polytope(transformed_polytope)
rounding_t2 = time.process_time()


## Sampling with Hopsy

sampling_hopsy_t1 = time.process_time()

problem = hopsy.Problem(rounded_polytope.A.to_numpy(), rounded_polytope.b.to_numpy())

proposal = hopsy.UniformCoordinateHitAndRunProposal

starting_point = hopsy.compute_chebyshev_center(problem)

#chain = hopsy.MarkovChain(problem, starting_point=starting_point)
#rng = hopsy.RandomNumberGenerator(seed=42)

number_of_chains = 1
markov_chains = [hopsy.MarkovChain(problem, proposal, starting_point=starting_point) for i in range(number_of_chains)]
rng = [hopsy.RandomNumberGenerator(seed=i*42) for i in range(number_of_chains)]
#for chain in markov_chains:
#    chain.proposal.stepsize = 0.2

d = rounded_polytope.A.shape[1]
accrate, states = hopsy.sample(markov_chains, rng, n_samples=num_of_samples, thinning=100*d)
sampling_hopsy_t2 = time.process_time()

rhat = hopsy.rhat(states)
ess = hopsy.ess(states)

print(model, d, rounded_polytope.A.shape[0],
    transforming_dingo_t2 - transforming_dingo_t1,
    sampling_dingo_t2 - sampling_dingo_t1,
    dingo_psrf,
    simpl_trnsf_t2 - simpl_trnsf_t1,
    rounding_t2 - rounding_t1,
    sampling_hopsy_t2 - sampling_hopsy_t1,
    rhat.max(),  ess.min())

sampling_hopsy_tB1 = time.process_time()
accrate2, states2 = hopsy.sample(markov_chains, rng, n_samples=num_of_samples, thinning=8*d*d)
sampling_hopsy_tB2 = time.process_time()

# when doing mcmc, assessing convergence is very important!
#rhat = hopsy.rhat(states, series=10)
#print('Acceptance rates (per chain):', *accrate)

# checks convergence and ESS
rhat2 = hopsy.rhat(states2)
#print(rhat.max())
#assert((rhat<=1.1).all()) # asserts that convergence has been reached. Here we use a strict limit of 1.01

ess2 = hopsy.ess(states2)
#print(ess.min())
#assert((ess >= 400).all()) # asserts that we have reached the minimum number of effective samples we want, in this case 400.

print(sampling_hopsy_tB2 - sampling_hopsy_tB1,
    rhat2.max(), ess2.min())

