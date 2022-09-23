#!/usr/bin/python3

# author: H.Z. 
# date: 2022-09-23


import os, sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo import plot_histogram
import time
import pickle
import multiprocessing


def dingo_preprocess_and_sampling(network):

   name = network
   current_directory = os.getcwd()
   model_stats_file = open(name + "_stats.log", "w")
   input_file_json = current_directory + "/ext_data/"
   input_file_json = input_file_json + network

   model = MetabolicNetwork.from_json(input_file_json)
   model.set_fast_mode()
   model.biomass_function = np.zeros(model.S.shape[1]) # comment out this line to maximize the biomass objective

   #-------------------------------remove redundant facets----------------------------------------#
   try:

      sampler = PolytopeSampler(model)
      sampler.set_fast_mode()
      print("PolytopeSampler instance for model " + name + " has initiated and get_polytope() is about to start.")

      start = time.time()
      A, b, N, N_shift = sampler.get_polytope()
      end = time.time()

      print("A polytope is now ready for model " + name + " and its shape is: " ) ; print(A.shape)

      prepro_runtime = end - start

      print("Model " + name + " took " + str(prepro_runtime) + " sec to preprocess having redundant facets removed")
      model_stats_file.write("Model " + name + " took " + str(prepro_runtime) + " sec to preprocess having redundant facets removed" + "\n")

      polytope_info = (
         sampler,
         name,
      )

      with open(
         "polytopes_redundancy_removals/polytope_" + name + ".pckl", "wb"
      ) as dingo_polytope_file:
         pickle.dump(polytope_info, dingo_polytope_file)

      start = time.time()
      print("A preprocessed polytope is now ready for sampling regarding the network" + name)
      samples = sampler.sample_from_polytope(A, b, ess=1000, psrf=True, parallel_mmcs=False, num_threads=4)
      end = time.time()

      sampling_runtime = end - start  ## You have to save runtime
      print("Model " + name + " took " + str(sampling_runtime) + " sec to sample having redundant facets removed")
      model_stats_file.write("Model " + name + " took " + str(sampling_runtime) + " sec to sample having redundant facets removed" + "\n")

      extra_2 = np.full((samples.shape[1], N.shape[0]), N_shift)
      steady_states = N.dot(samples) + extra_2.T
      with open("steady_states_redundancy_removals/steady_states_" + name + ".pckl", "wb") as dingo_steadystates_file: 
            pickle.dump(steady_states, dingo_steadystates_file)

   except ImportError as e:
      print("Removing reduntant facets for model " + model + " failed.")
      pass

   #--------------------------------DO NOT remove redundant facets-----------------------------------------#

   model = MetabolicNetwork.from_json(input_file_json)
   model.set_fast_mode()
   model.biomass_function = np.zeros(model.S.shape[1]) # comment out this line to maximize the biomass objective

   try:
      sampler = PolytopeSampler(model)
      sampler.set_fast_mode()

      # This was let as by default earlier 
      sampler.facet_redundancy_removal(False)

      print("PolytopeSampler instance for model " + name + " has initiated and get_polytope() is about to start. Redundant facets will not be removed this time.")

      start = time.time()
      A, b, N, N_shift = sampler.get_polytope()
      end = time.time()
      prepro_runtime = end - start  ## You have to save runtime

      print("Model " + name + " took " + str(prepro_runtime) + " sec to preprocess with redundant facets included")
      model_stats_file.write("Model " + name + " took " + str(prepro_runtime) + " sec to preprocess with redundant facets included" + "\n")

      polytope_info = (
            sampler,
            name,
      )

      with open(
            "polytopes_without_removals/polytope_" + name + ".pckl", "wb"
      ) as dingo_polytope_file:
            pickle.dump(polytope_info, dingo_polytope_file)

      start = time.time()

      try:
         samples_1 = sampler.sample_from_polytope(A, b, ess=1000, psrf=True, parallel_mmcs=False, num_threads=1)
      except:
          print("Model " + name + " was not sampled with redundant facets included.")
          return 1
   
      end = time.time()

      sampling_runtime = end - start  ## You have to save runtime
      print("Model " + name + " took " + str(sampling_runtime) + " sec to sample with redundant facets included")
      model_stats_file.write("Model " + name + " took " + str(sampling_runtime) + " sec to sample with redundant facets included" + "\n")

      extra_2 = np.full((samples_1.shape[1], N.shape[0]), N_shift)
      steady_states = N.dot(samples_1) + extra_2.T

      with open("steady_states_without_removals/steady_states_" + name + ".pckl", "wb") as dingo_steadystates_file:
            pickle.dump(steady_states, dingo_steadystates_file)

   except ImportError as e:
            pass

   print("Network " + name + " has been sampled. ^^^ \n\n")

if __name__ == '__main__':

   network = sys.argv[1]
   print("model ", network, " is about to be processed and sampled using dingo")
   dingo_preprocess_and_sampling(network)
   print("Network " + network + " completed.")




