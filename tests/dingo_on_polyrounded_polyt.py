#!/usr/bin/python3
import os, sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo import plot_histogram
from time import process_time
import pickle
import multiprocessing

def sample_on_polyround_processed_polytope(network):

   current_directory = os.getcwd()
   name              = network.split("/")[-1]
   polytope_file     = network

   with open(polytope_file, "rb") as f:
      obj = pickle.load(f)
   polytope = obj[0]

   polyround_A = polytope.A.to_numpy()
   polyround_b = polytope.b.to_numpy()


   # Sample from the polytope built using the parrallel MMCS 
   start = process_time() 
   steady_states = PolytopeSampler.sample_from_polytope(polyround_A,
                                                        polyround_b,
                                                        ess = 1000,
                                                        psrf = True,
                                                        parallel_mmcs = False
                                                        )
   end = process_time() 
   sampling_runtime = end - start

   print("Dingo for the " + name + " model took " + str(sampling_runtime) + " sec to sample using polyround processed polytope")

   with open("dingo_samples_on_polyrounded_polytopes/dingo_polyround" + name + ".pckl", "wb") as dingo_steadystates_file: 
         pickle.dump(steady_states, dingo_steadystates_file)

if __name__ == '__main__':

   file_name = sys.argv[1]
   current_directory = os.getcwd()
   path_to_net = current_directory + "/polyrounded_polytopes/" + file_name
   sample_on_polyround_processed_polytope(path_to_net)



