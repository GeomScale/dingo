#!/usr/bin/python3
import os, sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo import plot_histogram
from time import process_time
import pickle
import multiprocessing

from volestipy import HPolytope


def is_in(point: np.ndarray, polytope: HPolytope, loop: bool = True):
   if loop:
      for i in range(polytope.A.shape[0]):
         h = polytope.A.iloc[[i]].values
         product = h.dot(point)
         if h.dot(point)[0] > polytope.b[i]:
            return False
      loop = False
      return True
   else:
      return (polytope.A.dot(point) <= polytope.b).all()


def sample_on_polyround_processed_polytope(network, ess_value, psrf_status):

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
                                                        ess = int(ess_value),
                                                        psrf = psrf_status,
                                                        parallel_mmcs = False
                                                        )
   end = process_time() 
   sampling_runtime = end - start

   print("Dingo for the " + name + " model took " + str(sampling_runtime) + " sec to sample using simplified and transformed polytope with PolyRound.\n\n")

   with open("dingo_samples_on_polyrounded_polytopes/dingo_polyround" + name + ".pckl", "wb") as dingo_steadystates_file: 
         pickle.dump(steady_states, dingo_steadystates_file)

   return steady_states, polytope


if __name__ == '__main__':

   file_name = sys.argv[1]
   current_directory = os.getcwd()
   path_to_net = current_directory + "/" + file_name
   ess_asked = sys.argv[2]
   psrf_asked = sys.argv[3]
   dingo_samples, polytope_sampled = sample_on_polyround_processed_polytope(path_to_net, ess_asked, psrf_asked)

   dingo_samples_t = dingo_samples.T

   check = True
   for sample in dingo_samples_t: 
      control = is_in(sample, polytope_sampled)
      if control:
         continue
      else:
         print("huston we have a problem.")
         print(control)
         check = False
   if check:
      print("all samples in.")

