# This is a cython wrapper for the C++ library volesti
# volesti (volume computation and sampling library)
  
# Copyright (c) 2012-2021 Vissarion Fisikopoulos
# Copyright (c) 2018-2021 Apostolos Chalkis
# Copyright (c) 2020-2021 Pedro Zuidberg Dos Martires
# Copyright (c) 2020-2021 Haris Zafeiropoulos
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

#!python
#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False

# Global dependencies
import os
import sys
import numpy as np
cimport numpy as np
from cpython cimport bool

# For the read the json format BIGG files function
import json
import scipy.io
# ----------------------------------------------------------------------------------

from dingo.pyoptinterface_based_impl import inner_ball

# Set the time
def get_time_seed():
   import random
   import time
   return int(time.time())


################################################################################
#                  Classes for the volesti C++ code                            #
################################################################################

# Get classes from the bindings.h file
cdef extern from "bindings.h":

   # The HPolytopeCPP class along with its functions
   cdef cppclass HPolytopeCPP:

      # Initialization
      HPolytopeCPP() except +
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables) except +

      # Compute volume
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);

      # Random sampling
      double apply_sampling(int walk_len, int number_of_points, int number_of_points_to_burn, \
                              int method, double* inner_point, double radius, double* samples, double variance_value, double* bias_vector)
      
      # Initialize the parameters for the (m)ultiphase (m)onte (c)arlo (s)ampling algorithm
      void mmcs_initialize(unsigned int d, int ess, int psrf_check, int parallelism, int num_threads);

      # Perform a step of (m)ultiphase (m)onte (c)arlo (s)ampling algorithm
      double mmcs_step(double* inner_point_for_c, double radius, int &N);

      # Get the samples and the transformation matrices from (m)ultiphase (m)onte (c)arlo (s)ampling algorithm
      void get_mmcs_samples(double* T_matrix, double* T_shift, double* samples);

      void get_polytope_as_matrices(double* new_A, double* new_b);

      # Rounding H-Polytope
      void apply_rounding(int rounding_method, double* new_A, double* new_b, double* T_matrix, \
                          double* shift, double &round_value, double* inner_point, double radius);

   # The lowDimPolytopeCPP class along with its functions
   cdef cppclass lowDimHPolytopeCPP:

      # Initialization
      lowDimHPolytopeCPP() except +
      lowDimHPolytopeCPP(double *A, double *b, double *Aeq, double *beq, int n_rows_of_A, int n_cols_of_A, int n_row_of_Aeq, int n_cols_of_Aeq) except +

      # Get full dimensional polytope
      int full_dimensiolal_polytope(double* N_extra_trans, double* shift, double* A_full_extra_trans, double* b_full)

# Lists with the methods supported by volesti for volume approximation and random walk
volume_methods = ["sequence_of_balls".encode("UTF-8"), "cooling_gaussian".encode("UTF-8"), "cooling_balls".encode("UTF-8")]
walk_methods = ["uniform_ball".encode("UTF-8"), "CDHR".encode("UTF-8"), "RDHR".encode("UTF-8"), "gaussian_ball".encode("UTF-8"), \
                "gaussian_CDHR".encode("UTF-8"), "gaussian_RDHR".encode("UTF-8"), "uniform_ball".encode("UTF-8"), "billiard".encode("UTF-8")]
rounding_methods = ["min_ellipsoid".encode("UTF-8"), "svd".encode("UTF-8"), "max_ellipsoid".encode("UTF-8")]

# Build the HPolytope class
cdef class HPolytope:

   cdef HPolytopeCPP polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b

   # Set the specs of the class
   def __cinit__(self, double[:,::1] A, double[::1] b):
      self._A = A
      self._b = b
      n_hyperplanes, n_variables = A.shape[0], A.shape[1]
      self.polytope_cpp = HPolytopeCPP(&A[0,0], &b[0], n_hyperplanes, n_variables)

   # This is where the volesti functions are getting their python interface; first the compute_volume() function
   def compute_volume(self, walk_len = 2, epsilon = 0.05, vol_method = "sequence_of_balls", walk_method = "uniform_ball", \
      np.npy_int32 seed=get_time_seed()):

      vol_method = vol_method.encode("UTF-8")
      walk_method = walk_method.encode("UTF-8")

      if vol_method in volume_methods:
         if walk_method in walk_methods:
            return self.polytope_cpp.compute_volume(vol_method, walk_method, walk_len, epsilon, seed)
         else:
            raise Exception('"{}" is not implemented to walk methods. Available methods are: {}'.format(walk_method, walk_methods))
      else:
         raise Exception('"{}" is not implemented to compute volume. Available methods are: {}'.format(vol_method, volume_methods))

   # Likewise, the generate_samples() function
   def generate_samples(self, method, number_of_points, number_of_points_to_burn, walk_len, variance_value, bias_vector, solver = None):

      n_variables = self._A.shape[1]
      cdef double[:,::1] samples = np.zeros((number_of_points, n_variables), dtype = np.float64, order = "C")

      # Get max inscribed ball for the initial polytope
      temp_center, radius = inner_ball(self._A, self._b, solver)
      
      cdef double[::1] inner_point_for_c = np.asarray(temp_center)
      
      cdef double[::1] bias_vector_ = np.asarray(bias_vector)
      
      if method == 'cdhr':
         int_method = 1
      elif method == 'rdhr':
         int_method = 2
      elif method == 'billiard_walk':
         int_method = 3
      elif method == 'ball_walk':
         int_method = 4
      elif method == 'dikin_walk':
         int_method = 5
      elif method == 'john_walk':
         int_method = 6
      elif method == 'vaidya_walk':
         int_method = 7
      elif method == 'gaussian_hmc_walk':
         int_method = 8
      elif method == 'exponential_hmc_walk':
         int_method = 9
      elif method == 'hmc_leapfrog_gaussian':
         int_method = 10
      elif method == 'hmc_leapfrog_exponential':
         int_method = 11
      else:
         raise RuntimeError("Uknown MCMC sampling method")
      
      self.polytope_cpp.apply_sampling(walk_len, number_of_points, number_of_points_to_burn, \
                                       int_method, &inner_point_for_c[0], radius, &samples[0,0], variance_value, &bias_vector_[0])
      return np.asarray(samples)


   # The rounding() function; as in compute_volume, more than one method is available for this step
   def rounding(self, rounding_method = 'john_position', solver = None):

      # Get the dimensions of the items about to build
      n_hyperplanes, n_variables = self._A.shape[0], self._A.shape[1]

      # Set the variables of those items; notice that they are all cdef type, except for the last one, which is used
      # both as a C++ and a Python variable
      cdef double[:,::1] new_A = np.zeros((n_hyperplanes, n_variables), dtype=np.float64, order="C")
      cdef double[::1] new_b = np.zeros(n_hyperplanes, dtype=np.float64, order="C")
      cdef double[:,::1] T_matrix = np.zeros((n_variables, n_variables), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n_variables), dtype=np.float64, order="C")
      cdef double round_value
      
      # Get max inscribed ball for the initial polytope
      center, radius = inner_ball(self._A, self._b, solver)
      
      cdef double[::1] inner_point_for_c = np.asarray(center)
      
      if rounding_method == 'john_position':
         int_method = 1
      elif rounding_method == 'isotropic_position':
         int_method = 2
      elif rounding_method == 'min_ellipsoid':
         int_method = 3
      else:
         raise RuntimeError("Uknown rounding method")

      self.polytope_cpp.apply_rounding(int_method, &new_A[0,0], &new_b[0], &T_matrix[0,0], &shift[0], round_value, &inner_point_for_c[0], radius)

      return np.asarray(new_A),np.asarray(new_b),np.asarray(T_matrix),np.asarray(shift),np.asarray(round_value)
   
   
   # (m)ultiphase (m)onte (c)arlo (s)ampling algorithm to generate steady states of a metabolic network
   def mmcs(self, ess = 1000, psrf_check = True, parallelism = False, num_threads = 2, solver = None):

      n_hyperplanes, n_variables = self._A.shape[0], self._A.shape[1]

      cdef double[:,::1] new_A = np.zeros((n_hyperplanes, n_variables), dtype=np.float64, order="C")
      cdef double[::1] new_b = np.zeros(n_hyperplanes, dtype=np.float64, order="C")
      cdef double[:,::1] T_matrix = np.zeros((n_variables, n_variables), dtype=np.float64, order="C")
      cdef double[::1] T_shift = np.zeros((n_variables), dtype=np.float64, order="C")
      cdef int N_samples
      cdef int N_ess = ess
      cdef bint check_psrf = bool(psrf_check) # restrict variables to {0,1} using Python's rules
      cdef bint parallel = bool(parallelism)
      
      self.polytope_cpp.mmcs_initialize(n_variables, ess, check_psrf, parallel, num_threads)

      # Get max inscribed ball for the initial polytope
      temp_center, radius = inner_ball(self._A, self._b, solver)
      cdef double[::1] inner_point_for_c = np.asarray(temp_center)

      while True:

         check = self.polytope_cpp.mmcs_step(&inner_point_for_c[0], radius, N_samples)
         
         if check > 1.0 and check < 2.0:
            break

         self.polytope_cpp.get_polytope_as_matrices(&new_A[0,0], &new_b[0])
         new_temp_c, radius = inner_ball(np.asarray(new_A), np.asarray(new_b), solver)
         inner_point_for_c = np.asarray(new_temp_c)
      
      cdef double[:,::1] samples = np.zeros((n_variables, N_samples), dtype=np.float64, order="C")
      self.polytope_cpp.get_mmcs_samples(&T_matrix[0,0], &T_shift[0], &samples[0,0])
      self.polytope_cpp.get_polytope_as_matrices(&new_A[0,0], &new_b[0])

      return np.asarray(new_A), np.asarray(new_b), np.asarray(T_matrix), np.asarray(T_shift), np.asarray(samples)

   def A(self):
      return np.asarray(self._A)

   def b(self):
      return np.asarray(self._b)

   def dimension(self):
      return self._A.shape[1]
