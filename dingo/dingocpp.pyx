# This is a cython wrapper for the C++ library volesti
# volesti (volume computation and sampling library)
  
# Copyright (c) 2012-2021 Vissarion Fisikopoulos
# Copyright (c) 2018-2021 Apostolos Chalkis
# Copyright (c) 2020-2021 Pedro Zuidberg Dos Martires
# Copyright (c) 2020-2021 Haris Zafeiropoulos

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
from libcpp cimport bool
from cpython cimport bool

# For the read the json format BIGG files function
import json
import scipy.io
# ----------------------------------------------------------------------------------

from dingo.inner_ball import slow_inner_ball
from dingo.gurobi_based_implementations import fast_inner_ball


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
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, \
         bool cdhr, bool rdhr, bool gaussian, bool set_L, bool accelerated_billiard, bool billiard, bool ball_walk, \
         double a, double L, bool max_ball, double* inner_point, double radius, double* samples);

      # Rounding H-Polytope
      void rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_value, \
         bool max_ball, double* inner_point, double radius);

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

#  This is where the volesti functions are getting their python interface; first the compute_volume() function
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
   def generate_samples(self, walk_len = 1, number_of_points = 1000, number_of_points_to_burn = 0, boundary = False, cdhr=False, \
      rdhr = False, gaussian = False, set_L = False, accelerated_billiard = True, billiard = False, ball_walk = False, a = 0, \
      radius = 0, inner_point = [], L = 0):

      n_variables = self._A.shape[1]
      cdef double[:,::1] samples = np.zeros((number_of_points,  n_variables), dtype = np.float64, order = "C")
      cdef double[::1] inner_point_for_c = np.asarray(inner_point)
      
      # Check whether the user asks for a certai value of radius; this is of higher priority than having a radius from the corresponding function
      if radius <= 0:        
        max_ball = False
      else:
         max_ball = True
            
      if L <= 0:
         set_L = False
      else:
         set_L = True
      
      self.polytope_cpp.generate_samples(walk_len, number_of_points, number_of_points_to_burn, boundary, cdhr, rdhr, gaussian, set_L, \
                                 accelerated_billiard, billiard, ball_walk, a, L, max_ball, &inner_point_for_c[0], radius, &samples[0,0])
      return np.asarray(samples)      # we need to build a Python function for getting a starting point depending on the polytope


# The rounding() function; like the compute_volume; there are more than one methods for this step
   def rounding(self, rounding_method = 'max_ellipsoid', inner_point = [], radius = 0):

      # Get the dimensions of the items about to build
      n_hyperplanes, n_variables = self._A.shape[0], self._A.shape[1]

      # Set the variables of those items; notice that they are all cdef type except of the last one which is about to be used
      # both as a C++ and a Python variable
      cdef double[:,::1] new_A = np.zeros((n_hyperplanes, n_variables), dtype=np.float64, order="C")
      cdef double[::1] new_b = np.zeros(n_hyperplanes, dtype=np.float64, order="C")
      cdef double[:,::1] T_matrix = np.zeros((n_variables, n_variables), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n_variables), dtype=np.float64, order="C")
      cdef double round_value
      
      cdef double[::1] inner_point_for_c = np.asarray(inner_point)
      
      # Transform the rounding_method variable to UTF-8 coding
      rounding_method = rounding_method.encode("UTF-8")

      # Check whether a max ball has been given
      if radius > 0:
         max_ball = True
      else:
         max_ball = False
      
      # Check whether the rounding method the user asked for, is actually among those volestipy supports
      if rounding_method in rounding_methods:

         self.polytope_cpp.rounding(rounding_method, &new_A[0,0], &new_b[0], &T_matrix[0,0], &shift[0], round_value, max_ball, &inner_point_for_c[0], radius)

         np.save('A_rounded.npy', new_A) ; np.save('b_rounded.npy', new_b)
         np.save('T_rounded.npy', T_matrix) ; np.save('shift_rounded.npy', shift)
         np.save('round_value.npy', np.asarray(round_value))

         return np.asarray(new_A),np.asarray(new_b),np.asarray(T_matrix),np.asarray(shift),np.asarray(round_value)

      else:

         raise Exception('"{}" is not implemented to walk types. Available methods are: {}'.format(rounding_method, rounding_methods))

   def A(self):
      return np.asarray(self._A)

   def b(self):
      return np.asarray(self._b)

   def dimension(self):
      return self._A.shape[1]
