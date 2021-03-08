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

# For the preprocess step, we need the following dependencies
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

# For the read the json format BIGG files function
import json
import scipy.io
# ----------------------------------------------------------------------------------

#from .PySPQR.sparseqr import *
#from .PySPQR.sparseqr.sparseqr import *

#import subroutines
#from . import fba, fva, inner_ball, nullspace, scaling
from .fva import slow_fva, fast_fva
from .inner_ball import slow_inner_ball, fast_inner_ball
from .nullspace import nullspace_dense, nullspace_sparse
from .scaling import gmscale, apply_scaling, remove_almost_redundant_facets


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

   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def dimensions(self):
      return self._A.shape[1]

# Build the low_dim_polytope_cpp class
cdef class low_dim_HPolytope:

   cdef lowDimHPolytopeCPP low_dim_polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b
   cdef double [:,::1] _Aeq
   cdef double[::1] _beq

# Set the specs of the class
   def __cinit__(self, double[:,::1] A, double[::1] b, double[:,::1] Aeq, double[::1] beq):

      self._A = A
      self._b = b
      self._Aeq = Aeq
      self._beq = beq
      n_rows_of_A, n_cols_of_A = A.shape[0], A.shape[1]
      n_row_of_Aeq, n_cols_of_Aeq = Aeq.shape[0], Aeq.shape[1]

      # If statements to check whether the user's input is valid for the low_dim_HPolytope class to run
      if n_rows_of_A == b.shape[0]:

         if n_row_of_Aeq == beq.shape[0]:

            if n_cols_of_A == n_cols_of_Aeq:

               # Run the constructor
               self.low_dim_polytope_cpp = lowDimHPolytopeCPP(&A[0,0], &b[0], &Aeq[0,0], &beq[0], n_rows_of_A, n_cols_of_A, n_row_of_Aeq, n_cols_of_Aeq)

            else:
               raise Exception('The number of columns of A equals to "{}" while those of Aeq {}. \
                               A and Aeq need to have the same number of columns'.format(n_cols_of_A, n_cols_of_Aeq))
         else:
            raise Exception('The number of rows of Aeq equals to "{}" while the elements of the beq vector are {}. \
                            The beq vector needs to have length equal to the number of rows of Aeq.'.format(n_row_of_Aeq, beq.shape[0]))
      else:
         raise Exception('The number of rows of A equals to "{}" while the elements of b are {}. \
                         The b vector needs to have length equal to the number of rows of A.'.format(n_rows_of_A, b.shape[0]))

   # The get_full_dimensional_polytope() function(); that needs to run in case the user does not provide volestipy with a full dimensional polytope
   def full_dimensiolal_polytope(self):

      # Get dimensions of the initial S (Aeq) matrix
      m = self._Aeq.shape[0]
      n = self._Aeq.shape[1]
      k = self._A.shape[0]

      # Set the output variables
      # The number of lines in the transpose N (columns in the actual matrix) are at least n-m; but we do not know their exact number
      # So we initialize it with the maximum possible number of lines (n). the same is for the full A transpose matrix
      # Later, we will have to keep their actual dimension and remove these variables with the extra lines
      cdef double[:,::1] N_extra_trans = np.zeros((n, n), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n), dtype=np.float64, order="C")
      cdef double[:,::1] A_full_extra_trans = np.zeros((n,k), dtype=np.float64, order="C")
      cdef double[::1] b_full = np.zeros((k), dtype=np.float64, order="C")

      # We need to keep the final number of columns of the N / full_A matrices
      cpdef int n_of_cols_in_N

      # Call the C++ class to get the full_dimensional polytope
      n_of_cols_in_N = self.low_dim_polytope_cpp.full_dimensiolal_polytope(&N_extra_trans[0,0], &shift[0], &A_full_extra_trans[0,0], &b_full[0])

      # Get a matrix with exactly the number of lines and columns that N expands to and delete the one with the extra columns
      N = np.zeros((n, n_of_cols_in_N), dtype=np.float64, order="C")
      for i in range(n):
         for j in range(n_of_cols_in_N):
            N[i,j] = np.asarray(N_extra_trans[j,i])
      del N_extra_trans

      # Likewise, for the A matrix of the full dimensional polytope
      A_full = np.zeros((k, n_of_cols_in_N), dtype=np.float64, order="C")
      for i in range(k):
         for j in range(n_of_cols_in_N):
            A_full[i,j] = np.asarray(A_full_extra_trans[j,i])
      del A_full_extra_trans

      # Finally, we need to build an HP object for the full dumensional polytope we got
      full_dimensional_polytope = HPolytope(A_full,b_full)

      # Print all the output of the function in .npy files
      np.save('A_full_dim.npy', A_full) ; np.save('b_full_dim.npy', b_full)
      np.save('N_full_dim.npy', N) ; np.save('shift_full_dim.npy',shift)

      # Delete all non-needed vars
      del A_full
      del b_full

      # Return a tuple whith the full dimensional HPolytope object in the first position ([0]) the N matrix and the shift vector
      return full_dimensional_polytope, np.asarray(N), np.asarray(shift)


   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def Aeq(self):
      return np.asarray(self._Aeq)
   @property
   def beq(self):
      return np.asarray(self._beq)
   @property
   def dimensions(self):
      return self._A.shape[1]
   
