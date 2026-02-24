# This is a cython wrapper for the C++ library volesti
# volesti (volume computation and sampling library)

# Copyright (c) 2012-2021 Vissarion Fisikopoulos
# Copyright (c) 2018-2021 Apostolos Chalkis
# Copyright (c) 2020-2021 Pedro Zuidberg Dos Martires
# Copyright (c) 2020-2021 Haris Zafeiropoulos
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

# Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

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
   void _generate_cube_H "generate_cube_H"(int dim, double scale, double* A_out, double* b_out)
   void _generate_simplex_H "generate_simplex_H"(int dim, double* A_out, double* b_out)
   void _generate_birkhoff_H "generate_birkhoff_H"(int n, double* A_out, double* b_out)

   cdef struct SBDiagnostics:
        double minESS
        double maxPSRF
        long long N

   # The HPolytopeCPP class along with its functions
   cdef cppclass HPolytopeCPP:

      # Initialization
      HPolytopeCPP() except +
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables) except +

      # Compute volume
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);

      # Random sampling
      double apply_sampling(int walk_len, int number_of_points, int number_of_points_to_burn, \
                            char* method, double* inner_point, double radius, double* samples, \
                            double variance_value, double* bias_vector, int ess,int nreflections)

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

      int apply_boundary_sampling(int walk_len,int number_of_points,int number_of_points_to_burn,const char* sampler,int nreflections,double* samples)

      void set_sb_state_from_buffer(int d, int N, const double* samples, const SBDiagnostics& diag)
      void get_sb_samples(double* out) const
      void get_sb_diagnostics(double* out3) const
      void boundary_scaling_ratio(int d, int N, const double* samples,double tol, double min_ratio,double* scale_out,double* coverage_out,double* max_dev_out,double* avg_dev_out,int* zero_count_out,double* zero_pct_out) const

      @staticmethod
      SBDiagnostics sb_diagnostics(int d, int N, const double* samples)


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
                "gaussian_CDHR".encode("UTF-8"), "gaussian_RDHR".encode("UTF-8"), "uniform_ball".encode("UTF-8"), "billiard".encode("UTF-8"),"shake_and_bake".encode("UTF-8"),"billiard_shake_and_bake".encode("UTF-8")  ]
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
   def generate_samples(self, method, number_of_points, number_of_points_to_burn, walk_len,
                        variance_value, bias_vector, solver = None, ess = 1000, nreflections=None):

      n_variables = self._A.shape[1]
      cdef double[:,::1] samples = np.zeros((number_of_points, n_variables), dtype = np.float64, order = "C")

      # Get max inscribed ball for the initial polytope
      temp_center, radius = inner_ball(self._A, self._b, solver)

      cdef double[::1] inner_point_for_c = np.asarray(temp_center)

      cdef double[::1] bias_vector_ = np.asarray(bias_vector)

      self.polytope_cpp.apply_sampling(walk_len, number_of_points, number_of_points_to_burn, \
                                       method, &inner_point_for_c[0], radius, &samples[0,0], \
                                       variance_value, &bias_vector_[0], ess, nreflections)
      return np.asarray(samples)

   def boundary_sample(self, sampler="sb", number_of_points=10000,number_of_points_to_burn=0, walk_len=1, nreflections=0):
      cdef int d = <int> self._A.shape[1]           
      cdef int N   = <int> number_of_points
      cdef int burn = <int> number_of_points_to_burn
      cdef int wl   = <int> walk_len
      cdef int nref = <int> nreflections
      cdef bytes sampler_b = sampler.encode("UTF-8")

      cdef np.ndarray[np.float64_t, ndim=2, mode="fortran"] S = \
         np.zeros((d, N), dtype=np.float64, order="F")

      cdef int made = self.polytope_cpp.apply_boundary_sampling(wl, N, burn, sampler_b, nref, &S[0,0])
      if made < 0:
         raise RuntimeError("apply_boundary_sampling failed")

      cdef SBDiagnostics diag = HPolytopeCPP.sb_diagnostics(d, made, <const double*> &S[0,0])
      self.polytope_cpp.set_sb_state_from_buffer(d, made, <const double*> &S[0,0], diag)

      return np.asarray(S)[:, :made]
   
   def get_sb_diagnostics(self, out):
      """
      out: np.ndarray float64, shape (>=3,), order:
            [minESS, maxPSRF, N]
      """
      cdef np.ndarray[np.float64_t, ndim=1] out_arr = np.ascontiguousarray(out, dtype=np.float64)
      if out_arr.shape[0] < 3:
         raise ValueError("out must have length >= 3")
      self.polytope_cpp.get_sb_diagnostics(<double*> &out_arr[0])
      return np.asarray(out_arr)


   def get_sb_samples(self):
      """
      Returns a d x N matrix from the internal C++ sample buffer.
      """
      cdef int d
      cdef long long Nll
      cdef np.ndarray[np.float64_t, ndim=1] tmp = np.zeros(3, dtype=np.float64)
      # Retrieve N from diagnostics (tmp = [minESS, maxPSRF, N])
      self.polytope_cpp.get_sb_diagnostics(<double*> &tmp[0])
      Nll = <long long> tmp[2]
      d = <int> self._A.shape[1]

      cdef np.ndarray[np.float64_t, ndim=2, mode="fortran"] S = \
         np.zeros((d, Nll), dtype=np.float64, order="F")
      self.polytope_cpp.get_sb_samples(&S[0,0])
      return np.asarray(S)


   def boundary_diag(self, S):
      """
      Compute diagnostics directly from a given sample matrix S.
      Returns a Python dictionary with keys: minESS, maxPSRF, and N.
      """
      cdef np.ndarray[np.float64_t, ndim=2] Snp = np.array(S, dtype=np.float64, order="F")
      cdef int d = <int> Snp.shape[0]
      cdef int N = <int> Snp.shape[1]
      cdef SBDiagnostics diag = HPolytopeCPP.sb_diagnostics(d, N, <const double*> &Snp[0,0])
      return {"minESS": diag.minESS, "maxPSRF": diag.maxPSRF, "N": diag.N}
   
   def boundary_scaling_ratio(self, S, double tol=1e-10, double min_ratio=0.01):
      cdef np.ndarray[np.float64_t, ndim=2, mode="fortran"] Snp = \
         np.array(S, dtype=np.float64, order="F")

      cdef int d = <int> Snp.shape[0]
      cdef int N = <int> Snp.shape[1]
      cdef int m = <int> self._A.shape[0]
      cdef int K = 10

      cdef np.ndarray[np.float64_t, ndim=1] scale = np.zeros(K, dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] max_dev = np.zeros(m, dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] avg_dev = np.zeros(m, dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] coverage = np.zeros((m, K), dtype=np.float64, order="C")

      cdef int zero_count = 0
      cdef double zero_pct = 0.0

      self.polytope_cpp.boundary_scaling_ratio(d, N, <const double*> &Snp[0,0],tol, min_ratio,<double*> &scale[0],<double*> &coverage[0,0],<double*> &max_dev[0],<double*> &avg_dev[0],&zero_count,&zero_pct)

      return (np.asarray(scale),np.asarray(coverage),np.asarray(max_dev),np.asarray(avg_dev),int(zero_count),float(zero_pct))

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

def generate_cube(int dim, double scale=1.0):
   """
   Create an H-polytope for a hypercube in R^dim with side length 2*scale
   centered at the origin: [-scale, scale]^dim.

   """
   cdef int m = 2 * dim
   cdef int n = dim

   cdef np.ndarray[np.float64_t, ndim=2] A = \
      np.zeros((m, n), dtype=np.float64)
   cdef np.ndarray[np.float64_t, ndim=1] b = \
      np.zeros(m, dtype=np.float64)

   _generate_cube_H(dim, scale, &A[0, 0], &b[0])

   return np.asarray(A), np.asarray(b)


def generate_simplex(int dim):
   """
   Create an H-polytope for a standard simplex in R^dim.

   """
   cdef int m = dim + 1
   cdef int n = dim

   cdef np.ndarray[np.float64_t, ndim=2] A = \
      np.zeros((m, n), dtype=np.float64)
   cdef np.ndarray[np.float64_t, ndim=1] b = \
      np.zeros(m, dtype=np.float64)

   _generate_simplex_H(dim, &A[0, 0], &b[0])

   return np.asarray(A), np.asarray(b)

def generate_birkhoff(int n):
   """
   Creates an H-polytope for the Birkhoff polytope of size n.

   """
   cdef int m = n * n
   cdef int d = n * n - 2 * n + 1

   cdef np.ndarray[np.float64_t, ndim=2] A = \
      np.zeros((m, d), dtype=np.float64)
   cdef np.ndarray[np.float64_t, ndim=1] b = \
      np.zeros(m, dtype=np.float64)

   _generate_birkhoff_H(n, &A[0, 0], &b[0])

   return np.asarray(A), np.asarray(b)
