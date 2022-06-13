// This is binding file for the C++ library volesti
// volesti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Haris Zafeiropoulos
// Contributed and/or modified by Pedro Zuidberg Dos Martires

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H

#define DISABLE_LPSOLVE
#define DISABLE_NLP_ORACLES
#include <cmath>
// from SOB volume - exactly the same for CG and CB methods
#include <fstream>
#include <iostream>
#include "random_walks.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "sampling/mmcs.hpp"
#include "sampling/parallel_mmcs.hpp"
#include "diagnostics/univariate_psrf.hpp"

//from generate_samples, some extra headers not already included
#include <chrono>
#include "sampling/sampling.hpp"

// for rounding
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"
#include "preprocess/get_full_dimensional_polytope.hpp"

typedef double NT;
typedef Cartesian<NT>    Kernel;
typedef typename Kernel::Point    Point;
typedef HPolytope<Point> Hpolytope;
typedef typename Hpolytope::MT    MT;
typedef typename Hpolytope::VT    VT;
typedef BoostRandomNumberGenerator<boost::mt19937, double>    RNGType;


template <typename NT, typename MT, typename VT>
struct mmcs_parameters
{
public:

   mmcs_parameters() {}

   mmcs_parameters(int d, int ess, bool _psrf_check, bool _parallelism, int _num_threads)
         :  T(MT::Identity(d,d))
         ,  T_shift(VT::Zero(d))
         ,  store_ess(VT::Zero(50))
         ,  store_nsamples(VT::Zero(50))
         ,  skip_phase(0)
         ,  num_rounding_steps(20*d)
         ,  walk_length(1)
         ,  num_its(20)
         ,  Neff(ess)
         ,  fixed_Neff(ess)
         ,  phase(0)
         ,  window(100)
         ,  max_num_samples(100 * d)
         ,  round_it(1)
         ,  total_number_of_samples_in_P0(0)
         ,  total_neff(0)
         ,  num_threads(_num_threads)
         ,  psrf_check(_psrf_check)
         ,  parallelism(_parallelism)
         ,  complete(false)
         ,  request_rounding(true)
         ,  rounding_completed(false)
         ,  s_cutoff(NT(3))
   {
      req_round_temp = request_rounding;
   }

   MT T;
   MT samples;
   VT T_shift;
   VT store_ess;
   VT store_nsamples;
   unsigned int skip_phase;
   unsigned int num_rounding_steps;
   unsigned int walk_length;
   unsigned int num_its;
   int Neff;
   int fixed_Neff;
   unsigned int phase;
   unsigned int window;
   unsigned int max_num_samples;
   unsigned int total_samples;
   unsigned int nburns;
   unsigned int round_it;
   unsigned int total_number_of_samples_in_P0;
   unsigned int total_neff;
   unsigned int num_threads;
   bool psrf_check;
   bool parallelism;
   bool complete;
   bool request_rounding;
   bool rounding_completed;
   bool req_round_temp;
   NT s_cutoff;
};


// This is the HPolytopeCPP class; the main volesti class that is running the compute_volume(), rounding() and sampling() methods
class HPolytopeCPP{

   public:

      std::pair<Point,NT> CheBall;

      // regarding the rounding step
      typedef std::tuple<MT, VT, NT>    round_result;
      typedef mmcs_parameters<NT, MT, VT> mmcs_params;

      mmcs_params mmcs_set_of_parameters;

      // The class and its main specs
      HPolytopeCPP();
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables);

      Hpolytope HP;
      // Here we use the "~" destructor; this way we avoid a memory leak.
      ~HPolytopeCPP();

      // the compute_volume() function
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed) const;

      // the generate_samples() function
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, 
       bool cdhr, bool rdhr, bool gaussian, bool set_L, bool accelerated_billiard, bool billiard, bool ball_walk, double a, double L,  
       bool max_ball, double* inner_point, double radius,
       double* samples);

      void mmcs_initialize(int d, int ess, bool psrf_check, bool parallelism, int num_threads);

      double mmcs_step(double* inner_point_for_c, double radius, int &N);

      void get_mmcs_samples(double* T_matrix, double* T_shift, double* samples);

      void get_polytope_as_matrices(double* new_A, double* new_b) const;

      // the rounding() function
      void rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_value,
       bool max_ball, double* inner_point, double radius);
      
};


#endif 
