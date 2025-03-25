// This is binding file for the C++ library volesti
// volesti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Haris Zafeiropoulos
// Contributed and/or modified by Pedro Zuidberg Dos Martires

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H

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
#include "ode_solvers/ode_solvers.hpp"

// for rounding
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"

typedef double NT;
typedef Cartesian<NT>    Kernel;
typedef typename Kernel::Point    Point;
typedef HPolytope<Point> Hpolytope;
typedef typename Hpolytope::MT    MT;
typedef typename Hpolytope::VT    VT;
typedef BoostRandomNumberGenerator<boost::mt19937, double>    RNGType;

// This is the HPolytopeCPP class; the main volesti class that is running the compute_volume(), rounding() and sampling() methods
class HPolytopeCPP{

   public:

      std::pair<Point,NT> CheBall;

      // The class and its main specs
      HPolytopeCPP();
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables);

      Hpolytope HP;
      // Here we use the "~" destructor; this way we avoid a memory leak.
      ~HPolytopeCPP();

      // the compute_volume() function
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed) const;

      // the apply_sampling() function
      double apply_sampling(int walk_len, int number_of_points, int number_of_points_to_burn,
                            char* method, double* inner_point, double radius, double* samples,
                            double variance_value, double* bias_vector, int ess);

      void get_polytope_as_matrices(double* new_A, double* new_b) const;
};

#endif
