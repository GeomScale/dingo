// This is binding file for the C++ library volesti
// volesti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Haris Zafeiropoulos
// Contributed and/or modified by Pedro Zuidberg Dos Martires

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <math.h>
#include <stdexcept>
#include "bindings.h"
#include "hmc_sampling.h"
#include "sampling/mmcs.hpp"

using namespace std;

// >>> Main HPolytopeCPP class; compute_volume(), rounding() and generate_samples() volesti methods are included <<<

// Here is the initialization of the HPolytopeCPP class
HPolytopeCPP::HPolytopeCPP() {}
HPolytopeCPP::HPolytopeCPP(double *A_np, double *b_np, int n_hyperplanes, int n_variables){

   MT A;
   VT b;
   A.resize(n_hyperplanes,n_variables);
   b.resize(n_hyperplanes);

   int index = 0;
   for (int i = 0; i < n_hyperplanes; i++){
      b(i) = b_np[i];
      for (int j=0; j < n_variables; j++){
         A(i,j) = A_np[index];
         index++;
      }
   }

   HP = Hpolytope(n_variables, A, b);
}
// Use a destructor for the HPolytopeCPP object
HPolytopeCPP::~HPolytopeCPP(){}

//////////          Start of "compute_volume"          //////////
double HPolytopeCPP::compute_volume(char* vol_method, char* walk_method,
                                    int walk_len, double epsilon, int seed) const {

   double volume;

   if (strcmp(vol_method,"sequence_of_balls") == 0){
      if (strcmp(walk_method,"uniform_ball") == 0){
         volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"CDHR") == 0){
         volume = volume_sequence_of_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"RDHR") == 0){
         volume = volume_sequence_of_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   }
   else if (strcmp(vol_method,"cooling_gaussian") == 0){
      if (strcmp(walk_method,"gaussian_ball") == 0){
         volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_CDHR") == 0){
         volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_RDHR") == 0){
         volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   } else if (strcmp(vol_method,"cooling_balls") == 0){
       if (strcmp(walk_method,"uniform_ball") == 0){
         volume = volume_cooling_balls<BallWalk, RNGType>(HP, epsilon, walk_len).second;
       } else if (strcmp(walk_method,"CDHR") == 0){
         volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len).second;
       } else if (strcmp(walk_method,"RDHR") == 0){
         volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len).second;
       } else if (strcmp(walk_method,"billiard") == 0){
         volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, epsilon, walk_len).second;
       }
   }
   return volume;
}
//////////           End of "compute_volume()"            //////////


//////////         Start of "generate_samples()"          //////////
double HPolytopeCPP::apply_sampling(int walk_len,
                                    int number_of_points,
                                    int number_of_points_to_burn,
                                    char* method,
                                    double* inner_point,
                                    double radius,
                                    double* samples,
                                    double variance_value,
                                    double* bias_vector_,
                                    int ess){

   RNGType rng(HP.dimension());
   HP.normalize();
   int d = HP.dimension();
   Point starting_point;
   VT inner_vec(d);

   for (int i = 0; i < d; i++){
      inner_vec(i) = inner_point[i];
   }

   Point inner_point2(inner_vec);
   CheBall = std::pair<Point, NT>(inner_point2, radius);
   HP.set_InnerBall(CheBall);
   starting_point = inner_point2;
   std::list<Point> rand_points;

   NT variance = variance_value;

   if (strcmp(method, "cdhr")) { // cdhr
      uniform_sampling<CDHRWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                 starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "rdhr")) { // rdhr
      uniform_sampling<RDHRWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                 starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "billiard_walk")) { // accelerated_billiard
      uniform_sampling<AcceleratedBilliardWalk>(rand_points, HP, rng, walk_len,
                                                number_of_points, starting_point,
                                                number_of_points_to_burn);
   } else if (strcmp(method, "ball_walk")) { // ball walk
      uniform_sampling<BallWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                 starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "dikin_walk")) { // dikin walk
      uniform_sampling<DikinWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                  starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "john_walk")) { // john walk
      uniform_sampling<JohnWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                 starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "vaidya_walk")) { // vaidya walk
      uniform_sampling<VaidyaWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                   starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "mmcs") == 0) {
       // Use volesti's MMCS implementation
       MT S;
       int total_neff = 0;
       mmcs(HP, ess, S, total_neff, walk_len, rng);
       
       // Copy results to output array
       for (int i = 0; i < d; i++) {
           for (int j = 0; j < S.cols(); j++) {
               samples[i * S.cols() + j] = S(i, j);
           }
       }
       return total_neff;
   } else if (strcmp(method, "gaussian_hmc_walk")) { // Gaussian sampling with exact HMC walk
      NT a = NT(1)/(NT(2)*variance);
      gaussian_sampling<GaussianHamiltonianMonteCarloExactWalk>(rand_points, HP, rng, walk_len, number_of_points, a,
                                   starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "exponential_hmc_walk")) { // exponential sampling with exact HMC walk
      VT c(d);
      for (int i = 0; i < d; i++){
         c(i) = bias_vector_[i];
      }
      Point bias_vector(c);
      exponential_sampling<ExponentialHamiltonianMonteCarloExactWalk>(rand_points, HP, rng, walk_len, number_of_points, bias_vector, variance,
                                   starting_point, number_of_points_to_burn);
   } else if (strcmp(method, "hmc_leapfrog_gaussian")) { // HMC with Gaussian distribution
      rand_points = hmc_leapfrog_gaussian(walk_len, number_of_points, number_of_points_to_burn, variance, starting_point, HP);
   } else if (strcmp(method, "hmc_leapfrog_exponential")) { // HMC with exponential distribution
      VT c(d);
      for (int i = 0; i < d; i++) {
         c(i) = bias_vector_[i];
      }
      Point bias_vector(c);

      rand_points = hmc_leapfrog_exponential(walk_len, number_of_points, number_of_points_to_burn, variance, bias_vector, starting_point, HP);

   }

   else {
      throw std::runtime_error("This function must not be called.");
   }

   if (!strcmp(method, "mmcs")) {
    // The following block of code allows us to copy the sampled points
    auto n_si=0;
    for (auto it_s = rand_points.cbegin(); it_s != rand_points.cend(); it_s++){
        for (auto i = 0; i != it_s->dimension(); i++){
            samples[n_si++] = (*it_s)[i];
        }
    }
   }
   return 0.0;
}
//////////         End of "generate_samples()"          //////////


void HPolytopeCPP::get_polytope_as_matrices(double* new_A, double* new_b) const {

   int n_hyperplanes = HP.num_of_hyperplanes();
   int n_variables = HP.dimension();

   int n_si = 0;
   MT A_to_copy = HP.get_mat();
   for (int i = 0; i < n_hyperplanes; i++){
      for (int j = 0; j < n_variables; j++){
         new_A[n_si++] = A_to_copy(i, j);
      }
   }

   // create the new_b vector
   VT new_b_temp = HP.get_vec();
   for (int i=0; i < n_hyperplanes; i++){
      new_b[i] = new_b_temp[i];
   }
}

