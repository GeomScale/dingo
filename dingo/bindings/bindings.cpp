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
                                    int method,
                                    double* inner_point, 
                                    double radius,
                                    double* samples,
                                    double variance_value,
                                    double* bias_vector_){ 
                                                                    
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

   if (method == 1) { // cdhr
      uniform_sampling<CDHRWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                 starting_point, number_of_points_to_burn);
   } else if (method == 2) { // rdhr
      uniform_sampling<RDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, 
                                 starting_point, number_of_points_to_burn);
   } else if (method == 3) { // accelerated_billiard
      uniform_sampling<AcceleratedBilliardWalk>(rand_points, HP, rng, walk_len, 
                                                number_of_points, starting_point, 
                                                number_of_points_to_burn);
   } else if (method == 4) { // ball walk
      uniform_sampling<BallWalk>(rand_points, HP, rng, walk_len, number_of_points, 
                                 starting_point, number_of_points_to_burn);
   } else if (method == 5) { // dikin walk
      uniform_sampling<DikinWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                  starting_point, number_of_points_to_burn);
   } else if (method == 6) { // john walk
      uniform_sampling<JohnWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                 starting_point, number_of_points_to_burn);
   } else if (method == 7) { // vaidya walk
      uniform_sampling<VaidyaWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                   starting_point, number_of_points_to_burn);
   } else if (method == 8) { // Gaussian sampling with exact HMC walk
      NT a = NT(1)/(NT(2)*variance);
      gaussian_sampling<GaussianHamiltonianMonteCarloExactWalk>(rand_points, HP, rng, walk_len, number_of_points, a,
                                   starting_point, number_of_points_to_burn);
   } else if (method == 9) { // exponential sampling with exact HMC walk
      VT c(d);
      for (int i = 0; i < d; i++){
         c(i) = bias_vector_[i];
      }
      Point bias_vector(c);       
      exponential_sampling<ExponentialHamiltonianMonteCarloExactWalk>(rand_points, HP, rng, walk_len, number_of_points, bias_vector, variance,
                                   starting_point, number_of_points_to_burn);                                
   } else if (method == 10) { // HMC with Gaussian distribution         
      rand_points = hmc_leapfrog_gaussian(walk_len, number_of_points, number_of_points_to_burn, variance, starting_point, HP);                             
   } else if (method == 11) { // HMC with exponential distribution   
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

   // The following block of code allows us to copy the sampled points
   auto n_si=0;
   for (auto it_s = rand_points.cbegin(); it_s != rand_points.cend(); it_s++){
      for (auto i = 0; i != it_s->dimension(); i++){
         samples[n_si++] = (*it_s)[i];
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


void HPolytopeCPP::mmcs_initialize(int d, int ess, bool psrf_check, bool parallelism, int num_threads) {

   mmcs_set_of_parameters = mmcs_params(d, ess, psrf_check, parallelism, num_threads);

}


double HPolytopeCPP::mmcs_step(double* inner_point, double radius, int &N) {
   
   HP.normalize();
   int d = HP.dimension();

   VT inner_vec(d);
   NT max_s;
   
   for (int i = 0; i < d; i++){
      inner_vec(i) = inner_point[i];
   }
   
   Point inner_point2(inner_vec);   
   CheBall = std::pair<Point, NT>(inner_point2, radius);   
      
   HP.set_InnerBall(CheBall);   

   RNGType rng(d);

   
   if (mmcs_set_of_parameters.request_rounding && mmcs_set_of_parameters.rounding_completed) 
   {
      mmcs_set_of_parameters.req_round_temp = false;
   }

   if (mmcs_set_of_parameters.req_round_temp) 
   {
      mmcs_set_of_parameters.nburns = mmcs_set_of_parameters.num_rounding_steps / mmcs_set_of_parameters.window + 1;
   } 
   else 
   {
      mmcs_set_of_parameters.nburns = mmcs_set_of_parameters.max_num_samples / mmcs_set_of_parameters.window + 1;
   }

   NT L = NT(6) * std::sqrt(NT(d)) * CheBall.second;
   AcceleratedBilliardWalk WalkType(L);

   unsigned int Neff_sampled;
   MT TotalRandPoints;
   if (mmcs_set_of_parameters.parallelism)
   {
      mmcs_set_of_parameters.complete = perform_parallel_mmcs_step<AcceleratedBilliardWalkParallel>(HP, rng, mmcs_set_of_parameters.walk_length, 
                                                                                                    mmcs_set_of_parameters.Neff, 
                                                                                                    mmcs_set_of_parameters.max_num_samples, 
                                                                                                    mmcs_set_of_parameters.window, 
                                                                                                    Neff_sampled, mmcs_set_of_parameters.total_samples, 
                                                                                                    mmcs_set_of_parameters.num_rounding_steps, 
                                                                                                    TotalRandPoints, CheBall.first, mmcs_set_of_parameters.nburns, 
                                                                                                    mmcs_set_of_parameters.num_threads, 
                                                                                                    mmcs_set_of_parameters.req_round_temp, L);
   }
   else
   {
      mmcs_set_of_parameters.complete = perform_mmcs_step(HP, rng, mmcs_set_of_parameters.walk_length, mmcs_set_of_parameters.Neff,
                                                          mmcs_set_of_parameters.max_num_samples, mmcs_set_of_parameters.window,
                                                          Neff_sampled, mmcs_set_of_parameters.total_samples, 
                                                          mmcs_set_of_parameters.num_rounding_steps, TotalRandPoints, CheBall.first,
                                                          mmcs_set_of_parameters.nburns, mmcs_set_of_parameters.req_round_temp, WalkType);
   }

   mmcs_set_of_parameters.store_ess(mmcs_set_of_parameters.phase) = Neff_sampled;
   mmcs_set_of_parameters.store_nsamples(mmcs_set_of_parameters.phase) = mmcs_set_of_parameters.total_samples;
   mmcs_set_of_parameters.phase++;
   mmcs_set_of_parameters.Neff -= Neff_sampled;
   std::cout << "phase " << mmcs_set_of_parameters.phase << ": number of correlated samples = " << mmcs_set_of_parameters.total_samples << ", effective sample size = " << Neff_sampled;
   mmcs_set_of_parameters.total_neff += Neff_sampled;
   
   mmcs_set_of_parameters.samples.conservativeResize(d, mmcs_set_of_parameters.total_number_of_samples_in_P0 + mmcs_set_of_parameters.total_samples);
   for (int i = 0; i < mmcs_set_of_parameters.total_samples; i++)
   {
      mmcs_set_of_parameters.samples.col(i + mmcs_set_of_parameters.total_number_of_samples_in_P0) = 
                              mmcs_set_of_parameters.T * TotalRandPoints.row(i).transpose() + mmcs_set_of_parameters.T_shift;
   }

   N = mmcs_set_of_parameters.total_number_of_samples_in_P0 + mmcs_set_of_parameters.total_samples;   
   mmcs_set_of_parameters.total_number_of_samples_in_P0 += mmcs_set_of_parameters.total_samples;

   if (!mmcs_set_of_parameters.complete) 
   {
      if (mmcs_set_of_parameters.request_rounding && !mmcs_set_of_parameters.rounding_completed) 
      {
         VT shift(d), s(d);
         MT V(d, d), S(d, d), round_mat;
         for (int i = 0; i < d; ++i) 
         {
            shift(i) = TotalRandPoints.col(i).mean();
         }

         for (int i = 0; i < mmcs_set_of_parameters.total_samples; ++i) 
         {
            TotalRandPoints.row(i) = TotalRandPoints.row(i) - shift.transpose();
         }

         Eigen::BDCSVD<MT> svd(TotalRandPoints, Eigen::ComputeFullV);
         s = svd.singularValues() / svd.singularValues().minCoeff();

         if (s.maxCoeff() >= 2.0) 
         {
            for (int i = 0; i < s.size(); ++i) 
            {
               if (s(i) < 2.0) 
               {
                  s(i) = 1.0;
               }
            }
            V = svd.matrixV();
         } 
         else 
         {
            s = VT::Ones(d);
            V = MT::Identity(d, d);
         }
         max_s = s.maxCoeff();
         S = s.asDiagonal();
         round_mat = V * S;

         mmcs_set_of_parameters.round_it++;
         HP.shift(shift);
         HP.linear_transformIt(round_mat);
         mmcs_set_of_parameters.T_shift += mmcs_set_of_parameters.T * shift;
         mmcs_set_of_parameters.T = mmcs_set_of_parameters.T * round_mat;

         std::cout << ", ratio of the maximum singilar value over the minimum singular value = " << max_s << std::endl;

         if (max_s <= mmcs_set_of_parameters.s_cutoff || mmcs_set_of_parameters.round_it > mmcs_set_of_parameters.num_its) 
         {
            mmcs_set_of_parameters.rounding_completed = true;
         }
      }
      else 
      {
         std::cout<<"\n";
      }
   } 
   else if (!mmcs_set_of_parameters.psrf_check)
   {
      NT max_psrf = univariate_psrf<NT, VT>(mmcs_set_of_parameters.samples).maxCoeff();
      std::cout << "\n[5]total ess " << mmcs_set_of_parameters.total_neff << ": number of correlated samples = " << mmcs_set_of_parameters.samples.cols()<<std::endl;
      std::cerr << "[5]maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(mmcs_set_of_parameters.samples).maxCoeff() << std::endl;
      std::cout<<"\n\n";
      return 1.5;
   }
   else 
   {
      TotalRandPoints.resize(0, 0);
      NT max_psrf = univariate_psrf<NT, VT>(mmcs_set_of_parameters.samples).maxCoeff();

      if (max_psrf < 1.1 && mmcs_set_of_parameters.total_neff >= mmcs_set_of_parameters.fixed_Neff) {
         std::cout << "\n[4]total ess " << mmcs_set_of_parameters.total_neff << ": number of correlated samples = " << mmcs_set_of_parameters.samples.cols()<<std::endl;
         std::cerr << "[4]maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(mmcs_set_of_parameters.samples).maxCoeff() << std::endl;
         std::cout<<"\n\n";
         return 1.5;
      }
      std::cerr << "\n [1]maximum marginal PSRF: " <<  max_psrf << std::endl;

      while (max_psrf > 1.1 && mmcs_set_of_parameters.total_neff >= mmcs_set_of_parameters.fixed_Neff) {
         
         mmcs_set_of_parameters.Neff += mmcs_set_of_parameters.store_ess(mmcs_set_of_parameters.skip_phase);
         mmcs_set_of_parameters.total_neff -= mmcs_set_of_parameters.store_ess(mmcs_set_of_parameters.skip_phase);

         mmcs_set_of_parameters.total_number_of_samples_in_P0 -= mmcs_set_of_parameters.store_nsamples(mmcs_set_of_parameters.skip_phase);
         N -= mmcs_set_of_parameters.store_nsamples(mmcs_set_of_parameters.skip_phase);

         MT S = mmcs_set_of_parameters.samples;
         mmcs_set_of_parameters.samples.resize(d, mmcs_set_of_parameters.total_number_of_samples_in_P0);
         mmcs_set_of_parameters.samples = 
                     S.block(0, mmcs_set_of_parameters.store_nsamples(mmcs_set_of_parameters.skip_phase), d, mmcs_set_of_parameters.total_number_of_samples_in_P0);

         mmcs_set_of_parameters.skip_phase++;

         max_psrf = univariate_psrf<NT, VT>(mmcs_set_of_parameters.samples).maxCoeff();

         std::cerr << "[2]maximum marginal PSRF: " <<  max_psrf << std::endl;
         std::cerr << "[2]total ess: " <<  mmcs_set_of_parameters.total_neff << std::endl;

         if (max_psrf < 1.1 && mmcs_set_of_parameters.total_neff >= mmcs_set_of_parameters.fixed_Neff) {
            return 1.5;
         }
      }
      std::cout << "[3]total ess " << mmcs_set_of_parameters.total_neff << ": number of correlated samples = " << mmcs_set_of_parameters.samples.cols()<<std::endl;
      std::cerr << "[3]maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(mmcs_set_of_parameters.samples).maxCoeff() << std::endl;
      std::cout<<"\n\n";
      return 0.0;
   }

   return 0.0;
}


void HPolytopeCPP::get_mmcs_samples(double* T_matrix, double* T_shift, double* samples) {

   int n_variables = HP.dimension();

   int t_mat_index = 0;
   for (int i = 0; i < n_variables; i++){
      for (int j = 0; j < n_variables; j++){
         T_matrix[t_mat_index++] = mmcs_set_of_parameters.T(i, j);
      }
   }
   
   // create the shift vector
   for (int i = 0; i < n_variables; i++){
      T_shift[i] = mmcs_set_of_parameters.T_shift[i];
   }

   int N = mmcs_set_of_parameters.samples.cols();

   int t_si = 0;
   for (int i = 0; i < n_variables; i++){
      for (int j = 0; j < N; j++){
         samples[t_si++] = mmcs_set_of_parameters.samples(i, j);
      }
   }
   mmcs_set_of_parameters.samples.resize(0,0);
}


//////////         Start of "rounding()"          //////////
void HPolytopeCPP::apply_rounding(int rounding_method, double* new_A, double* new_b,
                                  double* T_matrix, double* shift, double &round_value,
                                  double* inner_point, double radius){

   // make a copy of the initial HP which will be used for the rounding step
   auto P(HP);
   RNGType rng(P.dimension());
   P.normalize();
   
   
      
   // read the inner point provided by the user and the radius
   int d = P.dimension();
   VT inner_vec(d);
      
   for (int i = 0; i < d; i++){
      inner_vec(i) = inner_point[i];
   }

   Point inner_point2(inner_vec);
   CheBall = std::pair<Point, NT>(inner_point2, radius);
   P.set_InnerBall(CheBall);
   
   // set the output variable of the rounding step
   round_result round_res;

   // walk length will always be equal to 2
   int walk_len = 2;

   // run the rounding method
   if (rounding_method == 1) { // max ellipsoid
      round_res = max_inscribed_ellipsoid_rounding<MT, VT, NT>(P, CheBall.first);
      
   } else if (rounding_method == 2) { // isotropization
      round_res = svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, CheBall, 1, rng);
   } else if (rounding_method == 3) { // min ellipsoid
      round_res = min_sampling_covering_ellipsoid_rounding<AcceleratedBilliardWalk, MT, VT>(P,
                                                                                            CheBall,
                                                                                            walk_len,
                                                                                            rng);
   } else {
      throw std::runtime_error("Unknown rounding method.");
   }

   // create the new_A matrix
   MT A_to_copy = P.get_mat();
   int n_hyperplanes = P.num_of_hyperplanes();
   int n_variables = P.dimension();

   auto n_si = 0;
   for (int i = 0; i < n_hyperplanes; i++){
      for (int j = 0; j < n_variables; j++){
         new_A[n_si++] = A_to_copy(i,j);
      }
   }

   // create the new_b vector
   VT new_b_temp = P.get_vec();
   for (int i=0; i < n_hyperplanes; i++){
      new_b[i] = new_b_temp[i];
   }

   // create the T matrix
   MT T_matrix_temp = get<0>(round_res);
   auto t_si = 0;
   for (int i = 0; i < n_variables; i++){
      for (int j = 0; j < n_variables; j++){
         T_matrix[t_si++] = T_matrix_temp(i,j);
      }
   }

   // create the shift vector
   VT shift_temp = get<1>(round_res);
   for (int i = 0; i < n_variables; i++){
      shift[i] = shift_temp[i];
   }

   // get the round val value from the rounding step and print int
   round_value = get<2>(round_res);

}
//////////         End of "rounding()"          //////////
