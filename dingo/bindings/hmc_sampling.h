#include "bindings.h" 
#include "ode_solvers/ode_solvers.hpp"
#include <iostream> 
#include <math.h> 
#include <stdexcept> 
using namespace std;

   

template<class Polytope> list<Point> hmc_leapfrog_gaussian(int walk_len,
                                    int number_of_points, 
                                    int number_of_points_to_burn, 
                                    double* samples,
                                    double variance,
                                    double* bias_vector_,
                                    double* inner_point,
                                    double radius,
                                    Polytope HP) {
                                    
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
   typedef GaussianFunctor::GradientFunctor<Point> NegativeGradientFunctor;
   typedef GaussianFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
   typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
   typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;
   unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
   RandomNumberGenerator rng(rng_seed); 

   GaussianFunctor::parameters<NT, Point> params(starting_point, 2 / (variance * variance), NT(-1));
   NegativeGradientFunctor F(params);
   NegativeLogprobFunctor f(params);
   HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, d);

   HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>hmc(&HP, starting_point, F, f, hmc_params);

   // burning points 
   for (int i = 0; i < number_of_points_to_burn ; i++) 
   hmc.apply(rng, walk_len);
 

   // actual sampling 
   for (int i = 0; i < number_of_points ; i++) { 
      hmc.apply(rng, walk_len); 
      rand_points.push_back(hmc.x); 
   }
   return rand_points;                                  
} 
 

template<class Polytope> list<Point> hmc_leapfrog_exponential(int walk_len,
                                    int number_of_points, 
                                    int number_of_points_to_burn, 
                                    double* samples,
                                    double variance,
                                    double* bias_vector_,
                                    double* inner_point,
                                    double radius,
                                    Polytope HP) {

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
   typedef ExponentialFunctor::GradientFunctor<Point> NegativeGradientFunctor;
   typedef ExponentialFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
   typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
   typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;
      
   unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
      RandomNumberGenerator rng(rng_seed); 
      
   VT c(d);
   for (int i = 0; i < d; i++) {
      c(i) = bias_vector_[i];
   }
   Point bias_vector(c);

   ExponentialFunctor::parameters<NT, Point> params(bias_vector, 2 / (variance * variance));
      
   NegativeGradientFunctor F(params);
   NegativeLogprobFunctor f(params);
   HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, d);
      
   HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>hmc(&HP, starting_point, F, f, hmc_params);

      
   // burning points 
   for (int i = 0; i < number_of_points_to_burn ; i++) 
      hmc.apply(rng, walk_len);
   // actual sampling 
   for (int i = 0; i < number_of_points ; i++) { 
   hmc.apply(rng, walk_len); 
   rand_points.push_back(hmc.x); 
   }
   return rand_points;                            
}

  
