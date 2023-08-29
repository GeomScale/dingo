#include "ode_solvers/ode_solvers.hpp"
#include "random_walks.hpp"

template<class NT, class Point, class Polytope> std::list<Point> hmc_leapfrog_gaussian(int walk_len,
                                    int number_of_points, 
                                    int number_of_points_to_burn, 
                                    NT variance,
                                    Point starting_point,
                                    Polytope HP) {
                               
   int d = HP.dimension();
   std::list<Point> rand_points;                                 
   typedef GaussianFunctor::GradientFunctor<Point> NegativeGradientFunctor;
   typedef GaussianFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
   typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
   typedef LeapfrogODESolver<Point, NT, Polytope, NegativeGradientFunctor> Solver;
   
   unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
   RandomNumberGenerator rng(rng_seed); 

   GaussianFunctor::parameters<NT, Point> params(starting_point, 2 / (variance * variance), NT(-1));
   NegativeGradientFunctor F(params);
   NegativeLogprobFunctor f(params);
   HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, d);

   HamiltonianMonteCarloWalk::Walk<Point, Polytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>hmc(&HP, starting_point, F, f, hmc_params);

   // burning points 
   for (int i = 0; i < number_of_points_to_burn ; i++) { 
      hmc.apply(rng, walk_len);
    }

   // actual sampling 
   for (int i = 0; i < number_of_points ; i++) { 
      hmc.apply(rng, walk_len); 
      rand_points.push_back(hmc.x); 
   }   
   return rand_points;                                  
} 
 
template<class NT, class Point, class Polytope> std::list<Point> hmc_leapfrog_exponential(int walk_len,
                                    int number_of_points, 
                                    int number_of_points_to_burn, 
                                    NT variance,
                                    Point bias_vector,
                                    Point starting_point,
                                    Polytope HP) {

   int d = HP.dimension();
   std::list<Point> rand_points;                                 
   typedef ExponentialFunctor::GradientFunctor<Point> NegativeGradientFunctor;
   typedef ExponentialFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
   typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
   typedef LeapfrogODESolver<Point, NT, Polytope, NegativeGradientFunctor> Solver;
      
   unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
   RandomNumberGenerator rng(rng_seed); 

   ExponentialFunctor::parameters<NT, Point> params(bias_vector, 2 / (variance * variance));
      
   NegativeGradientFunctor F(params);
   NegativeLogprobFunctor f(params);
   HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, d);
      
   HamiltonianMonteCarloWalk::Walk<Point, Polytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>hmc(&HP, starting_point, F, f, hmc_params);

      
   // burning points 
   for (int i = 0; i < number_of_points_to_burn ; i++) { 
      hmc.apply(rng, walk_len);
   }
   // actual sampling 
   for (int i = 0; i < number_of_points ; i++) { 
      hmc.apply(rng, walk_len); 
      rand_points.push_back(hmc.x); 
   }
   return rand_points;                            
}

  
