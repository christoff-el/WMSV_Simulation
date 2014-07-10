//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Centralised generation of pseudo-random numbers using the 
			 Mersenne Twister algorithm. All random numbers used in the entire 
			 simulation program must be provided by this class. 
			 By centralising the random sequence in this way, we ensure that 
			 all random numbers used in computations are as independent and 
			 identically distributed as the PRNG allows (for example, we prevent 
			 different threads/different parts of the program generating random 
			 numbers using the same seed, thereby using dupli- cate PRNs). 
			 
			 Includes framework for the generation of uniform, standard normal and 
			 chi-squared pseudo-random numbers.
             
*/

#ifndef RNG_MACHINERY_H
#define RNG_MACHINERY_H 1

#include <vector>
#include <unordered_map>
#include <tr1/random>
#include "assert.h"

using namespace std;


//Wrapper for chi2rnd
class Chi2rnd {
private:
	
	tr1::gamma_distribution<> chi2;
	tr1::variate_generator<tr1::mt19937, tr1::gamma_distribution<> > chi2rnd;
	
public:

	//Degenerate constructor (gives assertion error):
	Chi2rnd();
	
	//Proper constructor:
	Chi2rnd(tr1::mt19937 &engine, double dof);
	
	double
	gen();
	
};

//Common RNG machinery class:
class Rng_Machinery {
private:

	//RNG Machinery:
	double 			rngSeed;
	tr1::mt19937 	engine;						//(time(0))
	
	//Uniform:
	tr1::uniform_real<double> unif;				//(0,1)
	tr1::variate_generator<tr1::mt19937, tr1::uniform_real<double> > rand;				//(engine, unif)
	
	//Standard normal:
	tr1::normal_distribution<> std_norm;		//(0,1)
	tr1::variate_generator<tr1::mt19937, tr1::normal_distribution<> > randn;		//(engine,std_norm)
	
	//Chi squared:
	unordered_map<double,Chi2rnd> 	chi_map;
	
public:

	//Default constructor (uses time as seed):
	Rng_Machinery();
	
	//Constructor with specified seed:
	Rng_Machinery(double seed);
	
	//Generate uniform random variable:
	double
	gen_unif();
	
	//Generate standard normal random variable:
	double
	gen_std_norm();
	
	//Initialise chi squared framework for a particular dof (must be done before this dof can be simulated):
	void
	init_chi_dof(double dof);
	
	//Simulate chi squared for an already initialised dof:
	double
	gen_chi2(double dof);

};


#endif	//RNG_MACHINERY_H