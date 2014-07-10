//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RNG_MACHINERY_CPP
#define RNG_MACHINERY_CPP 1

#include "rng_machinery.h"


//Chi2rnd degenerate constructor:
Chi2rnd::Chi2rnd()
		:
			chi2rnd(tr1::mt19937(0), tr1::gamma_distribution<>(1))
{
	#ifndef NDEBUG
	double dof=0;
	double init=-1;
	assert(dof == init);
	#endif
}

//Chi2rnd Contructor:
Chi2rnd::Chi2rnd(tr1::mt19937 &engine, double dof)
		:
			chi2(0.5*dof),
			chi2rnd(engine,chi2)
{
}

//Simulate a chi-squared from a gamma:
double
Chi2rnd::gen() {

	return 2.0*chi2rnd();

}


//Rng_Machinery default Constructor:
Rng_Machinery::Rng_Machinery() 
		:
			engine(time(NULL)),
			unif(0,1),
			rand(engine,unif),
			std_norm(0,1),
			randn(engine,std_norm)
{
}

//Rng Machinery seed constructor:
Rng_Machinery::Rng_Machinery(double seed)
		:
			engine(seed),
			unif(0,1),
			rand(engine,unif),
			std_norm(0,1),
			randn(engine,std_norm)
{
}

//Generate uniform random number:
double
Rng_Machinery::gen_unif() {

	return rand();
	
}

//Generate standard normal random number:
double
Rng_Machinery::gen_std_norm() {

	return randn();
	
}

//Initialise a desired dof for chi2 random numbers:
void
Rng_Machinery::init_chi_dof(double dof) {

	chi_map.insert({dof,Chi2rnd(engine,dof)});

}

//Generate chi2 with desired degrees of freedom:
double
Rng_Machinery::gen_chi2(double dof) {

	return chi_map[dof].gen();

}
	
#endif	//RNG_MACHINERY_CPP