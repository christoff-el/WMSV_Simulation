//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Simulate a square Bessel process of the form 
             dX_t = (delta-r)*dt + 2*sqrt(X_t)*dW_t
             dX_0 = y
             
             at time t.
             
*/

#ifndef SIM_BES_EX_H
#define SIM_BES_EX_H 1

#include "../Machinery/rng_machinery.cpp"


class Sim_Bes {
private:

	//Rng_Machinery:
	Rng_Machinery &rng;
	
	//Delta of the wishart (fixed):
	double delta;
	
public:

	//Constructor:
	Sim_Bes(Rng_Machinery &_rng, double _delta, int d);
	
	//Copy constructor:
	Sim_Bes(const Sim_Bes &rhs, bool null);
	
	//Simulate:
	double
	simulate(int r, double y, double t);		
	//NOTE: assumes delta-r > d-1
	
};

/*Testing:

MC sims of t=1
   dX_t = 3.2dt + 2*sqrt(X_t)*dW_t
   X_0 = 1
   
i.e. Sim_Bes(rng,3.2,2).sim(0,1,1)

gives expected value about 4.19

*/

#endif	//SIM_BES_EX_H