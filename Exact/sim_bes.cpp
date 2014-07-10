//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef SIM_BES_CPP
#define SIM_BES_CPP 1

#include "sim_bes.h"


//Constructor:
Sim_Bes::Sim_Bes(Rng_Machinery &_rng, double _delta, int d)
		:
			rng(_rng),
			delta(_delta)
{
	
	//Initialise the potentially required chi2 degrees of freedom:
	for (int r=0; r<=d; ++r) {
	
		rng.init_chi_dof(delta - (double)r);
		
	}
	
}

//Null copy constructor (does nothing):
Sim_Bes::Sim_Bes(const Sim_Bes &rhs, bool null) 
		:
			rng(rhs.rng)
{
}

//Simulate square bessel with drift coeff delta-r:
double
Sim_Bes::simulate(int r, double y, double t) {

    return t * pow(rng.gen_std_norm()+sqrt(y/t),2.0) + rng.gen_chi2(delta-(double)r-1.0);

}


















#endif	//SIM_BES_CPP