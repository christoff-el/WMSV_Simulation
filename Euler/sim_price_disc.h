//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Simulates an asset price under the WMSV model
			 using the Euler-Maruyama discretisation method.
             
*/

#ifndef SIM_PRICE_DISC_H
#define SIM_PRICE_DISC_H 1

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <flens/flens.cxx>

#include "../utils/utils.h"
#include "../Machinery/bm_generator.cpp"
#include "../Machinery/rng_machinery.cpp"
#include "sim_wishart_disc.cpp"

using namespace std;


class Sim_Price_Disc {
private:

	//Typedefs:
	typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > 	matrix;
	typedef flens::DenseVector<flens::Array<double> > 						vec;
	typedef vector<matrix> 													vecmat;
	typedef vector<vec> 													vecvec;

	//RNG Machinery:
	Rng_Machinery &rng;
	
	//Inputs:
	const double T;			//Simulation end time
	const double h;			//Discretisation step size
	const int n;			//Wishart process dimension
	
	const double alpha;		//wishart parameters
	const matrix &Q;
	const matrix &K;
	const matrix &x_0;
	
	const double r;			//price process parameters
	const matrix &R;
	const double y_0;
	
	//Volatility (Wishart) process:
	Sim_Wishart_Disc X;
	
	//Outputs:
	vector<double> Y;
	
	//Brownian motion generator:
	BM_Generator Z;
	
	//Flag to say whether the simulation has already been run:
	bool simulated;

public:

	//Constructor:
	Sim_Price_Disc(Rng_Machinery &_rng, double _T, double _h, int _n, double _alpha,
					 matrix &_Q, matrix &_K, matrix &_x_0, double _r, 
					    matrix &_R, double _y_0);
	
	//Simulate the Wishart volatility process	
	void
	simulate_volatility(int retries = 0);
	
	//Simulate the WMSV process:
	void
	simulate(int retries = 0);
	
	//Access method for the asset price simulation:
	double										//<-- returns PRICE (not log price)
	get_end_price(bool sim_if_necessary = 0);
	
	//Get the number of simulations that didn't converge:
	int
	get_non_conv_count();
	
	//Write simulations to file:
	void
	write(std::string filename="out");

};


#endif	//SIM_PRICE_DISC_H