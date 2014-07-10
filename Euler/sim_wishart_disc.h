//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Simulates a Wishart process using the Euler-Maruyama
			 discretisation approximate method.
             
*/

#ifndef SIM_WISHART_DISC_H
#define SIM_WISHART_DISC_H 1

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <flens/flens.cxx>

#include "../utils/utils.h"
#include "../Machinery/bm_generator.cpp"

using namespace std;


/* Class to generate the process given by the dynamics:

	dX_t = (delta * sigma^T * sigma + H*X_t + X_t*H^T)dt
				+ sqrt(X_t)*dW_t*sigma + sigma^T*dW_t^T*sqrt(X_t)
				
	X_0 = x_0
	
*/

class Sim_Wishart_Disc {
private:
	
	//Typedefs:
	typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > 	matrix;
	typedef flens::DenseVector<flens::Array<double> > 						vec;
	typedef vector<matrix> 													vecmat;
	typedef vector<vec> 													vecvec;
	
	//RNG Machinery
	Rng_Machinery &rng;
	
	//Inputs:
	const double T;			//Simulation end time
	double h;				//Stepsize (may be altered)
	const double h_orig;	//Stepsize (original)
	const int n;			//Wishart process dimension
	
	//Model parameters:
	const double delta;	
	const matrix &sigma;
	const matrix &H;
	const matrix &x_0;

	//Outputs:
	vector<double> timestep;
	vecmat X;
	
	//Brownian motion generator/storage:
	BM_Generator W;
	
	//Flag to say whether the simulation has already been run:
	bool simulated;
	
	//Counter for non-convergence using refinement method:
	int ref_non_conv;
	
	//Private functions:
	void
	update_X(matrix &X_new, const matrix &X, const matrix &vola, 
				const matrix &drift, const matrix &drift_fix, const double h);

public:

	//Constructor:
	Sim_Wishart_Disc(Rng_Machinery &_rng, double _T, double _h, int _n, double _delta,
					 matrix &_sigma, matrix &_H, matrix &_x_0, double seed_adj=0);
		
	//Simulate (mesh refinement):
	void									//retry = #times to retry simulation after stepsize failure.
	simulate(int retry=0);					//retry < 0 means try until a solution is found.
	
	//Simulate (full truncation)
	void 
	simulate2();
	
	//Reset simulation:
	void
	reset();
	
	//Return whether process has been simulated:
	bool
	is_simulated();
	
	//Get the timestep:
	vector<double>
	get_timestep();
	
	//Get the simulation solution X indexed at i:
	matrix
	get_X_i(int i);
	
	//Get the Brownian motion increment used in simulation, indexed at i:
	matrix
	get_W_incr_i(int i);
	
	//Get the number of non-converging simulations:
	int
	get_non_conv_count();
	
	//Write simulations to file:
	void
	write(std::string filename="out");

};

#endif // SIM_WISHART_DISC_H