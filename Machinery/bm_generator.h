//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Class to simulate and store Brownian motions.
			 Includes ability to conditionally simulate sub-points in 
			 a chain.
             
*/

#ifndef BM_GENERATOR_H
#define BM_GENERATOR_H 1

#include <vector>

#include <flens/flens.cxx>

#include "rng_machinery.cpp"


class BM_Generator {
private:

	//Typedefs:
	typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > 	matrix;
	typedef flens::DenseVector<flens::Array<double> > 						vec;
	typedef std::vector<matrix> 											vecmat;
	typedef std::vector<vec> 												vecvec;

	//RNG Machinery
	Rng_Machinery &rng;
	
	//Inputs:
	const int n;				// <- each BM is n x n.
	
	//Storage:
	double t_now;
	vecmat W_inc; 				// <- stores increments of the BM. I.e. W_inc[i] contains B_i-B_(i-1) (W_inc[0]==0)
	std::vector<double> t;		// <- stores time points
	matrix B_next;				// <- for when we have done a mesh refinement and need to keep the next increment for later
	double t_next;				// <- time point of next BM point
	
	bool refined;				// <- flag to say whether a refinement has been performed
	
	//Private functions:
	matrix
	normrnd(int n, int m);
	
public:

	//Constructor:
	BM_Generator(Rng_Machinery &_rng, int n, double seed_adj = 0);
	
	//Simulate the next increment:
	matrix
	increment(double h, bool half_step = 0);
	
	//Return a particular increment:
	matrix
	get_incr_i(int i);
	
	//Reset the simulation:
	void
	reset();

};


#endif	//BM_GENERATOR