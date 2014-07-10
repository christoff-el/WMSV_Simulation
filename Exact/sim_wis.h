//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Simulate exactly a d-dimensional Wishart process with law
             X_t ~ WIS_d(x,delta,H,Sigma;t)
             
             (See Section 3,6).
             
*/

#ifndef SIM_WIS_H
#define SIM_WIS_H 1

#include <flens/flens.cxx>

#include "../utils/utils.h"
#include "../Machinery/rng_machinery.cpp"
#include "sim_bes.cpp"


class Sim_Wis {
private:

	//Typedefs:
	typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > matrix;
	typedef flens::DenseVector<flens::Array<double> > 					  vec;
	typedef vector<matrix> 												  vecmat;
	typedef vector<vec> 												  vecvec;
	typedef matrix::ConstView 											  ConstView;
	
	const flens::Underscore<matrix::IndexType>   _;
	
	//Rng machinery:
	Rng_Machinery &rng;
	
	//Parameters:
	int 	d;				//Dimension
	matrix 	x;				//Initial volatility
	double 	delta;			//Delta
	matrix 	H;				//H
	matrix 	Sigma;			//Sigma
	double	t;				//Simulation end time
	
	//Square Bessel simulation object:
	Sim_Bes sim_bes;
	
	//Derived calculations:
	int		rank_q;
	double	sqrt_t;
	matrix	y;				//y = (theta_t \ m_t) * x * m_t.' * (inv(theta_t)).'
	matrix	theta_t;		//theta_t = inv(p) \ [cn zeros(d-n,d-n); kn eye(d-n,d-n)] (from ext_chol of qt/t)
	vecmat  sim_p;
	
	//Permanently stored matrices:
	matrix	alg1_perm;
	matrix	alg1_x_tilde;
	matrix 	alg1_cr;
	matrix 	alg1_kr;
	matrix	alg1_cr_inv;
	matrix 	alg1_tmp1;
	matrix	alg1_tmp2;
	matrix	sim_pYp;
	matrix	sim_pyp;
	matrix	sim_Y;
	
	//To perform the precomputations:
	void
	precomp();
	
	//Sub-algorithm:
	void
	alg1(matrix &Y, const matrix &x);
	
public:

	//Constructor:
	Sim_Wis(Rng_Machinery &_rng, matrix _x, double _delta, matrix _H, matrix _Sigma, double _t);
	
	//Copy constructor:
	Sim_Wis(const Sim_Wis &rhs);
	
	//Copy constructor:
	Sim_Wis(const Sim_Wis &rhs, bool null);

	//Simulate:
	matrix
	simulate();
	
	//Reinitialise object with new parameters (updates precomputations):
	void
	reinitialise(matrix _x, double _delta, matrix _H, matrix _Sigma, double _t);
	
};


#endif	//SIM_WIS_H