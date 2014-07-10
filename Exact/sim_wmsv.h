//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Simulate exactly a batch of L log-asset prices
			 under the WMSV model.
             
*/

#ifndef SIM_WMSV_H
#define SIM_WMSV_H 1

#include <thread>

#include <flens/flens.cxx>

#include "../utils/utils.h"
#include "../Machinery/ode_machine.cpp"
#include "../Machinery/rng_machinery.cpp"
#include "sim_wis.cpp"


class Sim_WMSV {

	//Sub-class to handle POSIX threading:
	class PSX_wrapper {
	private:
	
		//Typedefs:
		typedef flens::DenseVector<flens::Array<double> >   fvec;
		
		//Reference to parent class:
		Sim_WMSV &parent;
		
		//Necessary thread-specific parameters:
		int  step;
		int  start, startN, startL;
		int  count, countN, countL;
		double h_in;
		
		//Output handling:
		fvec *params;
	
	public:

		//Constructors:
		PSX_wrapper(Sim_WMSV &_parent, int _step, int _start, int _count);
		PSX_wrapper(Sim_WMSV &_parent, int _step, int _start, int _count, fvec *_params);
		PSX_wrapper(Sim_WMSV &_parent, int _step, int _startN, int _countN,
						int _startL, int _countL, double _h_in);
	
		//To activate POSIX threads (upon initialisation):
		void
		operator()();

	};
	
private:

	//Define pi:
	const double PI = 3.14159265358979323846264338;
	
	//Typedefs:
	typedef complex<double>									z;
	typedef flens::GeMatrix<flens::FullStorage<double> > 	matrix;
	typedef flens::GeMatrix<flens::FullStorage<z> >  		zmatrix;
	typedef flens::DenseVector<flens::Array<int> >          fivec;
	typedef flens::DenseVector<flens::Array<double> > 		fvec;
	typedef flens::DenseVector<flens::Array<z> > 			zfvec;
	typedef vector<z>			zvec;
	typedef vector<matrix> 		vecmat;
	typedef vector<zmatrix>		zvecmat;
	typedef matrix::ConstView 	ConstView;

	const flens::Underscore<matrix::IndexType>   _;

	//Rng machinery:
	Rng_Machinery &rng;
	
	//Wishart dimension:
	int d;
	
	//Precision parameters:
	double epsi;
	int M;
	double tol;
	int maxit;
	double h_der;
	
	//Model parameters:
	matrix 	x;			//Initial volatility
	double  y;			//Initial log-asset price
	double 	delta;		//Wishart delta
	matrix 	H;			//Wishart H
	matrix 	Sigma;		//Wishart Sigma
	matrix  R;			//Wishart R
	double  r;			//Risk-free drift
	double	t;			//Simulation end time
	
	//Volatility simulator object and storage:
	Sim_Wis sim_wis;
	vecmat  xT;
	
	//Means/sigs of simulations storage:
	fvec mus, sigs;
	
	//ODE machinery object:
	Ode_machine ode;
	
	//ODE matrices:
	zmatrix psi_00, Psi_00, V_00, V_00_inv, psi_0h, Psi_0h, V_0h, V_0h_inv;
	zvecmat psi_0u, Psi_0u, V_0u, V_0u_inv;
	z		phi_00, phi_0h;
	zvec	phi_0u;
	
	//Characteristic function storage:
	zmatrix charfs;
	
	//Characteristic function precomputations:
	zfvec   mhg_p, mhg_q;
	z 		char_exppref_h;
	zmatrix char_hgpref_0, char_hgpref_h;
	zvec    char_exppref_u;
	zvecmat char_hgpref_u;
	zvec    char_hgpref_0_xT_eig_HG;
	
	//Distribution evaluation precomputations:
	fvec    F_sin_l_eps, F_cos_l_eps;
	
	//Derived accuracy parameters:
	double l_eps, h;
	int N;
	
	//Simulation storage:
	fvec Y;
	
	
	//---Private functions---//
	
	//Simulation steps:
	void step0(int start, int count);
	void step0MP(int start, int count);
	fvec step1(int start, int count);
	void step2(int start, int count);
	void step2ad(int start, int count, const fvec &ns);
	void step2MP(int start, int count);
	void step3(int start, int count);
	void step5(int start, int count);
	void step5ad(int start, int count, const fvec &ns);
	
	/*fvec
	grid_gen(const double l_eps, const double h, const double N);

	double
	nr(double lp, double up, double mu, double sig, double l_eps, double v);*/
	
	//Precomputation for eval_F:
	void
	F_precomp(int N, double h);
	
	//Evaluate distribution by Fourier series:
	double
	eval_F(double v, double h, int N, int l);

	//Evaluate derivative of distribution:
	double
	eval_Newton(double v, double h, int N, int l, double U);
	
	//Precomputation for characteristic function:
	void
	char_precomp(const int startN, const int countN, 
					const int startL, const int countL, const double h_in);

	//Characteristic function:
	z
	eval_char(const int inx_N, const int inx_L, const matrix &xT);
	
	//Characteristic function for derivative:
	z
	eval_char(const matrix &xT);
	
	
	
public:

	//Constructor:
	Sim_WMSV(Rng_Machinery &_rng, double _epsi, int _M, 
				matrix _x, double _y, double _delta, 
					matrix _H, matrix _Sigma, matrix _R,
						double _r, double _t, double _tol, 
							double _h_der);
	
	//Simulate:
	void
	simulate(const int L);
	
	//Simulate using multiple threads:
	void
	simulateMP(int L, const int MP);

	//Access function for array of simulates log-prices:
	fvec
	get_Y();
	
	//Compute and return Monte Carlo option price:
	double
	get_OP_call(const double K);
	
	//Access function for simulation standard error:
	double
	get_sim_sderr();
	
	//Access function for confidence interval:
	double
	get_CL(const double true_val, const double sim_val, const double sderr);
	
	//Write simulations to file:
	void
	write_sims(std::string);
	
};

#endif	//SIM_WMSV_H
