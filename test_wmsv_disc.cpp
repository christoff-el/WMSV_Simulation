//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Sample main function to price a European option using
             Monte Carlo simulation via the EULER MARUYAMA 
             simulation method.
             
*/

#ifndef PACIOLI
#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif
#endif

#include <iostream>
#include <complex>

#include <flens/flens.cxx>
#define EPS 2.220446049250313e-15

#include "utils/utils.h"
#include "Machinery/rng_machinery.cpp"
#include "Euler/sim_wishart_disc.cpp"
#include "Euler/sim_price_disc.cpp"
#include "Euler/sim_option_disc.cpp"

using namespace std;
using namespace flens;

/* Compilation commands:

OS X:
g++ -DUSE_CXXLAPACK -I ~/proj/FLENS/ -L /usr/local/Cellar/gsl/1.16/lib/ -lgsl -framework veclib -Wall -std=c++11  -DNDEBUG -O3 -o b.sim_disc   test_wmsv_disc.cpp 

Linux/Pacioli:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/davis/gsl/lib/
g++  -L../../lapack-3.5.0 -L/cm/shared/apps/blas/gcc/1/lib64 -L/cm/shared/apps/mpfr/3.1.0/lib -L../../gsl/lib  -I ../../FLENS/ -I /cm/shared/apps/mpfr/3.1.0/include/ -I../../lapack-3.5.0/lapacke/include/ -I../../gsl/include/ -O3 -Wall -fpermissive -DPACIOLI -DNDEBUG -O3 -std=c++11 -o b.sim_disc  test_wmsv_disc.cpp -lgsl -lgslcblas -ltmglib -llapacke -llapack -lblas -lgfortran -lmpfr -march=corei7-avx -mtune=core-avx-i -mno-avx -mno-aes -pthread -fopenmp

*/


int main(int argc, char *argv[]){

	//Typedefs:
	typedef flens::GeMatrix<flens::FullStorage<double> > matrix;
	
	//Check input parameters:
    if (argc!=5) {
        cerr << "usage: " << argv[0] << "  Sim_count, h^-1, threads t" << endl;
        exit(-1);
    }
    
    //Collect input parameters:
    int    L      = atoi(argv[1]);		//Simulation count
    double h_inv  = atof(argv[2]);		//Number of steps per unit of time = h^-1
    int    MP     = atoi(argv[3]);		//Number of threads to use
    double T	  = atof(argv[4]);		//Simulation end time

    //Display input parameters to the user:
    std::cout << std::endl;
    std::cout << "L=" << L << ", h^-1=" << h_inv << ", threads=" << MP << std::endl;

	//In debug mode, we set 0 to be the seed of the random number generator:
	#ifndef NDEBUG
	Rng_Machinery rng(0);
	#endif
	
	//Not in debug mode, we use the current time as the seed of the random number generator:
	#ifdef NDEBUG
	Rng_Machinery rng;
	#endif
	
	//Model parameters:
	double h = 1.0/h_inv;		//Stepsize
	int d = 2;					//Wishart process dimension
	
	matrix x(d,d);
	matrix sigma(d,d);
	matrix H(d,d);
	matrix R(d,d);
	
	x     = 0.0298, 0.0119, 0.0119, 0.0108;				//Initial volatility
	sigma = 0.3417, 0.3493, 0.1848, 0.309;				//Wishart Sigma
	H     = -1.2479, -0.8985, -0.0820, -1.1433;			//Wishart H
	R     = -0.2243, -0.1244, -0.2545, -0.7230;			//Wishart R

	double delta = 3.2; //1.1;							//Wishart delta

	double r     = 0; 									//Risk-free drift
	double y     = 0; //log(6420.54); //0;				//Initial asset price
	
	double strike = 1.0;								//Option strike price (not log)
	
	//Initialise and start timer:
	Timer timer;	
	timer.start();
	
	//Initialise the discrete asset price generation object with the model parameters:
	Sim_Price_Disc p2(rng,T,h,d,delta,sigma,H,x,r,R,y);
	
	//Initialise the option pricing object with the price generator
	Sim_Option_Disc o1(L,T,r,strike,p2);
	
	//Get the value of the option from the option pricing object (this initiates simulation):
	double o_val = o1.get_value(true,MP);
	
	//Stop the clock!
	timer.stop();

	//Output the results:
	std::cout << setprecision(10) << o_val << std::endl;
	std::cout << "elapsed: " << timer.elapsed() << std::endl;
	
	return 0;
	
}
	
	
	