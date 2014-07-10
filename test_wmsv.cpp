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
             Monte Carlo simulation via the EXACT simulation method.
             
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
#include "Exact/sim_wis.cpp"
#include "Exact/sim_wmsv.cpp"

using namespace std;
using namespace flens;

/* Compilation Commands:

OS X:
g++ -DUSE_CXXLAPACK -I ~/proj/FLENS/ -L /usr/local/Cellar/gsl/1.16/lib/ -lgsl -framework veclib -Wall -std=c++11  -DNDEBUG -O3 -o b.sim   test_wmsv.cpp -fopenmp

Linux/Pacioli:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/davis/gsl/lib/
g++  -L../../lapack-3.5.0 -L/cm/shared/apps/blas/gcc/1/lib64 -L/cm/shared/apps/mpfr/3.1.0/lib -L../../gsl/lib  -I ../../FLENS/ -I /cm/shared/apps/mpfr/3.1.0/include/ -I../../lapack-3.5.0/lapacke/include/ -I../../gsl/include/ -O3 -Wall -fpermissive -DPACIOLI -DNDEBUG -O3 -std=c++11 -o b.sim  test_wmsv.cpp -lgsl -lgslcblas -ltmglib -llapacke -llapack -lblas -lgfortran -lmpfr -march=corei7-avx -mtune=core-avx-i -mno-avx -mno-aes -pthread -fopenmp

*/


int main(int argc, char *argv[]){

	//Typedefs:
	typedef flens::GeMatrix<flens::FullStorage<double> > matrix;
	
	//Check input parameters:
    if (argc!=6) {
        cerr << "usage: " << argv[0] << "  Sim_count, threads, write, name, t" << endl;
        exit(-1);
    }
    
    //Collect input parameters:
    int L = atoi(argv[1]);					//Simulation count
    int MP = atoi(argv[2]);					//Number of threads

    bool write = atoi(argv[3]);				//Write simulations to file OnOff
    std::string outname = "out_";
    outname += argv[4];
    outname += std::string(".txt");
    
    double T = atof(argv[5]);				//Simulation end time
    
    //Fixed precision parameters (see Section 7):
    double h_der = 5e-3;
    double epsi  = 1e-4;
    double tol   = 5e-5;
    
    //Display input and precision parameters to user:
    std::cout << std::endl;
    std::cout << "L=" << L << ", threads=" << MP << ", epsi=" << epsi << 
    				", tol=" << tol << ", h_der=" << h_der << ", T="<< T << std::endl;

	//In debug mode, we set 0 to be the seed of the random number generator:
	#ifndef NDEBUG
	Rng_Machinery rng(0);
	#endif
	
	//Not in debug mode, we use the current time as the seed of the random number generator:
	#ifdef NDEBUG
	Rng_Machinery rng;
	#endif
	
	//Precisions:
	int M = 5; 
	
	//Model parameters:
	int d = 2;				//Wishart process dimension
	
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
	
	//Initialise the wmsv exact simulatino object using the model parameters:
	Sim_WMSV sim_wmsv(rng,epsi,M,x,y,delta,H,sigma,R,r,T,tol,h_der);
	
	//Simulate L simulations using MP threads:
	sim_wmsv.simulateMP(L,MP);
	
	//We are done, so stop the timer:
	timer.stop();
	
	//If writing to file is enabled, write simulations to file:
	if (write) {
		sim_wmsv.write_sims(outname);
	}
	
	//Error analysis:
	double true_val = 0.191575;								//Theoretical option price
	double sim_val = sim_wmsv.get_OP_call(strike);			//Simulation option price
	double sderr = sim_wmsv.get_sim_sderr();				//Simulation standard error
	
	//Output results:
	std::cout << setprecision(10) << "Price: " << sim_val << std::endl;
	std::cout << setprecision(10) << "StdErr: " << sderr << std::endl;
	std::cout << setprecision(10) << "CL: " << sim_wmsv.get_CL(true_val,sim_val,sderr) << std::endl;
	std::cout << "elapsed: " << timer.elapsed() << std::endl;
	
	return 0;
	
}
