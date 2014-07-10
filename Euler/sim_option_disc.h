//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Prices a European call option using the discrete
			 Euler-Maruyama simulation method, by successively simulating
			 asset prices using the provided price generator.
             
*/

#ifndef SIM_OPTION_DISC_H
#define SIM_OPTION_DISC_H 1

#include <math.h>
#include <thread>
#include <sstream>
#include <string.h>

#include <flens/flens.cxx>

#include "sim_price_disc.cpp"

using namespace std;


class Sim_Option_Disc {
private:

	//Typedefs:
	typedef flens::DenseVector<flens::Array<double> > 		fvec;

	//Inputs:
	const int N;
	const double T;
	const double r;
	const double strike;
	
	//Whether simulation has been run:
	bool simulated;
	
	//Price process generator (supplied):
	Sim_Price_Disc PPG;
	
	//Storage of values:
	double *values;
	
	//Storage of prices:
	fvec prices;
	
	//Output:
	double o_val;
	
	
public:

	//Constructor:
	Sim_Option_Disc(int _N, double _T, double _r, double _strike, Sim_Price_Disc _PPG);

	//Simulate:
	void
	simulate();
	
	//Simulate using multiple threads:
	void
	simulate_psx(int num_threads);
	
	//Write simulated prices to file:
	void
	write(std::string filename);
	
	//Access function to get option price:
	double
	get_value(bool sim_if_necessary=0, int num_threads=-1);
	
	//Access function for simulation standard error:
	double 
	get_sim_sderr();
	
	//Access function for simulation confidence interval:
	double 
	get_CL(const double true_val, const double sim_val, const double sderr);

};


//-----------POSIX Wrapper-----------//

//POSIX wrapper for simulation using multiple threads:
class Sim_Option_Disc_P {
private:
	
	//Thread id:
	int thread_id;
	
	//Option valuator (supplied):
	Sim_Option_Disc OV;
	
	//Output storage:
	double *value;
	double *sderr;
	double *CL;
	
	
public:

	//Constructor:
	Sim_Option_Disc_P(int _thread_id, double *_value, double *_sderr, double *_CL, Sim_Option_Disc _OV);
	
	//So that simulation starts directly after initialisation:
	void
	operator()();

};


#endif	//SIM_OPTION_DISC_H