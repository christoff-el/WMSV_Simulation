//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef SIM_OPTION_DISC_CPP
#define SIM_OPTION_DISC_CPP 1

#include "sim_option_disc.h"


//Constructor
Sim_Option_Disc::Sim_Option_Disc(int _N, double _T, double _r, 
									double _strike, Sim_Price_Disc _PPG) 
		:
			N(_N),
			T(_T),
			r(_r),
			strike(_strike),
			simulated(false),
			PPG(_PPG),			
			prices(N)
{

}

//Simulate option price:
void
Sim_Option_Disc::simulate() {

	o_val = 0;
	
	double val_i;
	
	for (int i=0; i<N; ++i) {

		//Simulate until price found (inf retries)
		PPG.simulate(-1);		
		
		//Add simulation to price storage:
		prices(i+1) = PPG.get_end_price(true);
		
		//Compute call option payoff:
		val_i = std::max(prices(i+1) - strike,0.0);
		
		//Add payoff to running sum of payoffs:
		o_val += val_i;

	}

	//To get the Monte Carlo price, divide by the simulation count and discount:
	o_val = (o_val/(double)N) * exp(-r*T);

}

//Simulate using posix threads:
void
Sim_Option_Disc::simulate_psx(int num_threads) {

	//Check for a sensible thread number!:
	if (num_threads < 1) {
		std::cout << "Error: invalid thread count!" << std::endl;
		return;
	}
	
	//Initialise threads:
	std::thread threads[num_threads];
	
	//Distribute work evenly:
	int elsPerThread = N / num_threads;					//Intentional truncation
	int leftovers = N - elsPerThread * num_threads;

	//Initialise return packages:
	double *values = new double[num_threads];
	double *sderrs = new double[num_threads];
	double *CLs = new double[num_threads];
	
	//Spawn threads:
	for (int i=0; i<num_threads; ++i) {
	
		int elCount = (i < leftovers) ? (elsPerThread+1) : (elsPerThread);	
		
		//Use POSIX wrapper:
		threads[i] = std::thread(Sim_Option_Disc_P(i, &values[i], &sderrs[i], &CLs[i], Sim_Option_Disc(elCount,T,r,strike,PPG)));
	
	}

	//Wait for all threads to complete:
	for (int i=0; i<num_threads; ++i) {
	
		threads[i].join();
		
	}
	
	//Collect and output results:
	double value = 0;
	for (int i=0; i<num_threads; ++i) {
	
		value += values[i];
		std::cout << "Thread " << i << ": sderr=" << sderrs[i] << std::endl;
		std::cout << "Thread " << i << ": CL=" << CLs[i] << std::endl;		
		
	}
	
	//Aggregate option value:
	o_val = value / (double)num_threads;
	
}

//Write simulations to file:
void
Sim_Option_Disc::write(std::string filename) {

	std::fstream f;
	f.open(filename.c_str(), std::ios::out);
	if (f.is_open()) {
		
		for (int i=1; i<=N; ++i) {
		
			f << std::setprecision(15) << log(prices(i)) << std::endl;
			
		}
		
		f.close();
		
	}
	
}

//Return calculated price (if calculated, else -1):
double
Sim_Option_Disc::get_value(bool sim_if_necessary, int num_threads) {

	double value;
	
	//If option price has been simulated..
	if (simulated) {
	
		//Just return the simulated value:
		value = o_val;
		
	}
	//otherwise..
	else {
		
		//Simulate if we are allowed:
		if (sim_if_necessary) {
			
			if (num_threads < 2) {
				Sim_Option_Disc::simulate();
			}
			else {
				Sim_Option_Disc::simulate_psx(num_threads);
			}
			
			value = o_val;
			
		}
		//otherwise, return -1
		else {
		
			value = -1;
			
		}
	}	
		
	return value;

}

//Get simulation standard error:
double
Sim_Option_Disc::get_sim_sderr() {

	double mean = 0;
	
	for (int i=1; i<=prices.length(); ++i) {
		
		mean += (prices(i));
		
	}
	
	mean /= (double)prices.length();
	
	double sd = 0;
	for (int i=1; i<=prices.length(); ++i) {
	
		sd += pow(((prices(i))-mean),2.0);
		
	}
	
	sd /= ((double)prices.length()-1.0);
	sd = sqrt(sd);
	
	return sd / sqrt((double)prices.length());

}

//Get simulation confidence interval:
double 
Sim_Option_Disc::get_CL(const double true_val, const double sim_val, const double sderr) {

	double alpha = std::abs(true_val - sim_val) / sderr;
	return 2*(1-stdNormalCDF(alpha));
	
}


//-----------POSIX Wrapper-----------//

//Constructor:
Sim_Option_Disc_P::Sim_Option_Disc_P(int _thread_id, double *_value, double *_sderr, double *_CL, Sim_Option_Disc _OV)
		:
			thread_id(_thread_id),
			OV(_OV),
			value(_value),
			sderr(_sderr),
			CL(_CL)
{

}

//What to do when initialised as thread:
void
Sim_Option_Disc_P::operator()() {
	
	//Simulate:
	double a = OV.get_value(true);
	
	//Store the value:
	value[0]=a;
	
	//Error analysis:
	sderr[0]=OV.get_sim_sderr();
	CL[0]=OV.get_CL(0.191575,a,sderr[0]);
	
	//Write to file (if enabled):
	std::ostringstream ofn;
	ofn << "out_disc_" << thread_id << ".txt";
	//OV.write(ofn.str());
	
	//Tell user that this thread has finished:
	std::ostringstream oss;
	oss << "Thread " << thread_id << " -- finished." << std::endl;
	std::cout << oss.str();
	
}


#endif	//SIM_OPTION_DISC_CPP
