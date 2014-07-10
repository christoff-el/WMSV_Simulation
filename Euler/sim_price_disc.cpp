//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef SIM_PRICE_DISC_CPP
#define SIM_PRICE_DISC_CPP 1

#include "sim_price_disc.h"


//Constructor:
Sim_Price_Disc::Sim_Price_Disc(Rng_Machinery &_rng, double _T, double _h, int _n, double _alpha,
					 			matrix &_Q, matrix &_K, matrix &_x_0, double _r, 
					    	       matrix &_R, double _y_0)
		:
		    rng(_rng),
			T(_T),
			h(_h),
			n(_n),
			alpha(_alpha),
			Q(_Q),
			K(_K),
			x_0(_x_0),
			r(_r),
			R(_R),
			y_0(_y_0),
			X(_rng,_T,_h,_n,_alpha,_Q,_K,_x_0),
			Y(1,y_0),
			Z(_rng,n),
			simulated(false)
{
}

//Simulate the volatility process:
void
Sim_Price_Disc::simulate_volatility(int retries) {

	X.simulate(retries);
	
}

//Simulate the price process:
void
Sim_Price_Disc::simulate(int retries) {

	//If we are re-simulating, we need to reset:
	if (simulated) {
	
		//X.simulate(retries);			//Re-simulate X
		X.simulate2();
		
		//Clear Y
		Y.clear();
		Y.push_back(y_0);
		
		//Reset the SBM Z:
		Z.reset();						
	
	}
	else {
	
		simulated = true;
		
		//Simulate volatility if necessary:
		if (X.is_simulated() == false) {
	
			//X.simulate(retries);
			X.simulate2();
		
		}
		
	}

	//Get the timestep vector from the volatility:
	vector<double> timestep = X.get_timestep();
	
	double h_t, drift;
	matrix vola(n,n), tmp(n,n);
	
	//For each discretisation step:
	for (unsigned int i=1; i<timestep.size(); ++i) {
	
		//Compute stepsize of the Brownian motion (in case this has changed):
		h_t = timestep[i] - timestep[i-1];

		//Even though we are not interested in eigenvalues anymore, we need the decomposition
		// for the matrix square root.
		matrix V_x(n,n), D_x(n,n);
		eig(V_x, D_x, X.get_X_i(i-1));
		
		//Compute next simulation point using the discretisation formula (see Section 5)
		matrix V_r(n,n), D_r(n,n);
		eig(V_r, D_r, (eye(n) - R * flens::transpose(R)));

		drift = (r - 0.5*trace(X.get_X_i(i-1))) * h_t;

		mat4mult(vola, Z.increment(h_t), V_r, matsqrt_diag(D_r), flens::transpose(V_r));
		tmp = X.get_W_incr_i(i) * flens::transpose(R) + vola;

		mat4mult(vola, V_x, matsqrt_diag(D_x), flens::transpose(V_x), tmp);

		//Add simulation point to storage:
    	Y.push_back(Y[i-1] + drift + trace(vola));

	}
	
}	

//Return last price, if simulated, else -1:
double
Sim_Price_Disc::get_end_price(bool sim_if_necessary) {

	double price;
	
	if (simulated) {
	
		price = exp(Y.back());
		
	}
	else {
	
		if (sim_if_necessary) {
		
			Sim_Price_Disc::simulate(-1);
			price = exp(Y.back());
			
		}
		else {
					
			price = -1;
			
		}
	}

	return price;
		
}

//Return non-conversion count of Wishart object:
int
Sim_Price_Disc::get_non_conv_count() {

	return X.get_non_conv_count();
	
}

//Write to file:
void
Sim_Price_Disc::write(std::string filename) {

	int ts = Y.size();
	
	//Check all dimensions are OK:
	if ((int)(X.get_timestep()).size() != ts) {
		std::cout << "Error: timestep vector does not match process dimensions." << std::endl;
		return;
	}
	
	std::ostringstream oss;
	oss << filename << "_Y.txt";
	
	std::fstream f;
	f.open(oss.str().c_str(), std::ios::out);
	
	if (f.is_open() == 0) {
		std::cout << "Error: couldn't open output file." << std::endl;
		return;
	}
	
	vector<double> timestep = X.get_timestep();
	
	for (int i=0; i<ts; ++i) {
		f << timestep[i] << " " << Y[i] << std::endl;
	}
	
	f.close();

}


#endif	//SIM_PRICE_DISC_CPP
