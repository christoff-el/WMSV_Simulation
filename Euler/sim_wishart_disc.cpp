//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef SIM_WISHART_DISC_CPP
#define SIM_WISHART_DISC_CPP 1

#include "sim_wishart_disc.h"


//Constructor:
Sim_Wishart_Disc::Sim_Wishart_Disc(Rng_Machinery &_rng, double _T, double _h, int _n, double _delta,
					 			   matrix &_sigma, matrix &_H, matrix &_x_0, double seed_adj)
		:
		    rng(_rng),
			T(_T),
			h(_h),
			h_orig(_h),
			n(_n),
			delta(_delta),
			sigma(_sigma),
			H(_H),
			x_0(_x_0),
			timestep(1,0),
			X(1,_x_0),
			W(_rng,n),
			simulated(false),
			ref_non_conv(0)
{
}

//Simulate the Wishart process using mesh refinement:
void
Sim_Wishart_Disc::simulate(int retry) {

	//If already simulated, reset before simulating:
	if (simulated) {
		
		Sim_Wishart_Disc::reset();
		
	}
	//Otherwise, set the status to be simulated:
	else {
	
		simulated = true;
		
	}
	
	//Reset to original h:
	h = h_orig;

	// drift_fix = delta*sigma^T*sigma :
	matrix drift_fix;
	drift_fix = delta * flens::transpose(sigma) * sigma;
	
	matrix vola(n,n);
	matrix drift(n,n);
	matrix X_new = x_0;
	
	//eigenvecs/eigenvals for sqrt(X):
	matrix V_new(n,n), D_new(n,n);
	eig(V_new, D_new, x_0);

	//Update t to our current time point (0):
	double t=0;
	
	//Initialise the flag that denotes a required mesh refinement:
	bool flag = false;
	
	//Initialise objects required during the loop:
	matrix X_t;
	vec eig_vec(n);
	double mineig;

	while (t+h < T) {
	
		//std::cout<<t+h<<std::endl;

		X_t = X_new;
		
		//Calculate the volatility in 2 steps, since we need sqrtm_S_t later:
		matrix sqrtm_X_t;		
		mat3mult(sqrtm_X_t, V_new, matsqrt_diag(D_new), flens::transpose(V_new));
		
		// vola = sqrt(X_t)*(B_(t+h) - B_t)*sigma :
		mat3mult(vola, sqrtm_X_t, W.increment(h), sigma);

		//Update drift: drift = H*(X_t) :
		drift = H * X_t;
		
		//Update X: X_new = X_(t+h) = X_t + vola + vola^T + drift*h + drift^T*h + drift_fix*h :
		update_X(X_new, X_t, vola, drift, drift_fix, h);
		
		//Update eigenvalues, and find minimum one (we are checking that they are all +ve):
		eig_vec = eig(V_new, D_new, X_new);
		mineig = (*std::min_element(eig_vec.data(), eig_vec.data()+eig_vec.length()));
		
		//If an eigenvalue is negative, then we need a mesh refinement:
		while (mineig < 0) {
			
			//Halve h:
			h = 0.5*h;

			//Raise the flag to say we need a refinement:
			flag = true;
			
			//If we have refined up to the limit, give an error and stop:
			if (h < EPS) {
				std::cout<<"Error: step size converges to zero!"<<std::endl;
				
				//Increment counter:
				ref_non_conv++;
				
				//Retry if desired:
				if (retry > 0 || retry < 0) {
					(*this).simulate(retry-1);
				}
				
				return;	
						
			}

			//Update volatility and X, using pre-computed sqrtm_X_t:
			mat3mult(vola, sqrtm_X_t, W.increment(h,true), sigma);
			
			//drift = H * X_t;
			update_X(X_new, X_t, vola, drift, drift_fix, h);
				
			//Update eigenvalues and find the minimum one:	
			eig_vec = eig(V_new, D_new, X_new);
			mineig = (*std::min_element(eig_vec.data(), eig_vec.data()+eig_vec.length()));
		
		}

		X.push_back(X_new);
		t += h;
		timestep.push_back(t);
		
		if (flag) {
			
			flag = false;	//Reset flag
			h = h_orig;		//This resets h to its original value after requiring a refinement.		
		
		}
		
	}
	
	double h_end = T-t;
	h = h_end;
	
	while (h_end > 0) {		//If there is a little bit of time left between the time we
							// reached and T (i.e. <h, due to some refinements)
		
		X_t = X_new;
		
		//Calculate the volatility in 2 steps, since we need sqrtm_S_t later:
		matrix sqrtm_X_t;		
		mat3mult(sqrtm_X_t, V_new, matsqrt_diag(D_new), flens::transpose(V_new));
		
		// vola = sqrt(X_t)*(B_(t+h) - B_t)*sigma :
		mat3mult(vola, sqrtm_X_t, W.increment(h), sigma);
		
		drift = H * X_t;
		update_X(X_new, X_t, vola, drift, drift_fix, h);
		
		//Update eigenvalues, and find minimum one (we are checking that they are all +ve):
		eig_vec = eig(V_new, D_new, X_new);
		mineig = (*std::min_element(eig_vec.data(), eig_vec.data()+eig_vec.length()));
		
		//If an eigenvalue is negative, then we need a mesh refinement:
		while (mineig < 0) {
			
			//Halve h:
			h = 0.5*h;
			
			//If we have refined up to the limit, give an error and stop:
			if (h < EPS) {
				std::cout<<"Error: step size converges to zero!"<<std::endl;
				
				ref_non_conv++;
				
				//Retry if desired:
				if (retry > 0 || retry < 0) {
					(*this).simulate(retry-1);
				}
				
				return;
			}

			//Update volatility and X, using pre-computed sqrtm_X_t:
			mat3mult(vola, sqrtm_X_t, W.increment(h,true), sigma);
			update_X(X_new, X_t, vola, drift, drift_fix, h);
				
			//Update eigenvalues and find the minimum one:	
			eig_vec = eig(V_new, D_new, X_new);
			mineig = (*std::min_element(eig_vec.data(), eig_vec.data()+eig_vec.length()));
		
		}

		X.push_back(X_new);
		t += h;
		timestep.push_back(t);
		
		h_end = T-t;
		h = h_end;
		
	}


}


//Simulate the Wishart process using truncation:
void
Sim_Wishart_Disc::simulate2() {

	//If already simulated, reset before simulating:
	if (simulated) {
		
		Sim_Wishart_Disc::reset();
		
	}
	//Otherwise, set the status to be simulated:
	else {
	
		simulated = true;
		
	}
	
	//Reset to original h:
	h = h_orig;

	// drift_fix = delta*sigma^T*sigma :
	matrix drift_fix;
	drift_fix = delta * flens::transpose(sigma) * sigma;
	
	matrix vola(n,n);
	matrix drift(n,n);
	matrix X_new = x_0;
	
	//eigenvecs/eigenvals for sqrt(X):
	matrix V_new(n,n), D_new(n,n);
	eig(V_new, D_new, x_0);

	//Update t to our current time point (0):
	double t=0;
	
	//Initialise objects required during the loop:
	matrix X_t;

	for (int i=0; i<T/h; ++i) {

		X_t = X_new;
		
		//Calculate the volatility in 2 steps, since we need sqrtm_S_t later:
		matrix sqrtm_X_t;		
		mat3mult(sqrtm_X_t, V_new, matsqrt_diag(D_new), flens::transpose(V_new));
		
		// vola = sqrt(X_t)*(B_(t+h) - B_t)*sigma :
		mat3mult(vola, sqrtm_X_t, W.increment(h), sigma);

		//Update drift: drift = H*(X_t) :
		drift = H * X_t;
		
		//Update X: X_new = X_(t+h) = X_t + vola + vola^T + drift*h + drift^T*h + drift_fix*h :
		update_X(X_new, X_t, vola, drift, drift_fix, h);
		
		//Update the eigendecomposition:
		eig(V_new, D_new, X_new);
		
		//Keep only the 'positive part:
		X.push_back(pos_part(X_new, V_new, D_new));
		
		t += h;
		timestep.push_back(t);
		
	}
	
	//Last bit due to funny discretisation:
	if (T-t>1e-4) {
	
		h = T-t;
		X_t = X_new;
		
		//Calculate the volatility in 2 steps, since we need sqrtm_S_t later:
		matrix sqrtm_X_t;		
		mat3mult(sqrtm_X_t, V_new, matsqrt_diag(D_new), flens::transpose(V_new));
		
		// vola = sqrt(X_t)*(B_(t+h) - B_t)*sigma :
		mat3mult(vola, sqrtm_X_t, W.increment(h), sigma);

		//Update drift: drift = H*(X_t) :
		drift = H * X_t;
		
		//Update X: X_new = X_(t+h) = X_t + vola + vola^T + drift*h + drift^T*h + drift_fix*h :
		update_X(X_new, X_t, vola, drift, drift_fix, h);
		
		//Update the eigendecomposition:
		eig(V_new, D_new, X_new);
		
		//Keep only the 'positive part:
		X.push_back(pos_part(X_new, V_new, D_new));
		
		t += h;
		timestep.push_back(t);
		
	}

}

//Reset the process:
void
Sim_Wishart_Disc::reset() {

	if (simulated) {
		//vector<double>().swap(timestep);	//Clear timestep
		timestep.clear();
		timestep.push_back(0);
		
		//vecmat().swap(X);					//Clear X
		X.clear();
		X.push_back(x_0);
		
		W.reset();							//Reset the SBM
		
	}
	
}

//Return whether process has been simulated yet:
bool
Sim_Wishart_Disc::is_simulated() {

	bool status;
	
	if (simulated) {
		status = true;
	}
	else {
		status = false;
	}
	
	return status;
	
}

//Return timestep vector:
vector<double>
Sim_Wishart_Disc::get_timestep() {

	return timestep;
	
}

//Return i-th entry in X:
matrix
Sim_Wishart_Disc::get_X_i(int i) {

	return X[i];
	
}

//Apply the regularly used X_new = X_t + V + V' + (drift+drift'+drift_fix)*h
void
Sim_Wishart_Disc::update_X(matrix &X_new, const matrix &X, const matrix &vola, 
							const matrix &drift, const matrix &drift_fix, const double h) {

	X_new  = X + vola + flens::transpose(vola) + drift*h + drift_fix*h + flens::transpose(drift*h);

}


//Return i-th W increment:
matrix
Sim_Wishart_Disc::get_W_incr_i(int i) {

	return W.get_incr_i(i);
	
}


//Return non-conversion count of refinement method:
int
Sim_Wishart_Disc::get_non_conv_count() {

	return ref_non_conv;
	
}

//Write to file:
void
Sim_Wishart_Disc::write(std::string filename) {

	int n = X[0].numRows();
	int m = X[0].numCols();
	int ts = X.size();
	
	//Check all dimensions are OK:
	if ((int)timestep.size() != ts) {
		std::cout << "Error: timestep vector does not match process dimensions." << std::endl;
		return;
	}
	
	//File handlers:
	for (int i=1; i<=n; ++i) {
		for (int j=1; j<=m; ++j) {
		
			std::ostringstream oss;
			oss << filename << "_" << i << "_" << j << ".txt";
			
			std::fstream f;
			f.open(oss.str().c_str(), std::ios::out);
			
			//Check file opened properly:
			if (f.is_open() == 0) {
				std::cout << "Error: couldn't open output files." << std::endl;
				return;
			}
			
			for (int k=0; k<ts; ++k) {
				f << timestep[k] << " " << X[k](i,j) << std::endl;
			}
			
			f.close();
			
		}
	}

}


#endif	// SIM_WISHART_DISC_CPP