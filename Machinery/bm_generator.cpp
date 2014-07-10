//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef BM_GENERATOR_CPP
#define BM_GENERATOR_CPP 1

#include "bm_generator.h"


//Constructor:
BM_Generator::BM_Generator(Rng_Machinery &_rng, int _n, double seed_adj)
		:
			rng(_rng),
			n(_n),
			t_now(0),
			W_inc(1,matrix(n,n)),
			t(1,0),
			refined(false)
{
}

//Simulate the next increment:
matrix
BM_Generator::increment(double h, bool half_step) {
	
	//Usual increment = sqrt(h)*Z, Z~N(0,1)
	if (!half_step) {
	
		if (refined) {
		
			//Check requested point is after stored B_next:
			if (t_now + h - t_next + EPS < 0) {
				std::cerr << "Error in BM generator: wrong h given." << std::endl;
				std::cerr << "h: " << h << ", t_now: " << t_now << ", t_next: " << t_next << std::endl;
			}
			
			//Calculate next increment using the already stored B_next as the previous one:
			W_inc.push_back(B_next + sqrt(t_now + h - t_next + EPS) * normrnd(n,n));
			
			refined = false;
					
		}
		else {
		
			W_inc.push_back(sqrt(h) * (normrnd(n,n)));
			
		}	

		t_now += h;
		t.push_back(t_now);
		

	}
	//'Inbetween' BM point = 0.5*[most recent increment] + sqrt(0.5*h)*Z, Z~N(0,1)
	// Note that this overwrites inserts the in between BM point. The old BM point is now at the back
	//  of the vector, and needs to be used as the next BM point.
	// NOTE: the h supplied to this function should be the time increment desired.
	//		e.g. if the original time increment is h, and we want a sub increment h/2, then
	//			enter h/2, not h!
	else {

		//Store old BM point as temporary:
		matrix BM_inc_old = W_inc.back();
			
		//Calculate in between point using the last increment:
		W_inc.back() = 0.5*BM_inc_old + sqrt(0.5*h) * normrnd(n,n);
	
		//Check if already refined before:
		if (refined) {
		
			//Refresh the last increment wrt the in between increment, and store it for later:
			B_next += (BM_inc_old - W_inc.back());
			
		}
		else {
		
			//Store time point of next increment:
			t_next = t_now;
		
			//Refresh the last increment wrt the in between increment, and store it for later:
			B_next = BM_inc_old - W_inc.back();
			
			refined = true;

		}
		
		//Update time:
		t_now -= h;
		t.back() = t_now;
		
	}
			
	return W_inc.back();

}

//Fill a given matrix with standard normal random numbers:
matrix
BM_Generator::normrnd(int n, int m) {

	matrix R(n,m);
	
	for (int i=R.firstRow(); i<=R.lastRow(); ++i) {
		for (int j=R.firstCol(); j<=R.lastCol(); ++j) {

			R(i,j) = rng.gen_std_norm();
			
		}
	}
	
	return R;

}

//Get the i-th increment:
matrix
BM_Generator::get_incr_i(int i) {

	return W_inc[i];
	
}

//Reset the SBM
void
BM_Generator::reset() {

	//Reset W_inc
	W_inc.clear();
	W_inc.push_back(matrix(n,n));
	
	//Reset t
	t.clear();
	t.push_back(0);
	
	t_now = 0;
	refined = false;
	
}


#endif	//BM_GENERATOR_CPP