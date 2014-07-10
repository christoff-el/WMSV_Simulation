//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Computes the 1-norm of a matrix A.
			 I.e. computes max(sum(A,1))
             
*/

#ifndef NORM1_CPP
#define NORM1_CPP 1

#include <flens/flens.cxx>


template<typename MATRIX>
double
norm1(const MATRIX &A) {

	int n = A.numRows();
	
	double max = 0;
	double run_sum;
	
	for (int j=1; j<=n; ++j) {
	
		run_sum = 0;
		
		for (int i=1; i<=n; ++i) {
			run_sum += A(i,j);
		}
		
		if (run_sum > max) {
			max = run_sum;
		}
		
	}
	
	return max;

}

#endif	//NORM1_CPP