//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Computes the elementwise square-root of a diagonal matrix.
             
*/

#ifndef UTILS_MATSQRT_DIAG_CPP
#define UTILS_MATSQRT_DIAG_CPP 1


template<typename MATRIX>
MATRIX
matsqrt_diag(const MATRIX &R) {

	MATRIX SR(R.numRows(),R.numCols());
	
	//Only touch the diagonal (no point doing sqrt(0)!):
	for (int i=R.firstRow(); i<=R.lastRow(); ++i) {

			SR(i,i) = sqrt(R(i,i));

	}
	
	return SR;	
	
}


#endif	//UTILS_MATSQRT_DIAG_CPP