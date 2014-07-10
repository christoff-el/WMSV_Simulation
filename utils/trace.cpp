//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Compute the trace of a matrix M.

*/

#ifndef UTILS_TRACE_CPP
#define UTILS_TRACE_CPP 1


template<typename MATRIX>
typename MATRIX::ElementType
trace(const MATRIX &M) {

	typename MATRIX::ElementType tr = 0;

	for (int i=1; i<=M.numRows(); ++i) {
	
		tr += M(i,i);
		
	}
	
	return tr;
	
}


#endif	//UTILS_TRACE_CPP