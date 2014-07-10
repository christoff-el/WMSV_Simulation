//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Compute trace[A*B] efficiently.
             
*/

#ifndef TRACEMUL_CPP
#define TRACEMUL_CPP 1


template<typename MATRIX>
typename MATRIX::ElementType
tracemul(const MATRIX &A, const MATRIX &B) {

	typename MATRIX::ElementType tr = 0;

	for (int i=1; i<=A.numRows(); ++i) {
		for (int j=1; j<=B.numCols(); ++j) {
		
			tr += A(i,j)*B(j,i);
			
		}		
	}
	
	return tr;
	
}


#endif	//UTILS_TRACEMUL_CPP