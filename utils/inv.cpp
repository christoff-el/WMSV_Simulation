//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: wrapper to invert a matrix using the trf (LU decomposition)
			 function of LAPACK. Returns inverse of A, stored in A
			 (i.e., A is mangled).
             
*/

#ifndef INV_CPP
#define INV_CPP 1

#include <flens/flens.cxx> 


template<typename MATRIX>
void
inv(MATRIX &A) {

	//Pivot object:
	flens::DenseVector<flens::Array<int> > piv(A.numRows());
	
	flens::lapack::trf(A, piv);
	flens::lapack::tri(A, piv);

}


#endif	//INV_CPP