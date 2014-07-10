//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Computes the determinant of a matrix by using
			 the trf LAPACK function 
			 (i.e. LU decomposition, then take product of diagonal).
             
*/

#ifndef DET_CPP
#define DET_CPP 1

#include <flens/flens.cxx>


template<typename MATRIX>
typename MATRIX::ElementType
det(const MATRIX &A) {

	//Dimension:
	int n = A.numCols();
	
	//Since A gets mangled:
	MATRIX B=A;
	
	//Pivot thing:
	flens::DenseVector<flens::Array<int> > piv(n);
	
	//Compute LU decomposition:
	flens::lapack::trf(B, piv);
	
	//Compute determinant (product of diagonal):
	typename MATRIX::ElementType d = 1;
	for (int i=1; i<=n; ++i) {
		d *= B(i,i);
	}
	
	return d;
	
}

#endif	//DET_CPP