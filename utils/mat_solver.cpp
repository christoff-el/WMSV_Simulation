//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Wrapper to solve AX=B using LAPACK.
             
*/

#ifndef MAT_SOLVER_CPP
#define MAT_SOLVER_CPP 1

#include <flens/flens.cxx>

using namespace flens;


void
mat_solver(GeMatrix<FullStorage<double> > &X, 
     const GeMatrix<FullStorage<double> > &A, 
     const GeMatrix<FullStorage<double> > &B) {

	typedef FullStorage<double, ColMajor>  FS;
    typedef GeMatrix<FS>                   GeMatrix;
	typedef GeMatrix::View                 GeMatrixView;
	
	//B gets mangled to X, so copy it:
	X = B;
	
	//A gets mangled
	GeMatrix Atemp;
	
	//Note: no checking that A,B,X all nxn.
	int n = A.numRows();
	
	//Unrequired pivot array:
	DenseVector<Array<int> >	piv(n);
	
	const Underscore<GeMatrix::IndexType>   _;
	
	for (int i=1; i<=A.numRows(); ++i) {
	
		GeMatrixView       			Xi = X(_(1,n),_(i,i));
		
		Atemp = A;
		lapack::sv(Atemp, piv, Xi);

	}

}

#endif	//MAT_SOLVER_CPP