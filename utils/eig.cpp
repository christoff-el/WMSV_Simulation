//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Compute the eigenvalues and eigenvectors of a matrix.
			 Either uses FLENS, or for 2x2 matrices computes eigenvalues directly.
             
*/

#ifndef UTILS_EIG_CPP
#define UTILS_EIG_CPP 1

#ifndef EPS
#define EPS 2.220446049250313e-15
#endif

#include <cmath>

#include <flens/flens.cxx>


typedef flens::DenseVector<flens::Array<double> > vec;
typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > matrix;


//Port of matlab eig function, with eigendecomposition:
vec
eig(matrix &V, matrix &D, const matrix &A) {

	vec eigX(A.numCols());
	
	vec D_i, work(0);
	matrix V_L;

	//The input matrix (A) is overwritten with the eigenvalue matrix (D) during ev.
	D = A;
	
	//D_0i is imaginary part of eigenvals, V_0L is left eigen vectors. 
	flens::lapack::ev(false, true, D, eigX, D_i, V_L, V, work);

	//So D already set.
	
	//We need to set small values to be zero, otherwise we get BIG problems (non-spdness)
	int dim = D.numRows();
	for (int i=1; i<=dim; ++i) {
		for (int j=1; j<=dim; ++j) {
			if (std::abs(D(i,j))<EPS) {
				D(i,j)=0.0;
			}
		}
	}

	return eigX;

}


//Port of matlab eig, returning vector of eigenvalues only:
template<typename MAT>
flens::DenseVector<flens::Array<typename MAT::ElementType> >
eig(const MAT &A) {

	typedef flens::DenseVector<flens::Array<typename MAT::ElementType> > VEC;
	
	int n = A.numCols();
	VEC eigX(n);

	MAT V_L(n,n), V(n,n);
	MAT D = A;
	
	//D_0i is imaginary part of eigenvals, V_0L is left eigen vectors. 
	flens::lapack::ev(false, false, D, eigX, V_L, V);

	return eigX;

}


//Direct calculation of eigenvalues of a 2x2 matrix, without lapack overhead:
template<typename MAT>
flens::DenseVector<flens::Array<typename MAT::ElementType> >
eig2(const MAT &A) {

	typedef typename MAT::ElementType			  EL;
	typedef flens::DenseVector<flens::Array<EL> > VEC;
	
	VEC eigX(2);
	EL  apd, sq;
	
	apd = A(1,1) + A(2,2);
	sq  = sqrt(apd*apd - 4.0*(A(1,1)*A(2,2)-A(1,2)*A(2,1)));
	
	eigX(1) = 0.5 * (apd + sq);
	eigX(2) = 0.5 * (apd - sq);
	
	return eigX;
	
}

#endif	//UTILS_EIG_CPP