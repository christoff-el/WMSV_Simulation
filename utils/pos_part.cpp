//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Takes the positive part of a matrix M.
			 I.e. V D+ V^T, where VDV^T is the eigen-decomposition of M,
			 and D+ means take min(0,D_ii) for all values of the diagonal matrix D.
             
*/

#ifndef POS_PART_H
#define POS_PART_H

#include <cmath>

#include <flens/flens.cxx>

#include "eig.cpp"
#include "matmult.cpp"

typedef flens::DenseVector<flens::Array<double> > vec;


template<typename MATRIX>
MATRIX
pos_part(const MATRIX &M) {

	MATRIX V, D;
	vec eig_vec;
	double mineig;
	
	//Get eigenvalues:
	eig_vec = eig(V, D, M);
	mineig = (*std::min_element(eig_vec.data(), eig_vec.data()+eig_vec.length()));
	
	//If nothing needs to be done:
	if (mineig >= 0) {
		return M;
	}
	
	//Otherwise:
	for (int i=1; i<=D.numRows(); ++i) {
		D(i,i) = std::max(D(i,i),0.0);//std::abs(D(i,i));
	}
	
	matrix out;
	mat3mult(out, V, D, flens::transpose(V));
	
	return out;

}

//As above, but eigendecomposition is supplied:
template<typename MATRIX>
MATRIX
pos_part(const MATRIX &M, const MATRIX &V, MATRIX &D) {

	//Note that D is changed here!

	//Could be improved?
	for (int i=1; i<=D.numRows(); ++i) {
		D(i,i) = std::max(D(i,i),0e-10);
	}
	
	matrix out;
	mat3mult(out, V, D, flens::transpose(V));
	
	return out;

}

#endif	//POS_PART_H