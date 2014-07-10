//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: return an identity matrix multiplied by a
			 (default a=1).
             
*/

#ifndef UTILS_EYE_CPP
#define UTILS_EYE_CPP 1

#include <flens/flens.cxx>


template<typename ELTYPE=double>
flens::GeMatrix<flens::FullStorage<ELTYPE> >
eye(const int n, const ELTYPE a = 1.0) {

	flens::GeMatrix<flens::FullStorage<ELTYPE> > id(n,n);
	
	for (int i=1; i<=n; ++i) {
		
		//Set directly to a (no need to do 1*a!):
		id(i,i) = a;
		
	}
	
	return id;
	
}


#endif	//UTILS_EYE_CPP