//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Wrapper for the mhg_worker function to compute the hypergeometric
			 function of matrix argument (i.e. FLENS --> C arrays).
             
*/

#ifndef MHG_CPP
#define MHG_CPP 1


#include <flens/flens.cxx>

#include "mhg_worker.cpp"

using namespace flens;


template<typename VEC>
typename VEC::ElementType
mhg(const int MAX, const double alpha, const VEC &p, const VEC &q, const VEC &x) {

	typename VEC::ElementType ans = mhg_worker(MAX, alpha, 
												p.length(), p.data(), 
												q.length(), q.data(), 
												x.length(), x.data() );
	
	return ans;
		
}          

#endif	//MHG_CPP