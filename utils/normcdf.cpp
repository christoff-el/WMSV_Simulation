//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Computes the density and distribution functions of the standard normal
			 distribution at some point d_x, using the density formula
			 and the erf error function, respectively.
             
*/

#ifndef NORMCDF_CPP
#define NORMCDF_CPP 1

#include <cmath>
#include "erfs.h"


//CDF:
double stdNormalCDF(double d_x) {
	
	//return 0.5 * std::erfc(-d_x / sqrt(2.));	//Note: requires c++11 for erfc
	return 0.5 * erfc_sun(-d_x / sqrt(2.));		//Using Sun Microsystems code.
	
}

//CDF requires evaluation of standard normal PDF:
double stdNormalPDF(double d_x) {

	double PI = 3.14159265358979323846264338;
	
	return (1./sqrt(2.*PI))*exp(-pow(d_x,2)/2.0);

}


#endif	//NORMCDF_CPP