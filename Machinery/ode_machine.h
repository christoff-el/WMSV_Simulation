//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Framework to solve the system of ODEs decribed in
			 Section 6.
             
*/

#ifndef ODE_MACHINE_H
#define ODE_MACHINE_H 1

#include <iostream>
#include <iomanip>
#include <complex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <flens/flens.cxx>

#include "../utils/utils.h"


class Ode_machine {
public:				
	//Prototypes
	struct Params;
private:

	//Typedefs:
	typedef complex<double>  				z;
	typedef	GeMatrix<FullStorage<double> >	matrix;
	typedef GeMatrix<FullStorage<z> > 		zmatrix;

	size_t sys_size;	//Size of the scalar ODE system (6d^2 +2)
	Params *params;		//Storage of the model parameters
	double *y;			//Storage of the solution vector
	double T;			//Time to solve to
	
	//GSL solver machinery:
	gsl_odeiv2_system sys;
	gsl_odeiv2_driver *driver;
	
	//Reset the solution vector:
	void
	reset_y();
	
public:

	//Constructor:
	Ode_machine(int d, const zmatrix &sigma, const zmatrix &H, const zmatrix &R,
						double r, double delta, double _T);
	
	//Solve the initialised system with the variable u:		
	void
	solve(z u);
	
	//Access functions to acquire the solutions:
	zmatrix get_psi();
	z		get_phi();
	zmatrix get_Psi();
	zmatrix get_V();

};


#endif	//ODE_MACHINE