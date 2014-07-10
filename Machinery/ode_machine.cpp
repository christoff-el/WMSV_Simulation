//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


#ifndef ODE_MACHINE_CPP
#define ODE_MACHINE_CPP 1

#include "ode_machine.h"

//Prototype for integrand function:
int
func(double, const double*, double*, void *);

//Parameter storage:
struct
Ode_machine::Params{

  	//Params:
  	int d;
  	zmatrix sig_tilde;
  	zmatrix H;
  	zmatrix Htr;
  	zmatrix Q;
	z      u;		//<-- complex
	double r;
	double delta;  	
	
	//Temps:
	zmatrix tmp;
	zmatrix psi;
	z		phi;
	zmatrix Psi;
	zmatrix	V;
	zmatrix dpsi;
	z		dphi;
	zmatrix dPsi;
	zmatrix dV;
	
	//Indexes:
	int p_z;
	int p_psi;
	int p_phi;
	int p_Psi;
	int p_V;
	
};

//Constructor:
Ode_machine::Ode_machine(int d, const zmatrix &sigma, const zmatrix &H, const zmatrix &R,
							double r, double delta, double _T) 
		:
			sys_size(6*d*d+2),
			T(_T)
{
	
	//Initialise parameter storage:					
	params = new Ode_machine::Params;
	
	//Package parameters into storage:
	(*params).d 		= d;
	(*params).sig_tilde = flens::transpose(sigma) * sigma;
	(*params).H 		= H;
	(*params).Q			= R * sigma;
	(*params).r			= r;
	(*params).delta		= delta;
	
	//Temporary initialisation:
	(*params).tmp  = zmatrix(d,d);
	(*params).psi  = zmatrix(d,d);
	(*params).phi  = z(0,0);
	(*params).Psi  = zmatrix(d,d);
	(*params).V    = zmatrix(d,d);
	(*params).dpsi = zmatrix(d,d);
	(*params).dphi = z(0,0);
	(*params).dPsi = zmatrix(d,d);
	(*params).dV   = zmatrix(d,d);
	
	(*params).p_z = 3*d*d+1;
	(*params).p_psi = -1;
	(*params).p_phi = d*d;
	(*params).p_Psi = d*d;
	(*params).p_V = 2*d*d;
	
	//Initialise y:
	y = new double[sys_size];
	
	//Set up ODE driver:
	gsl_odeiv2_system _sys = {func, NULL, sys_size, params};
	sys = _sys;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
				  				           1e-12, 1e-12, 1e-12);
	
}

//Solve ODE:
void
Ode_machine::solve(z u) {

	//Need to reset first:
	Ode_machine::reset_y();
	
	//Update parameter storage with the supplied value of u:
	(*params).u = u;
	
	//Activate the GSL solver:
	double t0=0;
	gsl_odeiv2_driver_apply(driver, &t0, T, y);
	
}
	
//Reset y to initial conditions:
void
Ode_machine::reset_y() {

	//See Section 6.1 for details:
	int d=(*params).d;
	for (unsigned int i=0; i<sys_size; ++i) {
		y[i]=0;
	}
	for (int i=(d*d+1); i<(2*d*d+1); i+=(d+1)) {
		y[i]=1;
	}
	
}

//Access function for psi:
Ode_machine::zmatrix
Ode_machine::get_psi() {

	int d=(*params).d;
	for (int i=1; i<=d; ++i) {
		for (int j=1; j<=d; ++j) {
			(*params).psi(i,j) = z(y[(*params).p_psi + (i-1)*d+j] , y[(*params).p_z+(*params).p_psi + (i-1)*d+j]);
		}
	}
	return (*params).psi;
}

//Access function for phi:
Ode_machine::z
Ode_machine::get_phi() {
	return z(y[(*params).p_phi] , y[(*params).p_z+(*params).p_phi]);
}

//Access function for Psi:
Ode_machine::zmatrix
Ode_machine::get_Psi() {

	int d=(*params).d;
	for (int i=1; i<=d; ++i) {
		for (int j=1; j<=d; ++j) {
			(*params).Psi(i,j) = z(y[(*params).p_Psi + (i-1)*d+j] , y[(*params).p_z+(*params).p_Psi + (i-1)*d+j]);
		}
	}
	return (*params).Psi;
}

//Access function for V:
Ode_machine::zmatrix
Ode_machine::get_V() {

	int d=(*params).d;
	for (int i=1; i<=d; ++i) {
		for (int j=1; j<=d; ++j) {
			(*params).V(i,j) = z(y[(*params).p_V + (i-1)*d+j] , y[(*params).p_z+(*params).p_V + (i-1)*d+j]);
		}
	}
	return (*params).V;
}

//Derivative/integrand function:
int
func(double t, const double y[], double f[], void *_func_params)
{
	
	//Typedefs:
	typedef complex<double>  				z;
	
	//Collect parameters out of storage:
	Ode_machine::Params* func_params = static_cast<Ode_machine::Params*>(_func_params);
  	int d     = (*func_params).d;
  	z u       = (*func_params).u;
  	double r  = (*func_params).r;
  	double delta = (*func_params).delta;
  	
  	int p_z = 	(*func_params).p_z;
  	int p_psi = (*func_params).p_psi;
  	int p_phi = (*func_params).p_phi;
  	int p_Psi = (*func_params).p_Psi;
  	int p_V = 	(*func_params).p_V;
	
	//Form matrices from the scalar decomposition:
	for (int i=1; i<=d; ++i) {
		for (int j=1; j<=d; ++j) {
		
			(*func_params).psi(i,j) = z(y[p_psi + (i-1)*d+j] , y[p_z+p_psi + (i-1)*d+j]);
			(*func_params).Psi(i,j) = z(y[p_Psi + (i-1)*d+j] , y[p_z+p_Psi + (i-1)*d+j]);
			(*func_params).V(i,j)   = z(y[p_V   + (i-1)*d+j] , y[p_z+p_V   + (i-1)*d+j]);
			
		}
	}
	
	(*func_params).phi = z(y[p_phi] , y[p_z+p_phi]);
	
	//Compute derivatives (see Section 6.1):
	
	//dpsi = -2*psi*sig_tilde*psi + (H.'-u*Q)*psi + psi*(H-u*(Q.')) - 0.5*u*(u+1)*eye(d);
	mat3mult((*func_params).dpsi, (*func_params).psi, (*func_params).sig_tilde, (*func_params).psi);
	(*func_params).dpsi *= (-2.0);
	(*func_params).tmp = - u*(*func_params).Q;
	(*func_params).tmp += flens::transpose((*func_params).H);
	zmatrix a= (*func_params).tmp * (*func_params).psi;
	(*func_params).dpsi += a;
	
	(*func_params).tmp = (zmatrix)(- u*flens::transpose((*func_params).Q));
	(*func_params).tmp += (*func_params).H;
	(*func_params).dpsi += (*func_params).psi * (*func_params).tmp;
	(*func_params).dpsi -= eye(d,0.5*u*(u+1.0));
	
	//dphi = delta*trace(psi*sig_tilde)+u*r;
	(*func_params).tmp = (*func_params).psi * (*func_params).sig_tilde;
	(*func_params).dphi = delta*trace((*func_params).tmp) - u*r;
	
	//dPsi = (H.' - u*Q - 2*psi*sig_tilde)*Psi;
	(*func_params).tmp = -2.0 * (*func_params).psi * (*func_params).sig_tilde;
	(*func_params).tmp -= u * (*func_params).Q;
	(*func_params).tmp += (zmatrix)flens::transpose((*func_params).H);
	(*func_params).dPsi = (*func_params).tmp * (*func_params).Psi;
	
	//dV = (Psi.')*sig_tilde*Psi;
	mat3mult((*func_params).dV, (zmatrix)flens::transpose((*func_params).Psi), (*func_params).sig_tilde, (*func_params).Psi);
	
	for (int i=1; i<=d; ++i) {
		for (int j=1; j<=d; ++j) {
		
			//real
			f[p_psi + (i-1)*d+j] = (*func_params).dpsi(i,j).real();
			f[p_Psi + (i-1)*d+j] = (*func_params).dPsi(i,j).real();
			f[p_V   + (i-1)*d+j] = (*func_params).dV(i,j).real();
			
			//imag
			f[p_z+p_psi + (i-1)*d+j] = (*func_params).dpsi(i,j).imag();
			f[p_z+p_Psi + (i-1)*d+j] = (*func_params).dPsi(i,j).imag();
			f[p_z+p_V   + (i-1)*d+j] = (*func_params).dV(i,j).imag();
			
		}
	}
		
	f[p_phi]     = (*func_params).dphi.real();
	f[p_z+p_phi] = (*func_params).dphi.imag();

	//Required by GSL:
  	return GSL_SUCCESS;
  	
}


#endif	//ODE_MACHINE_CPP