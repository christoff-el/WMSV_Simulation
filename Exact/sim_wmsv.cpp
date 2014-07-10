//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef SIM_WMSV_CPP
#define SIM_WMSV_CPP 1

#include "sim_wmsv.h"


//Constructor:
Sim_WMSV::Sim_WMSV(Rng_Machinery &_rng, double _epsi, int _M,
						matrix _x, double _y, double _delta, 
							matrix _H, matrix _Sigma, matrix _R,
								double _r, double _t, double _tol, double _h_der)
		:
			rng(_rng),
			d(_H.numRows()),
			epsi(_epsi),
			M(_M),
			tol(_tol),
			maxit(20),
			h_der(_h_der),
			x(_x),
			y(_y),
			delta(_delta),
			H(_H),
			Sigma(_Sigma),
			R(_R),
			r(_r),
			t(_t),
			sim_wis(rng, x, delta, H, Sigma, t),
			ode(d, Sigma, H, R, r, delta, t),
			psi_00(d,d), Psi_00(d,d), V_00(d,d), V_00_inv(d,d),
			psi_0h(d,d), Psi_0h(d,d), V_0h(d,d), V_0h_inv(d,d),
			psi_0u(1,matrix(d,d)), Psi_0u(1,matrix(d,d)), V_0u(1,matrix(d,d)), V_0u_inv(1,matrix(d,d)),
			mhg_p(0), mhg_q(1),
			char_hgpref_0(d,d), char_hgpref_h(d,d)
{

	//Perform necessary pre-computations:
	
	//We always send this to mhg as the q-vec of stuff:
	mhg_q(1) = z(0.5*delta,0);
	
	//The ODE components at u=0 and u=h_der (for the derivatives) are model-specific:
	ode.solve(z(0,0));			//<-- u=0;
	psi_00 = ode.get_psi();
	phi_00 = ode.get_phi();
	Psi_00 = ode.get_Psi();
	V_00   = ode.get_V();
	
	ode.solve(z(0,-h_der));		//<-- u=h_der;
	psi_0h = ode.get_psi();
	phi_0h = ode.get_phi();
	Psi_0h = ode.get_Psi();
	V_0h   = ode.get_V();
	
	//Obviously, so are their inverses/determinants:
	V_00_inv = V_00; inv(V_00_inv);
	V_0h_inv = V_0h; inv(V_0h_inv);
	
	//Stuff for general characteristic function:
	//char_hgpref_0:
	mat5mult(char_hgpref_0, V_00_inv, flens::transpose(Psi_00), x, Psi_00, V_00_inv);
	char_hgpref_0 *= 0.25;
	
	
	//Stuff for derivative characteristic function:
	z char_dets_h = pow( (det(V_00) / det(V_0h)) , 0.5*delta);
	
	zmatrix tmp1(d,d);
	zmatrix tmp2(d,d);
	
	//char_exppref_h:
	mat3mult(tmp1, Psi_0h, V_0h_inv, flens::transpose(Psi_0h));
	mat3mult(tmp2, Psi_00, V_00_inv, flens::transpose(Psi_00));
	tmp1 += 2.0*psi_0h;
	tmp1 -= tmp2;
	tmp2 = tmp1*(zmatrix)x;
	char_exppref_h = char_dets_h * exp( -phi_0h - z(0,-h_der)*y -0.5*trace( tmp2 ));
	
	//char_hgpref_h:
	mat5mult(char_hgpref_h, V_0h_inv, flens::transpose(Psi_0h), x, Psi_0h, V_0h_inv);
	char_hgpref_h *= 0.25;
	
	
	/* Reference guide (See Appendix 1):
		char_hgpref =  C^phi
	    char_dets =    D^phi
	    char_exppref = E^phi
	*/

}

//Simulate!:
void
Sim_WMSV::simulate(const int L) {
	
	//Simulate volatilities (Step0):
	std::cout << "Generating volatilities.." << std::endl;
	xT.resize(L);
	Sim_WMSV::step0(0,L);
	
	//Step 1:
	std::cout << "Step 1.." << std::endl;
	mus.resize(L);
	sigs.resize(L);
	
	fvec params = Sim_WMSV::step1(1,L);
	l_eps = params(1);
	h = params(2);
	N = params(3);

	//Output derived precision parameters to user:
	std::cout << setprecision(10) << "l_eps=" << l_eps << std::endl;
	std::cout << setprecision(10) <<"h=" << h << std::endl;
	std::cout << "N=" << N << std::endl;

	//Step 2:
	std::cout << "Step 2.." << std::endl;
	
	//Resize psi_0u etc to the derived N:
	psi_0u.resize(N);
	phi_0u.resize(N);
	Psi_0u.resize(N);
	V_0u.resize(N);
	V_0u_inv.resize(N);
	
	Sim_WMSV::step2(0,N);
	
	//Step 3:
	std::cout << "Step 3 Precomps.." << std::endl;
	
	//Resize characteristic function et al. storage to derived size N:
	charfs.resize(L,N);
	char_exppref_u.resize(N);
	char_hgpref_u.resize(N);
	char_hgpref_0_xT_eig_HG.resize(L);
	
	//Perform pre-computations:
	Sim_WMSV::char_precomp(0,N,0,L,h);
	
	std::cout << "Step 3.." << std::endl;
	
	//Perform Step 3:
	Sim_WMSV::step3(1,L);
	
	//Step 4 is just simulating uniform random variables - we can do this
	// on the fly in Step 5.
	
	//Step 5:
	std::cout << "Step 5.." << std::endl;
	
	//Perform F precomputations:
	Sim_WMSV::F_precomp(N,h);
	
	//Resize storage:
	Y.resize(L);
	
	//Perform Step 5:
	Sim_WMSV::step5(1,L);
    
}

//Simulate with multiple processors!:
void
Sim_WMSV::simulateMP(const int L, const int MP) {
	
	//Initialise threads and work distribution:
	std::thread threads[MP];
	
	int elsPerThreadL = L / MP;					//Intentional truncation
	int leftoversL = L - elsPerThreadL * MP;
	
	//Distribute work evenly:
	fivec countsL(MP);
	fivec startsL(MP);
	startsL(1) = 1;
	for (int i=1; i<=leftoversL; ++i) {
		countsL(i) = elsPerThreadL + 1;
	}
	for (int i=(leftoversL+1); i<=MP; ++i) {
		countsL(i) = elsPerThreadL;
	}
	for (int i=2; i<=MP; ++i) {
		startsL(i) = startsL(i-1)+countsL(i-1);
	}

	//Simulate volatilities (Step0):
	std::cout << "Generating volatilities.." << std::endl;
	xT.resize(L);
	for (int i=0; i<MP; ++i) {
		threads[i] = std::thread(Sim_WMSV::PSX_wrapper(*this,0,startsL(i+1)-1,countsL(i+1)));
	}
	for (int i=0; i<MP; ++i) {
		threads[i].join();
	}

	//Step 1:
	std::cout << "Step 1.." << std::endl;
	mus.resize(L);
	sigs.resize(L);
	l_eps = 10e10;
	h = 10e10;
	N = -1;
	
	fvec *params = new fvec[MP];
	for (int i=0; i<MP; ++i) {
		threads[i] = std::thread(Sim_WMSV::PSX_wrapper(*this,1,startsL(i+1),countsL(i+1),&params[i]));
	}
	for (int i=0; i<MP; ++i) {
		threads[i].join();
		l_eps = std::min(l_eps,(params[i])(1));
		h = std::min(h,(params[i])(2));
		N = std::max((double)N,(params[i])(3));
	}

	//Output derived precision parameters to user:
	std::cout << setprecision(10) << "l_eps=" << l_eps << std::endl;
	std::cout << setprecision(10) <<"h=" << h << std::endl;
	std::cout << "N=" << N << std::endl;

	//Step 2:
	std::cout << "Step 2.." << std::endl;
	
	//Resize psi_0u etc to N:
	psi_0u.resize(N);
	phi_0u.resize(N);
	Psi_0u.resize(N);
	V_0u.resize(N);
	V_0u_inv.resize(N);

	//Distribute the next work evenly amongst the threads:
	int elsPerThreadN = N / MP;					//Intentional truncation
	int leftoversN = N - elsPerThreadN * MP;
	fivec countsN(MP);
	fivec startsN(MP);
	startsN(1) = 1;
	for (int i=1; i<=leftoversN; ++i) {
		countsN(i) = elsPerThreadN + 1;
	}
	for (int i=(leftoversN+1); i<=MP; ++i) {
		countsN(i) = elsPerThreadN;
	}
	for (int i=2; i<=MP; ++i) {
		startsN(i) = startsN(i-1)+countsN(i-1);
	}

	for (int i=0; i<MP; ++i) {
		threads[i] = std::thread(Sim_WMSV::PSX_wrapper(*this,2,startsN(i+1)-1,countsN(i+1)));
	}
	for (int i=0; i<MP; ++i) {
		threads[i].join();
	}
		
	//Step 3:
	std::cout << "Step 3 Precomps.." << std::endl;
	charfs.resize(L,N);
	char_exppref_u.resize(N);
	char_hgpref_u.resize(N);
	char_hgpref_0_xT_eig_HG.resize(L);
	Sim_WMSV::char_precomp(0,N,0,L,h);
	
	std::cout << "Step 3.." << std::endl;
	for (int i=0; i<MP; ++i) {
		threads[i] = std::thread(Sim_WMSV::PSX_wrapper(*this,3,startsL(i+1),countsL(i+1)));
	}
	for (int i=0; i<MP; ++i) {
		threads[i].join();
	}
	
	//Step 5:
	std::cout << "Step 5.." << std::endl;
	Sim_WMSV::F_precomp(N,h);
	Y.resize(L);
	for (int i=0; i<MP; ++i) {
		threads[i] = std::thread(Sim_WMSV::PSX_wrapper(*this,5,startsL(i+1),countsL(i+1)));
	}
	for (int i=0; i<MP; ++i) {
		threads[i].join();
	}
    
}


//------------Simulation Steps------------//

//Step0
void
Sim_WMSV::step0(int start, int count) {

	//Wishart process simulation:
	for (int i=start; i<(start+count); ++i) {
		xT[i] = sim_wis.simulate();
	}
	
}

//Step0 for multithreaded:
void
Sim_WMSV::step0MP(int start, int count) {

	//Need to copy the Wishart simulation object to avoid conflicts:
	Sim_Wis s1(sim_wis);
	
	//Then simulate each random variable:
	for (int i=start; i<(start+count); ++i) {
			xT[i] = s1.simulate();
	}
	
}

//Step 1 (single- and multi-threaded):
Sim_WMSV::fvec
Sim_WMSV::step1(int start, int count) {

	fvec params(3);
	z    fx_h;
	
	//Compute mu and sigma for each Wishart simulation:
	for (int i=start; i<(start+count); ++i) {
		
		fx_h = Sim_WMSV::eval_char(xT[i-1]);

		mus(i)  = (fx_h.imag())/h_der;
		sigs(i) = -(2*(fx_h.real()-1)) / (h_der*h_der) - mus(i)*mus(i);

	}

	//Find the global l_eps, h, N in the conservative fashion:
	double l_eps_tmp = 1e10;
	double h_tmp = 1e10;
	double Ntmp = -1e10;
	double ninv   = norminv(epsi/4.0);
	double Nconst = 1 + 0.5*sqrt(trace((matrix)(flens::transpose(R)*R)))*sqrt(-2*log(0.25*PI*epsi));
	for (int i=start; i<(start+count); ++i) {
		l_eps_tmp = std::min(l_eps_tmp, mus(i) + sigs(i)*ninv);
		h_tmp     = std::min(h_tmp, (0.1*std::abs(mus(i)) + sigs(i)) * ninv);
		Ntmp      = std::max(Ntmp, Nconst / sigs(i));
	}
	h_tmp = -PI / h_tmp;
	int N_tmp = (int)(Ntmp/h_tmp) + 1;
	
	//Package parameters for return to the main simulation function:
	params(1) = l_eps_tmp;
	params(2) = h_tmp;
	params(3) = N_tmp;
	
	return params;
	
}

//Step2:
void
Sim_WMSV::step2(int start, int count) {

	//For each lambda = -inh, n=1,...,N:
	for (int i=start; i<(start+count); ++i) {
	
		//Solve the ode system for the particular lambda:
		ode.solve(z(0,-(i+1)*h));
		
		//Collect and store the solutions:
		psi_0u[i] = ode.get_psi();
		phi_0u[i] = ode.get_phi();
		Psi_0u[i] = ode.get_Psi();
		V_0u[i]   = ode.get_V();
		
		//Precompute the inverse of V_0u:
		V_0u_inv[i] = V_0u[i]; inv(V_0u_inv[i]);
		
	}
	
}

//Step2 (for multi-threaded):
void
Sim_WMSV::step2MP(int start, int count) {

	//Must copy the ODE machinery to avoid conflicts:
	Ode_machine odeMP(d, Sigma, H, R, r, delta, t);

	//For each lambda = -inh, n=1,...,N:
	for (int i=start; i<(start+count); ++i) {
	
		//Solve the ode system for the particular lambda:
		odeMP.solve(z(0,-(i+1)*h));
		
		//Collect and store the solutions:
		psi_0u[i] = odeMP.get_psi();
		phi_0u[i] = odeMP.get_phi();
		Psi_0u[i] = odeMP.get_Psi();
		V_0u[i]   = odeMP.get_V();
		
		//Precompute the inverse of V_0u:
		V_0u_inv[i] = V_0u[i]; inv(V_0u_inv[i]);
		
	}
	
}

//Step3 (for both single- and multi-threaded):
void
Sim_WMSV::step3(int start, int count) {

	//Compute and store the characteristic function for 
	// l=1,...,L and lambda=-ih,...,-iNh:
	for (int i=start; i<(start+count); ++i) {
		for (int j=1; j<=N; ++j) {
		
			charfs(i,j) = eval_char((j-1),(i-1),xT[i-1]);

		}

	}
	
}

//Step5 (for both single- and multi-threaded):
void
Sim_WMSV::step5(int start, int count) {

	double U, vold;
	int it;
	
	//For each l=1,...,L:
	for (int i=start; i<(start+count); ++i) {
	
		//Generate uniform random number:
		U = rng.gen_unif();

		//Use the mean as best guess:
		Y(i) = mus(i) + sigs(i)*norminv(U);
		
		//Newton:
		vold = -100;
		it = 0;
		while ((std::abs(Y(i)-vold) > tol) && (it < maxit)) {
		
			vold = Y(i);

			Y(i) -= eval_Newton(Y(i),h,N,i,U);
			it++;
			
		}
		
		//Bisection if no convergence:
    	if (it==maxit) {
    	
	        double va=-10;
        	double vb=3;
        	double fc = 100;
	
	        while ( (std::abs(vb-va)>tol) && (std::abs(fc)>tol) ) {
	        
	            Y(i) = 0.5*(vb+va);
    	        fc = eval_F(Y(i),h,N,i) - U;
    	        
        	    if (fc > 0) {
            	    vb = Y(i);
            	}
	            else {
    	            va = Y(i);
        	    }
        	    
        	}
        		
        }  

    }
    
}

//Perform precomputations for distribution function evaluation:
void
Sim_WMSV::F_precomp(int N, double h) {

	F_sin_l_eps.resize(N);
	F_cos_l_eps.resize(N);
	
	double lh = l_eps*h;
	for (int i=1; i<=N; ++i) {
	
		F_sin_l_eps(i) = sin(i*lh);
		F_cos_l_eps(i) = cos(i*lh);
		
	}
	
}
	
//Evaluate distribution by Fourier series:
double
Sim_WMSV::eval_F(double v, double h, int N, int l) {

	double F = (h*(v-l_eps)) / 2.0;
	
	for (int i=1; i<=N; ++i) {
	
		F += ( (charfs(l,i).real() * (sin(v*i*h) - F_sin_l_eps(i)))
					- (charfs(l,i).imag() * (cos(v*i*h) - F_cos_l_eps(i))) )
							/ (double)i;
							
	}

	F /= PI;
	
	return F;
	
}

//Evaluate F-U divided by its derivative (saves on expensive trig functions):
double
Sim_WMSV::eval_Newton(double v, double h, int N, int l, double U) {

	double F = (h*(v-l_eps)) / 2.0;
	double f = 0.5;
	
	double vh = v*h;
	double s,c;
	for (int i=1; i<=N; ++i) {
	
		s = sin(i*vh);
		c = cos(i*vh);
		
		F += ( (charfs(l,i).real() * (s - F_sin_l_eps(i)))
					- (charfs(l,i).imag() * (c - F_cos_l_eps(i))) )
							/ (double)i;
							
		f += ( charfs(l,i).real() * c ) 
				+ ( charfs(l,i).imag() * s );
				
	}
	
	F /= PI;
	f *= (h/PI);

	return (F-U)/f;
	
}

//Perform precomputations for characteristic function:
void
Sim_WMSV::char_precomp(const int startN, const int countN, 
							const int startL, const int countL, const double h_in) {

	zmatrix tmp_fixed(d,d);
	mat3mult(tmp_fixed, Psi_00, V_00_inv, flens::transpose(Psi_00));
	
	#pragma omp parallel for
	for (int i=startN; i<(startN+countN); ++i) {
	
		z       char_dets_u;
		zmatrix tmp1, tmp2(d,d), V0u_invPsi0u(d,d);
		
		char_dets_u = pow( (det(V_00) / det(V_0u[i])) , 0.5*delta);

		V0u_invPsi0u = V_0u_inv[i] * flens::transpose(Psi_0u[i]);

		//char_exppref_u:
		tmp1 = Psi_0u[i] * V0u_invPsi0u;
		tmp1 += 2.0*psi_0u[i];
		tmp1 -= tmp_fixed;
		tmp2 = tmp1*(zmatrix)x;
		char_exppref_u[i] = char_dets_u * exp( -phi_0u[i] - z(0,-(i+1)*h_in)*y -0.5*trace( tmp2 ));
	
		//char_hgpref_u:
		mat4mult(char_hgpref_u[i], V0u_invPsi0u, x, Psi_0u[i], V_0u_inv[i]);
		char_hgpref_u[i] *= 0.25;
		
		//Set V_0u_inv to V_0u_inv-V_00_inv, since it is only needed as this from now on:
		V_0u_inv[i] -= V_00_inv;
		
	}
	
	//Eigenvalues of char_hgpref_0 * xT:
	#pragma omp parallel for
	for (int i=startL; i<(startL+countL); ++i) {
	
		zmatrix tmp1 = char_hgpref_0 * xT[i];
		char_hgpref_0_xT_eig_HG[i] = mhg(M,1.0,mhg_p, mhg_q, eig2(tmp1));
		
	}

}	

//Evaluate characteristic function using precomputed stuff:
Sim_WMSV::z
Sim_WMSV::eval_char(const int inx_N, const int inx_L, const matrix &xT) {

	zmatrix tmp1 = char_hgpref_u[inx_N] * xT;
	
	z h1 = mhg(M, 1.0, mhg_p, mhg_q, eig2(tmp1)) / char_hgpref_0_xT_eig_HG[inx_L];
	z charf = char_exppref_u[inx_N]
				* exp(-0.5*tracemul(V_0u_inv[inx_N],(zmatrix)xT))
				* (h1);

	return charf;

}

//Evaluate characteristic function for derivative:
Sim_WMSV::z
Sim_WMSV::eval_char(const matrix &xT) {

	zmatrix tmp1 = char_hgpref_h * xT;
	z h1 = mhg(M, 1.0, mhg_p, mhg_q, eig2(tmp1));
	
	tmp1 = char_hgpref_0 * xT;
	z h2 = mhg(M, 1.0, mhg_p, mhg_q, eig2(tmp1));
	
	tmp1 = (V_0h_inv - V_00_inv);
	zmatrix tmp2 = tmp1*xT;
	z charf = char_exppref_h
				* exp(-0.5*trace(tmp2))
				* (h1/h2);
				
	return charf;

}

//Get simulated values:
Sim_WMSV::fvec
Sim_WMSV::get_Y() {

	return Y;
	
}

//Get option price based on simulated values:
double
Sim_WMSV::get_OP_call(const double K) {

	double price = 0;
	
	for (int i=1; i<=Y.length(); ++i) {
	
		//European call option:
		price += std::max(0.0, exp(Y(i)) - K);
		
	}

	return (price / (double)Y.length())*exp(-r*t);
	
}

//Get standard error:
double
Sim_WMSV::get_sim_sderr() {

	double mean = 0;
	
	for (int i=1; i<=Y.length(); ++i) {
		
		mean += exp(Y(i));
		
	}
	
	mean /= (double)Y.length();
	
	double sd = 0;
	for (int i=1; i<=Y.length(); ++i) {
	
		sd += pow((exp(Y(i))-mean),2.0);
		
	}
	
	sd /= ((double)Y.length()-1.0);
	sd = sqrt(sd);
	
	return sd / sqrt((double)Y.length());

}

//Get confidence interval:
double 
Sim_WMSV::get_CL(const double true_val, const double sim_val, const double sderr) {

	double alpha = std::abs(true_val - sim_val) / sderr;
	return 2*(1-stdNormalCDF(alpha));
	
}

//Write simulations to file:
void
Sim_WMSV::write_sims(std::string filename) {

	std::fstream f;
	f.open(filename.c_str(), std::ios::out);
	
	if (f.is_open()) {
		
		for (int i=1; i<=Y.length(); ++i) {
		
			f << std::setprecision(15) << Y(i) << std::endl;
			
		}
		
		f.close();
		
	}
	
}
		
		
//-----------Posix wrapper-------------//

//Constructors:
Sim_WMSV::PSX_wrapper::PSX_wrapper(Sim_WMSV &_parent, int _step, int _start, int _count)
		:
			parent(_parent),
			step(_step),
			start(_start),
			count(_count)
{
}

Sim_WMSV::PSX_wrapper::PSX_wrapper(Sim_WMSV &_parent, int _step, int _start, int _count, Sim_WMSV::PSX_wrapper::fvec *_params)
		:
			parent(_parent),
			step(_step),
			start(_start),
			count(_count),
			params(_params)
{
}

Sim_WMSV::PSX_wrapper::PSX_wrapper(Sim_WMSV &_parent, int _step, int _startN, int _countN,
												int _startL, int _countL, double _h_in)
		:
			parent(_parent),
			step(_step),
			startN(_startN),
			startL(_startL),
			countN(_countN),
			countL(_countL),
			h_in(_h_in)
{
}

//What to do when initialised as thread:
void
Sim_WMSV::PSX_wrapper::operator()() {
	
	//Execute the appropriate Step function from the parent class:
	if (step == 0) {
		parent.step0MP(start,count);
	}
	else if (step == 1) {
		(*params) = parent.step1(start,count);
	}
	else if (step == 2) {
		parent.step2MP(start,count);
	}
	else if (step == 3) {
		parent.step3(start,count);
	}
	else if (step == 5) {
		parent.step5(start,count);
	}

}

#endif	//SIM_WMSV_CPP
