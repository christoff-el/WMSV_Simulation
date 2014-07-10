//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef SIM_WIS_CPP
#define SIM_WIS_CPP 1

#include "sim_wis.h"


//Constructor:
Sim_Wis::Sim_Wis(Rng_Machinery &_rng, matrix _x, double _delta, matrix _H, matrix _Sigma, double _t)
		:
			rng(_rng),
			d(_H.numRows()),
			x(_x),
			delta(_delta),
			H(_H),
			Sigma(_Sigma),
			t(_t),
			sim_bes(rng,delta,d),
			y(d,d),
			theta_t(d,d),
			sim_p(d,matrix(d,d)),
			alg1_perm(d,d),
			alg1_x_tilde(d,d),
			alg1_cr(d-1,d-1),
			alg1_kr(d-1,d-1),
			alg1_cr_inv(d-1,d-1),
			alg1_tmp1(d,d),
			alg1_tmp2(d,d),
			sim_pYp(d,d),
			sim_pyp(d,d),
			sim_Y(d,d)
{
	
	Sim_Wis::precomp();

}

//Copy constructor (no need to precomp):
Sim_Wis::Sim_Wis(const Sim_Wis &rhs) 
		:
			rng(rhs.rng),
			d(rhs.d),
			x(rhs.x),
			delta(rhs.delta),
			H(rhs.H),
			Sigma(rhs.Sigma),
			t(rhs.t),
			sim_bes(rhs.sim_bes),
			rank_q(rhs.rank_q),
			sqrt_t(rhs.sqrt_t),
			y(rhs.y),
			theta_t(rhs.theta_t),
			sim_p(rhs.sim_p),
			alg1_perm(rhs.alg1_perm),
			alg1_x_tilde(d,d),
			alg1_cr(d-1,d-1),
			alg1_kr(d-1,d-1),
			alg1_cr_inv(d-1,d-1),
			alg1_tmp1(d,d),
			alg1_tmp2(d,d),
			sim_pYp(d,d),
			sim_pyp(d,d),
			sim_Y(d,d)
{
}

//Null copy constructor (does nothing):
Sim_Wis::Sim_Wis(const Sim_Wis &rhs, bool null) 
		: 
			rng(rhs.rng),
			sim_bes(rhs.sim_bes,0)
{
}

//Simulator:
matrix
Sim_Wis::simulate() {

	sim_pYp = y;
	
	for (int k=1; k<=rank_q; ++k) {

		//y=yPy where Y is sampled according to WIS_d(pyp,delta,0,e_d^1;t)
		mat3mult(sim_pyp,sim_p[k-1],sim_pYp,sim_p[k-1]);

		Sim_Wis::alg1(sim_Y, sim_pyp);
		
		mat3mult(sim_pYp, sim_p[k-1], sim_Y, sim_p[k-1]);
		
	}

	mat3mult(sim_Y,theta_t,sim_pYp,flens::transpose(theta_t));
	return sim_Y;

}

//Sub-algorithm:
void
Sim_Wis::alg1(matrix &Y, const matrix &x) {

	matrix p2(d-1,d-1);
	
	//Perform extended Cholesky decomposition:
	int rank_pyp_trim = ext_chol(alg1_cr, alg1_kr, p2, x(_(2,d),_(2,d)));

	//Set perm = [1 0;0 p]:
	//perm(1,1) = 1;
	alg1_perm(_(2,d),_(2,d)) = p2;

	//Set x_tilde = perm*x*perm.':
	mat3mult(alg1_x_tilde,alg1_perm,x,flens::transpose(alg1_perm));
	
	//Set initial value u:
	alg1_cr_inv = alg1_cr(_(1,rank_pyp_trim),_(1,rank_pyp_trim));
	inv(alg1_cr_inv);
	vec u(rank_pyp_trim+1);
	
	//u(_(2,rank_pyp_trim+1)) = cr_inv * x_tilde(_(2,rank_pyp_trim+1),_(1,1));
	for (int i=2; i<=rank_pyp_trim+1; ++i) {
		for (int j=1; j<=rank_pyp_trim; ++j) {
			u(i) += alg1_cr_inv(i-1,j) * alg1_x_tilde(j+1,1);
		}
	}			
	
	u(1) = alg1_x_tilde(1,1);
	for (int j=2; j<=rank_pyp_trim+1; ++j) {
		u(1) -= u(j)*u(j);
	}

	//Set (U_t):
	vec U(rank_pyp_trim+1);

	//Sample (U_t)_1,1 as a square bessel process
	// d(U_t)_1,1 = (delta - r)*dt + 2*sqrt((U_t)_1,1)*dW_t
	//U(1) = sim_bes.simulate(rank_pyp_trim,u(1),t);

	for (int i=2; i<=rank_pyp_trim+1; ++i) {
		U(i) = u(i) + sqrt_t * rng.gen_std_norm();
	}
		
	alg1_tmp1(1,1) = 1;
	alg1_tmp2(1,1) = sim_bes.simulate(rank_pyp_trim,u(1),t);

	for (int j=2; j<=rank_pyp_trim+1; ++j) {
		alg1_tmp2(1,1) += U(j)*U(j);
	}
	for (int j=2; j<=rank_pyp_trim+1; ++j) {
			alg1_tmp2(1,j) = U(j);
			alg1_tmp2(j,1) = U(j);
	}
		
	//If x is SPSD, we compile with -DSPSD, to perform some extended Cholesky computations:
	#ifdef SPSD
	if (rank_pyp_trim<d-1) {
		
		alg1_tmp1(_(2,rank_pyp_trim+1),_(2,rank_pyp_trim+1)) = alg1_cr(_(1,rank_pyp_trim),_(1,rank_pyp_trim));
		alg1_tmp1(_(rank_pyp_trim+2,d),_(2,rank_pyp_trim+1)) = alg1_kr(_(1,d-rank_pyp_trim-1),_(1,rank_pyp_trim));
		alg1_tmp1(_(rank_pyp_trim+2,d),_(rank_pyp_trim+2,d)) = eye(d-rank_pyp_trim-1);
		alg1_tmp2(_(2,rank_pyp_trim+1),_(2,rank_pyp_trim+1)) = eye(rank_pyp_trim);
			
	}
	//Otherwise, things are simpler:
	else {
		
		alg1_tmp1(_(2,d),_(2,d)) = alg1_cr(_(1,rank_pyp_trim),_(1,rank_pyp_trim));
		alg1_tmp2(_(2,d),_(2,d)) = eye(rank_pyp_trim);
		
	}
	#endif
	
	#ifndef SPSD
	alg1_tmp1(_(2,d),_(2,d)) = alg1_cr(_(1,rank_pyp_trim),_(1,rank_pyp_trim));
	alg1_tmp2(_(2,d),_(2,d)) = eye(rank_pyp_trim);
	#endif
	
	mat5mult(Y, flens::transpose(alg1_perm), alg1_tmp1, alg1_tmp2, flens::transpose(alg1_tmp1), alg1_perm);
	
}

//Reinitialise with new parameters:
void
Sim_Wis::reinitialise(matrix _x, double _delta, matrix _H, matrix _Sigma, double _t) {

	d = _H.numRows();
	x = _x;
	delta = _delta;
	H = _H;
	Sigma = _Sigma;
	t = _t;
	
	y.resize(d,d);
	theta_t.resize(d,d);
	
	Sim_Wis::precomp();
	
}


//Precomputations for speed:
void
Sim_Wis::precomp() {
	
	//Square root of t:
	sqrt_t = sqrt(t);
	
	//Top left of perm matrix:
	alg1_perm(1,1) = 1;
	
	//Calculate qt = \int_0^t expm(s*H)*Sigma.'*Sigma*expm(s*H.') ds using \cite{Loan:1978}:
	matrix qt(d,d);
	matrix tmp(2*d,2*d);
	matrix Qc = flens::transpose(Sigma) * Sigma;

	tmp(_(1,d),_(1,d)) = -1.0*H;
	tmp(_(1,d),_(d+1,2*d)) = Qc;
	tmp(_(d+1,2*d),_(d+1,2*d)) = flens::transpose(H);

	matrix e = expm(t*tmp);

	qt = flens::transpose(e(_(d+1,2*d),_(d+1,2*d))) * (e(_(1,d),_(d+1,2*d)));
	
	//Extended cholesky decomposition of qt/t:
	qt = (1.0/t)*qt;
	
	matrix cr(d,d);
	matrix kr(d,d);
	matrix p(d,d);
	rank_q = ext_chol(cr,kr,p,qt);
	
	//Set theta_t = p^-1[cn 0; kn I_(d-n)]:
	inv(p);			//mangles p, but we don't need it anymore.

	if (rank_q<d) {
	
		matrix tmpd(d,d);
		tmpd(_(1,rank_q),_(1,rank_q)) = cr;
		tmpd(_(rank_q+1,d),_(1,rank_q)) = kr;
		tmpd(_(rank_q+1,d),_(rank_q+1,d)) = eye(d-rank_q);
		
		theta_t = p * tmpd;
		
	}
	else {
	
		theta_t = p * cr;
	
	}
	
	//Precompute necessary sim_p matrices for the given rank_q:
	sim_p.resize(rank_q);
	
	//Set p_k,1 = p_1,k = p_i,i = 1 for i not 1 or k, and 0 otherwise:
	for (int k=1; k<=rank_q; ++k) {
		sim_p[k-1](1,1) = 0;
		sim_p[k-1](k,k) = 0;
		sim_p[k-1](k,1) = 1;
		sim_p[k-1](1,k) = 1;
		
		for (int i=2; i<=k-1; ++i) {
			sim_p[k-1](i,1) = 0;
			sim_p[k-1](1,i) = 0;
		}
		for (int i=k+1; i<=d; ++i) {
			sim_p[k-1](i,1) = 0;
			sim_p[k-1](1,i) = 0;
		}

		for (int i=2; i<=k-1; ++i) {
			sim_p[k-1](i,i) = 1;
		}
		for (int i=k+1; i<=d; ++i) {
			sim_p[k-1](i,i) = 1;
		}
		for (int i=2; i<=d-1; ++i) {
			for (int j=i+1; j<=d; ++j) {
				sim_p[k-1](i,j) = 0;
			}
		}
	}

	//m_t = exp(tH):
	matrix m_t = expm(t*H);

	//y = inv(theta_t) * m_t * x * m_t.' * (inv(theta_t)).';
	matrix theta_inv = theta_t;		//Don't want to mangle theta_t, so take a copy.
	inv(theta_inv);

	mat5mult(y,theta_inv, m_t, x, flens::transpose(m_t), flens::transpose(theta_inv));
	
}

#endif	//SIM_WIS_CPP
