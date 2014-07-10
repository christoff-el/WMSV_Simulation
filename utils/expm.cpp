//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Compute the matrix exponential using the algorithm in
			 Appendix 1.
             
*/

#ifndef EXPM_CPP
#define EXPM_CPP 1

#include <flens/flens.cxx>

using namespace flens;


//Constants:
double B0 =  64764752532480000;
double B1 =  32382376266240000;
double B2 =  7771770303897600;
double B3 =  1187353796428800;
double B4 =  129060195264000;
double B5 =  10559470521600;
double B6 =  670442572800;
double B7 =  33522128640;
double B8 =  1323241920;
double B9 =  40840800;
double B10 = 960960;
double B11 = 16380;
double B12 = 182;
double B13 = 1;

double th3 =  1.495585217958292e-2;
double th5 =  2.539398330063230e-1;
double th7 =  9.504178996162932e-1;
double th9 =  2.097847961257068e0;
double th13 = 5.371920351148152e0;

typedef GeMatrix<FullStorage<double> >	matrix;


matrix
expm(matrix A) {

	int n = A.numRows();
	
	double mu = (trace(A)/(double)n);
	A = A - mu*eye(n);
	
	double A1norm = norm1(A);
	
	matrix U(n,n), Utmp(n,n);
	matrix V(n,n), Vtmp(n,n);
	
	
	//m=13
	if (A1norm > th9) {
	
		int s = log2(A1norm/th13);
		s++;
	
		A = A * (1/pow(2.0,(double)s));
	
		matrix A2=A*A;
		matrix A4=A2*A2;
		matrix A6=A2*A4;
		matrix A8=A4*A4;
	
    	Utmp 	= B13*A4 + B11*A2+ eye(n,B9);
	    U 		= A8 * Utmp;
    	Utmp 	= U + B7*A6+ B5*A4+ B3*A2 + eye(n,B1);
	    U	 	= A*Utmp;
		//U = A*(A6*(B13*A6+B11*A4+B9*A2) + B7*A6 + B5*A4 + B3*A2+B1*eye(n));
	
		V	 	= B12*A4 + B10*A2 + eye(n,B8);
		Vtmp 	= A8*V;
		V	 	= Vtmp + B6*A6 + B4*A4+ B2*A2 + eye(n,B0);
		//V = A6*(B12*A6 + B10*A4 + B8*A2) + B6*A6 + B4*A4 + B2*A2 + B0*eye(n);
		
		mat_solver(A,(V-U),(U+V));
	
		for (int i=1; i<=s; ++i) {
	
			U = A;
			V = A;
			A = U*V;
		
		}
	
		return exp(mu)*A;
	
	}
	//Else:
	
	matrix A2 = A*A;
	
	//m = 3
	if (A1norm < th3) {
	
		Utmp 	= eye(n,B1) + B3*A2;
		U 		= Utmp * A;
		V = eye(n,B0) + B2*A2;
		
		mat_solver(A,(V-U),(U+V));
		
		return exp(mu)*A;
		
	}
	
	matrix A4 = A2*A2;
	
	//m=5
	if (A1norm < th5) {
	
		Utmp	= eye(n,B1) + B3*A2 + B5*A4;
		U		= Utmp * A;
		V = eye(n,B0) + B2*A2 + B4*A4;
		
		mat_solver(A,(V-U),(U+V));
		
		return exp(mu)*A;
		
	}
	
	matrix A6 = A2*A4;
	
	//m=7
	if (A1norm < th7) {
	
		Utmp	= eye(n,B1) + B3*A2 + B5*A4 + B7*A6;
		U		= Utmp * A;
		V = eye(n,B0) + B2*A2 + B4*A4 + B6*A6;
		
		mat_solver(A,(V-U),(U+V));
		
		return exp(mu)*A;
		
	}
	
	matrix A8 = A4*A4;
	//m=9
	//if (A1norm < th9) {
	
		Utmp	= eye(n,B1) + B3*A2 + B5*A4 + B7*A6 + B9*A8;
		U		= Utmp * A;
		V = eye(n,B0) + B2*A2 + B4*A4 + B6*A6 + B8*A8;
		
		mat_solver(A,(V-U),(U+V));
		
		return exp(mu)*A;
		
	//} <-- remove if to suppress warning.
	
}


#endif	//EXPM_CPP