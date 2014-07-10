//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Enables chaining matrix multiplications, since FLENs will not allow
			 temporary matrices (unavoidable when chaining).
			 Operator overloading allows defining real and complex versions independently.
             
*/

#ifndef UTILS_MATMULT_CPP
#define UTILS_MATMULT_CPP 1

#include <complex>

#include <flens/flens.cxx>

typedef flens::GeMatrix<flens::FullStorage<double> > matrix;
typedef flens::GeMatrix<flens::FullStorage<std::complex<double> > > zmatrix;


//To multiply 3 matrices:
void
mat3mult(matrix &out, const matrix &M1, const matrix &M2, 
								const matrix &M3) {
			
	matrix tmp;
	
	tmp = M1 * M2;
	out = tmp * M3;

}

//To multiply 3 complex matrices:
void
mat3mult(zmatrix &out, const zmatrix &M1, const zmatrix &M2, 
								const zmatrix &M3) {
			
	zmatrix tmp;
	
	tmp = M1 * M2;
	out = tmp * M3;

}

//To multiply 4 matrices:
void
mat4mult(matrix &out, const matrix &M1, const matrix &M2, 
								const matrix &M3, const matrix &M4) {
	
	matrix tmp;
							
	out = M1 * M2;
	tmp = out * M3;
	
	out = tmp * M4;
	
}

//To multiply 4 complex matrices:
void
mat4mult(zmatrix &out, const zmatrix &M1, const zmatrix &M2, 
								const zmatrix &M3, const zmatrix &M4) {
	
	zmatrix tmp;
							
	out = M1 * M2;
	tmp = out * M3;
	
	out = tmp * M4;
	
}

//To multiply 5 matrices:
void
mat5mult(matrix &out, const matrix &M1, const matrix &M2, 
								const matrix &M3, const matrix &M4, const matrix &M5) {
	
	matrix tmp;					
	tmp = M1 * M2;
	out = tmp * M3;
	
	tmp = out * M4;
	out = tmp * M5;

}

//To multiply 5 complex matrices:
void
mat5mult(zmatrix &out, const zmatrix &M1, const zmatrix &M2, 
								const zmatrix &M3, const zmatrix &M4, const zmatrix &M5) {
		
	zmatrix tmp;					
	tmp = M1 * M2;
	out = tmp * M3;
	
	tmp = out * M4;
	out = tmp * M5;

}

#endif	//UTILS_MATMULT_CPP