#ifndef EXT_CHOL_CPP
#define EXT_CHOL_CPP 1

using namespace flens;


#ifndef PACIOLI

#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif
#include <flens/flens.cxx>

//Extended cholesky decomposition via LAPACK dpstrf:
double
ext_chol(GeMatrix<FullStorage<double> > &cr, GeMatrix<FullStorage<double> > &kr, 
GeMatrix<FullStorage<double> > &p, const GeMatrix<FullStorage<double> > &A) {

	//cr, kr p, A should all have same dimension, d x d.
	//Assumes p is full of zeros.
	//Returns rank. cr, kr, p all returned with d x d dimension. 
	// Only top left 1:rank,1:rank is relevant for cr
	// Only top left 1:(d-rank),1:(d-rank) is relevant for kr
	
	typedef typename GeMatrix<FullStorage<double> >::IndexType          IndexType;
	
	typedef DenseVector<Array<double> >         Vector;
    typedef DenseVector<Array<IndexType> >      IndexVector;
    
    //Copy since input to dpstrf is mutilated (we only need lower part)
	GeMatrix<FullStorage<double> > C = A.lower();
	
	const IndexType N = A.numRows();
	
	IndexType 		Rank;
	IndexVector		piv(N);
	Vector			work(2*N);
	
	cxxlapack::pstrf('L', 
    				  N, 
    				  C.data(), 
    				  C.leadingDimension(), 
    				  piv.data(), 
    				  Rank, 
    				  -1, 
    				  work.data());
   
    //Extract cr:			  
	for (int i=1; i<=Rank; ++i) {
		for (int j=1; j<=i; ++j) {
		
			cr(i,j) = C(i,j);
			
		}
	}
	
	//cr.resize(Rank,Rank);
	
	//Extract kr:
	if (Rank < N) {
	
		for (int i=Rank+1; i<=N; ++i) {
			for (int j=1; j<=Rank; ++j) {
		
				kr(i-Rank,j) = C(i,j);
			
			}
		}
		
		//kr.resize(N-Rank,Rank);
		
	}
	//else {
	
		//kr.resize(0,0);
		
	//}
	
	//Build p:
	for (int i=1; i<=N; ++i) {
	
		p(piv(i),i) = 1;
		
	}
	
	return Rank;
	
}

#endif	//PACIOLI

#ifdef PACIOLI

#include <flens/flens.cxx>

#include <lapacke.h>

//Extended cholesky decomposition via LAPACK dpstrf:
double
ext_chol(GeMatrix<FullStorage<double> > &cr, GeMatrix<FullStorage<double> > &kr, 
GeMatrix<FullStorage<double> > &p, const GeMatrix<FullStorage<double> > &A) {

	//cr, kr p, A should all have same dimension, d x d.
	//Assumes p is full of zeros.
	//Returns rank. cr, kr, p all returned with d x d dimension. 
	// Only top left 1:rank,1:rank is relevant for cr
	// Only top left 1:(d-rank),1:(d-rank) is relevant for kr
	
	typedef lapack_int		            		IndexType;
	
	typedef DenseVector<Array<double> >         Vector;
    typedef DenseVector<Array<IndexType> >      IndexVector;
    
    //Copy since input to dpstrf is mutilated (we only need lower part)
	GeMatrix<FullStorage<double> > C = A.lower();
	
	const IndexType N = A.numRows();
	
	IndexType 		Rank;
	IndexVector		piv(N);
	
	LAPACKE_dpstrf(LAPACK_COL_MAJOR,
				   'L', 
    				N, 
    				C.data(), 
    				C.leadingDimension(), 
    				piv.data(), 
    				&Rank, 
    				-1.0);
    				
   
    //Extract cr:			  
	for (int i=1; i<=Rank; ++i) {
		for (int j=1; j<=i; ++j) {
		
			cr(i,j) = C(i,j);
			
		}
	}
	
	//cr.resize(Rank,Rank);
	
	//Extract kr:
	if (Rank < N) {
	
		for (int i=Rank+1; i<=N; ++i) {
			for (int j=1; j<=Rank; ++j) {
		
				kr(i-Rank,j) = C(i,j);
			
			}
		}
		
		//kr.resize(N-Rank,Rank);
		
	}
	//else {
	
		//kr.resize(0,0);
		
	//}
	
	//Build p:
	for (int i=1; i<=N; ++i) {
	
		p(piv(i),i) = 1;
		
	}
	
	return Rank;
	
}

#endif	//PACIOLI

#endif	//UTILS_TRACE_CPP