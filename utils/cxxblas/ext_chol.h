#ifndef EXT_CHOL_H
#define EXT_CHOL_H 1

using namespace flens;


//Extended cholesky decomposition via LAPACK dpstrf:
double
ext_chol(GeMatrix<FullStorage<double> > &cr, GeMatrix<FullStorage<double> > &kr, 
GeMatrix<FullStorage<double> > &p, const GeMatrix<FullStorage<double> > &A);


#endif	//EXT_CHOL_H