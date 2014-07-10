#include <iostream>


#include <flens/flens.cxx>

#include "ext_chol.h"

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef Matrix::IndexType                   IndexType;

    const IndexType n = 3;
	
	Matrix orig(n,n);
    orig =  
    23,3,7,
    3,11,9,
    7,9,13;
  	
	Matrix A = orig.lower();
	Matrix cr(n,n);
	Matrix kr(n,n);
	Matrix p(n,n);
	
    cerr << "orig = " << A << endl;

    int Rank = ext_chol(cr,kr,p,A);

	cout << "Rank = " << Rank << endl;
	cout << "cr = " << cr << endl;
	cout << "kr = " << kr << endl;
	cout << "p = " << p << endl;
    
    
    
    
    
}