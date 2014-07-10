//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Header to include all utility code files in this folder.
             
*/

#ifndef UTILS_H
#define UTILS_H 1

#ifndef EPS
#define EPS 2.220446049250313e-15
#endif


#include "det.cpp"
#include "eig.cpp"
#include "erfs.h"
#include "eye.cpp"
#include "inv.cpp"
#include "mat_solver.cpp"
#include "matmult.cpp"
#include "matsqrt_diag.cpp"
#include "mhg.cpp"
#include "norm1.cpp"
#include "normcdf.cpp"
#include "norminv.cpp"
#include "pos_part.cpp"
#include "timer.h"
#include "trace.cpp"
#include "tracemul.cpp"

//The matrix exponential algorithm requires the above functions:
#include "expm.cpp"

//Ext_chol requires cxxblas:
#include "cxxblas/ext_chol.cpp"


#endif	//UTILS_H