/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

/* 
 * svd.h - Double precision SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine was adapted from Numerical Recipes SVD routines in C
 * to be called with C++ std::vectors by Michael Rogers.
 *
 * Input to svd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

#ifndef SRC_SVD_H_
#define SRC_SVD_H_

#include <vector>

namespace Portage {

// Perform singular value decomposition

int svd(std::vector<std::vector<double> > & a, std::vector<double> & w, std::vector<std::vector<double> > & v);

// Singular value backsubstitution.

void svd_solve(const std::vector<std::vector<double> > & u, const std::vector<double> & w, \
	       const std::vector<std::vector<double> > & v, \
	       std::vector<double> & b, std::vector<double> & x);
 // Solves A*x=b for vector x.  A is specified by the arrays
 //   u[1...m][1...n], w[1...n],v[1...n[1...n], as returned by
 //   svd.  m and n are the dimensions of a. b[1...m] is the rhs.
 //   x[1...n] is the solution vector (output). Nothing is destroyed
 //   so it can be called sequentially with different b vectors. 

}  // namespace Portage

#endif
