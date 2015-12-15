/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef PORTAGE_MATRIX_H
#define PORTAGE_MATRIX_H

#include <cassert>
#include <vector>
#include <iostream>


/*!
  \brief matrix helper functions for portage 
*/

/*! \todo Note: We could convert these to a functional form but that
  entails a copy out of a potentially bulky data structure */

namespace Portage {

/// Matrix-Vector multiply

void matvec_multiply(std::vector<std::vector<double>> A, std::vector<double> X,
                     std::vector<double> *AX) { 
  int m = A.size();
  int n = A[0].size(); // assume that all other rows are of same size
  assert(n == X.size());

  (*AX).resize(m,0.0);  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j)
      (*AX)[i] += A[i][j]*X[j];
  }
}

/// Matrix-Matrix multiply

void matmat_multiply(std::vector<std::vector<double>> A, 
                     std::vector<std::vector<double>> B,
                     std::vector<std::vector<double>> *AB) {
  int ma = A.size();
  int na = A[0].size(); // assume all other rows are of same size
  int mb = B.size();
  int nb = B[0].size(); // assume all other rows are of same size 

  assert(na == mb);

  (*AB).resize(ma);
  for (int i = 0; i < ma; ++i)
    (*AB)[i].resize(nb,0.0);

  for (int i = 0; i < ma; ++i)
    for (int j = 0; j < nb; ++j)
      for (int k = 0; k < na; ++k)
        (*AB)[i][j] += A[i][k]*B[k][j];
}


/// Transpose of a matrix

void mat_transpose(std::vector<std::vector<double>> A,
                   std::vector<std::vector<double>> *AT) {
  int m = A.size();
  int n = A[0].size(); // assume all other rows are of the same size

  (*AT).resize(n);
  for (int i = 0; i < n; ++i)
    (*AT)[i].resize(m);

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      (*AT)[j][i] = A[i][j];
}

/// Inverse of a matrix 

// Code from DSP Design Performance Page by Dr. Jeffrey Tafts
// http://www.nauticom.net/www/jdtaft/FortranMatrix.htm 

void mat_invert(std::vector<std::vector<double>> A,
                std::vector<std::vector<double>> *Ainv) {

  int m = A.size();
  assert(m == A[0].size());

  (*Ainv).resize(m);
  for (int i = 0; i < m; ++i)
    (*Ainv)[i].resize(m);

  // Temporary matrix

  std::vector<std::vector<double>> D(m);
  for (int i = 0; i < m; ++i)
    D[i].resize(2*m);

  // initialize the reduction matrix 
  int m2 = 2*m;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      D[i][j] = A[i][j];
      D[i][m+j] = 0.0;
    }
    D[i][m+i] = 1.0;
  }

  //  do the reduction
  for (int i = 0; i < m; i++) {
    double alpha = D[i][i];
    if (alpha == 0.0) {
      std::cerr << "invert: Singular Matrix" << std::endl;
      return;
    }
    
    for (int j = 0; j < m2; j++)
      D[i][j] = D[i][j]/alpha;

    for (int k = 0; k < m; k++) {
      if ((k-i)== 0) 
	continue;

      double beta = D[k][i];
      for (int j = 0; j < m2; j++)
	D[k][j] = D[k][j] - beta*D[i][j];
    }
  }

  // copy result into output matrix

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      (*Ainv)[i][j] = D[i][j+m];

}

} // namespace Portage

#endif
