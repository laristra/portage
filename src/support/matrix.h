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
  @file matrix.h
  @brief Matrix helper functions for Portage 
*/

/*! \todo Note: We could convert these to a functional form but that
  entails a copy out of a potentially bulky data structure */

namespace Portage {

/*!
  @brief  Matrix-Vector multiply
  @param[in] A The input matrix, as a vector of vectors
  @param[in] X The input vector
  @param[in,out] AX Pointer to the resulting vector
 */
void matvec_multiply(const std::vector<std::vector<double>> & A, const std::vector<double> & X,
                     std::vector<double> * const AX) { 
  int m = A.size();
  int n = A[0].size(); // assume that all other rows are of same size
  assert(n == X.size());

  (*AX).resize(m,0.0);  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j)
      (*AX)[i] += A[i][j]*X[j];
  }
}

/*!
  @brief  Matrix-Matrix multiply
  @param[in] A First input matrix, as a vector of vectors
  @param[in] B Second input matrix, as a vector of vectors
  @param[in,out] AB Pointer to the resulting matrix
 */
void matmat_multiply(const std::vector<std::vector<double>> & A, 
                     const std::vector<std::vector<double>> & B,
                     std::vector<std::vector<double>> * const AB) {
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

/*!
  @brief Transpose of a matrix
  @param[in] A Input matrix
  @param[in,out] AT Pointer to the resulting matrix
 */
void mat_transpose(const std::vector<std::vector<double>> & A,
                   std::vector<std::vector<double>> * const AT) {
  int m = A.size();
  int n = A[0].size(); // assume all other rows are of the same size

  (*AT).resize(n);
  for (int i = 0; i < n; ++i)
    (*AT)[i].resize(m);

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      (*AT)[j][i] = A[i][j];
}

/*!
  @brief Inverse of a matrix 
  @param[in] A Input matrix
  @param[in,out] Ainv Pointer to the resulting matrix

  Code from DSP Design Performance Page by Dr. Jeffrey Tafts
  http://www.nauticom.net/www/jdtaft/FortranMatrix.htm 
 */
void mat_invert(const std::vector<std::vector<double>> & A,
                std::vector<std::vector<double>> * const Ainv) {

  int m = A.size();
  assert(m == A[0].size());

  (*Ainv).resize(m);
  for (int i = 0; i < m; ++i)
    (*Ainv)[i].resize(m);

  // Create a temporary matrix
  std::vector<std::vector<double>> D(m);
  for (int i = 0; i < m; ++i)
    D[i].resize(2*m);

  // Initialize the reduction matrix 
  int m2 = 2*m;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      D[i][j] = A[i][j];
      D[i][m+j] = 0.0;
    }
    D[i][m+i] = 1.0;
  }

  // Do the reduction
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

  // Copy result into output matrix
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      (*Ainv)[i][j] = D[i][j+m];

}

} // namespace Portage

#endif
