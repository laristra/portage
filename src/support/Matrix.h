/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef PORTAGE_MATRIX_H
#define PORTAGE_MATRIX_H

#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <type_traits>

// #include <Vector.h>   // Portage::Vector

/*!
  @file matrix.h
  @brief Matrix class for Portage 

  @tparam Rows     number of rows
  @tparam Columns  number of columns (defaults to number of rows if omitted)
*/

namespace Portage {

template<unsigned int Rows, unsigned int Columns=Rows>
class Matrix {
 public:
  Matrix() {}  // Default constructor - no initialization

  // Initialize to some constant value
  explicit Matrix(double initval) {
    for (unsigned int i = 0; i < Rows; ++i)
      for (unsigned int j = 0; j < Columns; ++j)
        A[i][j] = initval;
  }

  explicit Matrix(std::array<std::array<double, Columns>, Rows> const& inarray)
      : A(inarray) {
    for (unsigned int i = 0; i < Rows; ++i)
      for (unsigned int j = 0; j < Columns; ++j)
        A[i][j] = inarray[i][j];
  }

  // Rely on default copy and assignment constructors

  // Destructor
  ~Matrix() {}

  //! Return a row of values
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of the matrix

  std::array<double, Columns>& operator[](unsigned int i) {
    return A[i];
  }

  //! Return a row of values that cannot be modified
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of a const matrix

  std::array<double, Columns> const& operator[](unsigned int i) const {
    return A[i];
  }

  //! Get a transpose
  Matrix<Columns, Rows> transpose() const {
    Matrix<Columns, Rows> AT;
    for (unsigned int i = 0; i < Rows; ++i)
      for (unsigned int j = 0; j < Columns; ++j)
        AT[j][i] = A[i][j];
    return AT;
  }

  //! Get Inverse - only if its a square matrix

  Matrix<Rows, Rows> inverse() const {
    Matrix<Rows, Rows> Ainv(0.0);
    if (Rows != Columns) {
      std::cerr << "Matrix is rectangular" << std::endl;
      throw std::exception();
    }
    
    // Create a temporary matrix with twice as many columns
    Matrix<Rows, 2*Rows> D;
    
    // Initialize the reduction matrix
    unsigned int Rows2 = 2*Rows;
    for (unsigned int i = 0; i < Rows; i++) {
      for (unsigned int j = 0; j < Rows; j++) {
        D[i][j] = A[i][j];
        D[i][Rows+j] = 0.0;
      }
      D[i][Rows+i] = 1.0;
    }
    
    // Do the reduction
    for (unsigned int i = 0; i < Rows; i++) {
      double alpha = D[i][i];
      if (alpha == 0.0) {
        std::cerr << "invert: Singular Matrix" << std::endl;
        return Ainv;
      }
      
      for (unsigned int j = 0; j < Rows2; j++)
        D[i][j] = D[i][j]/alpha;
      
      for (unsigned int k = 0; k < Rows; k++) {
        if ((k-i) == 0) 
          continue;
        
        double beta = D[k][i];
        for (unsigned int j = 0; j < Rows2; j++)
          D[k][j] = D[k][j] - beta*D[i][j];
      }
    }
    
    // Copy result unsigned into output matrix
    for (unsigned int i = 0; i < Rows; i++)
      for (unsigned int j = 0; j < Rows; j++)
        Ainv[i][j] = D[i][j + Rows];

    return Ainv;
  }

/*!
  @brief  Matrix-Vector multiply
  @param[in] X The vector to post-multiply with
  
 */

Vector<Rows> operator*(Vector<Columns> const& X) {
  Vector<Rows> AX;
  for (unsigned int i = 0; i < Rows; ++i) {
    AX[i] = 0.0;
    for (unsigned int j = 0; j < Columns; ++j)
      AX[i] += A[i][j]*X[j];
  }
  return AX;
}
  
/*!
  @brief  Matrix-Matrix multiply
  @param[in] B   matrix to post-multiply with

  Not able to make B const& since it complains about 
 */

template<unsigned int Rows2 = Columns, unsigned int Columns2> 
Matrix<Rows, Columns2> operator*(Matrix<Rows2, Columns2> const& B) {
  assert(Rows2 == Columns);

  Matrix<Rows, Columns2> AB(0.0);

  for (unsigned int i = 0; i < Rows; ++i)
    for (unsigned int j = 0; j < Columns2; ++j)
      for (unsigned int k = 0; k < Columns; ++k)
        AB[i][j] += A[i][k]*B[k][j];

  return AB;
}


 private:
  std::array<std::array<double, Columns>, Rows> A;
};

}  // namespace Portage

#endif
