/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/Matrix.h"
#include "portage/support/Vector.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

/*!
  @file test_matfuncs.cc
  @brief Tests for matrix operations in matrix.h
 */

/*!
  @brief Test the matrix-vector multiply functionality
 */
TEST(Matrix, MatVecMult) {

  // Use rectangular matrix to catch problems with dimensions and indices

  Portage::Matrix A({{1, 2.5, -1, 3.0},
                     {2, -3, -2, 1.0},
                     {0, 1.1, 2.3, -3.0}});

  std::vector<double> X({3.0, 2.0, 1.0, -1.0});

  std::vector<double> AX_expected({4.0, -3.0, 7.5});

  // Check multiplication with std::vector
  std::vector<double> AX = A*X;

  for (int i = 0; i < 3; ++i)
    ASSERT_EQ(AX_expected[i], AX[i]);

  // This needs to be a square matrix
  Portage::Matrix B({{1, 2.5, -1},
                     {2, -3, -2},
                     {0, 1.1, 2.3}});

  // Check multiplication with Portage::Vector
  Portage::Vector<3> Y, BY;
  for (int i = 0; i < 3; ++i) Y[i] = X[i];
  BY = B*Y;

  std::vector<double> BY_expected({7.0, -2.0, 4.5});

  for (int i = 0; i < 3; ++i)
    ASSERT_EQ(BY_expected[i], BY[i]);
}

/*!
  @brief Test the matrix-matrix multiply functionality
 */
TEST(Matrix, MatMatMult) {

  // Use rectangular matrices to catch problems with dimensions and indices
  Portage::Matrix A({{1, 2.5, -1, 3.0},
                     {2, -3, -2, 1.0},
                     {0, 1.1, 2.3, -3.0}});

  Portage::Matrix B({{3.0, 2.0, 1.0},
                     {2.0, -3.0, -5.0},
                     {2.5, 2.6, 0.0},
                     {5.5, -2.0, 3.5}});

  Portage::Matrix AB_expected({{22, -14.1, -1},
                               {0.5, 5.8, 20.5},
                               {-8.55, 8.68, -16.0}});

  
  Portage::Matrix AB = A*B;

  for (int i = 0; i < AB_expected.rows(); ++i)
    for (int j = 0; j < AB_expected.columns(); ++j)
      ASSERT_EQ(AB_expected[i][j], AB[i][j]);

}

/*!
  @brief Test the matrix transpose functionality
 */
TEST(Matrix, MatTranspose) {

  // Use rectangular matrix to catch problems with dimensions and indices
  Portage::Matrix A({{1, 2.5, -1, 3.0},
                     {2, -3, -2, 1.0},
                     {0, 1.1, 2.3, -3.0}});


  auto AT = A.transpose();

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 4; ++j)
      ASSERT_EQ(A[i][j], AT[j][i]);

}

/*! 
  @brief Test the matrix inverse functionality
 */
TEST(Matrix, MatInverse) {

  Portage::Matrix A({{1, 2.5, -1, 3.0},
                     {2, -3, -2, 1.0},
                     {0, 1.1, 2.3, -3.0},
                     {-2, 1, -1, 2}});

  Portage::Matrix Ainv = A.inverse();

  Portage::Matrix Ainv_expected({{0.213592, 0.097087, 0.048544, -0.296117},
                                 {0.407767, 0.912621, 1.456311, 1.116505},
                                 {-0.524272, -2.601942, -3.300971, -2.864078},
                                 {-0.252427, -1.660194, -2.330097, -1.786408}});

  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      ASSERT_NEAR(Ainv_expected[i][j], Ainv[i][j], 1.0e-6);

}

