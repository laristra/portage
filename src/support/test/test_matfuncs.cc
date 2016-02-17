/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/matrix.h"

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

  std::vector<std::vector<double>> A = {{1, 2.5, -1, 3.0},
                                        {2, -3, -2, 1.0},
                                        {0, 1.1, 2.3, -3.0}};


  std::vector<double> X = {3.0, 2.0, 1.0, -1.0};

  std::vector<double> AX;
  std::vector<double> AX_expected = {4.0,-3.0,7.5};
  
  Portage::matvec_multiply(A,X,&AX);

  ASSERT_EQ(A.size(),AX.size());

  for (int i = 0; i < AX_expected.size(); ++i)
    ASSERT_EQ(AX_expected[i],AX[i]);

}

/*!
  @brief Test the matrix-matrix multiply functionality
 */
TEST(Matrix, MatMatMult) {

  // Use rectangular matrices to catch problems with dimensions and indices
  std::vector<std::vector<double>> A = {{1, 2.5, -1, 3.0},
                                        {2, -3, -2, 1.0},
                                        {0, 1.1, 2.3, -3.0}};

  std::vector<std::vector<double>> B = {{3.0, 2.0, 1.0},
                                        {2.0,-3.0,-5.0},
                                        {2.5,2.6,0.0},
                                        {5.5,-2.0,3.5}};

  std::vector<std::vector<double>> AB;
  std::vector<std::vector<double>> AB_expected = {{22,-14.1,-1},
                                                  {0.5,5.8,20.5},
                                                  {-8.55,8.68,-16.0}};

  
  Portage::matmat_multiply(A,B,&AB);

  for (int i = 0; i < AB_expected.size(); ++i)
    for (int j = 0; j < AB_expected[i].size(); ++j)
      ASSERT_EQ(AB_expected[i][j],AB[i][j]);

}

/*!
  @brief Test the matrix transpose functionality
 */
TEST(Matrix, MatTranspose) {

  // Use rectangular matrix to catch problems with dimensions and indices
  std::vector<std::vector<double>> A = {{1, 2.5, -1, 3.0},
                                        {2, -3, -2, 1.0},
                                        {0, 1.1, 2.3, -3.0}};


  std::vector<std::vector<double>> AT;

  Portage::mat_transpose(A,&AT);
   
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < A[i].size(); ++j)
      ASSERT_EQ(A[i][j],AT[j][i]);

}

/*! 
  @brief Test the matrix inverse functionality
 */
TEST(Matrix, MatInverse) {

  // Use rectangular matrix to catch problems with dimensions and indices
  std::vector<std::vector<double>> A = {{1, 2.5, -1, 3.0},
                                        {2, -3, -2, 1.0},
                                        {0, 1.1, 2.3, -3.0},
                                        {-2, 1, -1, 2}};

  std::vector<std::vector<double>> Ainv;
  std::vector<std::vector<double>> Ainv_expected = 
      {{0.213592,0.097087,0.048544,-0.296117},
       {0.407767,0.912621,1.456311,1.116505},
       {-0.524272,-2.601942,-3.300971,-2.864078},
       {-0.252427,-1.660194,-2.330097,-1.786408}};

  Portage::mat_invert(A,&Ainv);
   
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < A.size(); ++j)
      ASSERT_NEAR(Ainv_expected[i][j],Ainv[i][j],1.0e-6);

}

