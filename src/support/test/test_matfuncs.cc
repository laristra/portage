/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include "portage/support/Vector.h"
#include "portage/support/Matrix.h"

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
  @brief Test the matrix-scalar multiply functionality
 */
TEST(Matrix, MatScalarMult) {

  // Use rectangular matrices to catch problems with dimensions and indices
  Portage::Matrix A({{1, 2.5, -1, 3.0},
                     {2, -3, -2, 1.0},
                     {0, 1.1, 2.3, -3.0}});

  Portage::Matrix Ab_expected({{2., 5, -2, 6.0},
                               {4., -6., -4., 2.},
                               {0, 2.2, 4.6, -6.0}});


  Portage::Matrix Ab = A*2.0;

  for (int i = 0; i < Ab_expected.rows(); ++i)
    for (int j = 0; j < Ab_expected.columns(); ++j)
      ASSERT_EQ(Ab_expected[i][j], Ab[i][j]);

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

/*!
  @brief Test the matrix solve functionality
 */
TEST(Matrix, SolveWithInverse) {

  Portage::Matrix A({{1, 2.5, -1, 3.0},
                     {2, -3, -2, 1.0},
                     {0, 1.1, 2.3, -3.0},
                     {-2, 1, -1, 2}});

  Portage::Matrix B({{-0.7296116504854356, -1.313592233009709},
         {13.111650485436888,-1.207766990291261},
      {-27.286407766990273, 1.1242718446601963},
     {-19.37864077669902, 2.9524271844660195}});

  Portage::Matrix AinvB = A.solve(B);

  Portage::Matrix AinvB_expected({{5.530888867942306, -1.2175181449712502},
    {-49.705372796682, 3.2958148741634483},
    {111.84033367895168, -8.335950607974366},
    {76.61442171740964, -5.557187293807148}});

  for (int i = 0; i < B.rows(); ++i)
    for (int j = 0; j < B.columns(); ++j)
      ASSERT_NEAR(AinvB_expected[i][j], AinvB[i][j], 1.0e-6);

}

/*!
  @brief Test the matrix solve functionality
 */
TEST(Matrix, SolveWithJunk) {

  /* Mathematica code to generate this example:
   *
   * a = RotationMatrix[{{1, 2.5, -1., 3.0}, {2., -3, -2, 1.0}}];
   * eval = DiagonalMatrix[{1., 2., 3., 4.}];
   * aa = a.eval.Transpose[a];
   *
   * The matrix aa is symmetric positive-definite with eigenvalues 1,2,3,4.
   */

  Portage::Matrix A({{1.5909092362921544, 0.3023197488806776, -0.0008986359725132367,
    0.5719175069965953}, {0.3023197488806776, 3.684725833270652, 0.22450059345997342,
    -0.7267204571061852}, {-0.0008986359725132367, 0.22450059345997342,
    2.410888035652873, 0.6054854449257954}, {0.5719175069965953, -0.7267204571061852,
    0.6054854449257954, 2.3134768947843307}});

  Portage::Matrix B({{1.2, 3.4}, {-5.6, 1.7}, {9.8, -7.6}, {3.1, 6.2}});

  try {
    Portage::Matrix AinvB = A.solve(B, "la-di-da");
    ASSERT_TRUE(false);
  }
  catch (...) {
    ASSERT_TRUE(true);
  }

}

/*!
  @brief Test the matrix solve functionality
 */
TEST(Matrix, SolveWithPOSV) {

  /* Mathematica code to generate this example:
   *
   * a = RotationMatrix[{{1, 2.5, -1., 3.0}, {2., -3, -2, 1.0}}];
   * eval = DiagonalMatrix[{1., 2., 3., 4.}];
   * aa = a.eval.Transpose[a];
   *
   * The matrix aa is symmetric positive-definite with eigenvalues 1,2,3,4.
   */

  Portage::Matrix A({{1.5909092362921544, 0.3023197488806776, -0.0008986359725132367,
    0.5719175069965953}, {0.3023197488806776, 3.684725833270652, 0.22450059345997342,
    -0.7267204571061852}, {-0.0008986359725132367, 0.22450059345997342,
    2.410888035652873, 0.6054854449257954}, {0.5719175069965953, -0.7267204571061852,
    0.6054854449257954, 2.3134768947843307}});

  Portage::Matrix B({{1.2, 3.4}, {-5.6, 1.7}, {9.8, -7.6}, {3.1, 6.2}});

  Portage::Matrix AinvB = A.solve(B, "lapack-posv");

  Portage::Matrix AinvB_expected({{1.4543699572914695, 0.3263211034432376},
    {-2.077583421788211, 1.5322821892984466},
      {4.470408685674243, -4.354857599247892},
     {-0.8421823492688412, 4.2203641515254}});

  Portage::Matrix resB = A*AinvB;

  for (int i = 0; i < B.rows(); ++i)
    for (int j = 0; j < B.columns(); ++j) {
      ASSERT_NEAR(AinvB_expected[i][j], AinvB[i][j], 1.0e-12);
      ASSERT_NEAR(resB[i][j], B[i][j], 1.e-12);
    }

}

/*!
  @brief Test the matrix solve functionality
 */
TEST(Matrix, SolveWithSYSV) {

  /* Mathematica code:
   *
   */

  Portage::Matrix A({{1.5909092362921544, 0.3023197488806776, -0.0008986359725132367,
    0.5719175069965953}, {0.3023197488806776, 3.684725833270652, 0.22450059345997342,
    -0.7267204571061852}, {-0.0008986359725132367, 0.22450059345997342,
    2.410888035652873, 0.6054854449257954}, {0.5719175069965953, -0.7267204571061852,
    0.6054854449257954, 2.3134768947843307}});

  Portage::Matrix B({{1.2, 3.4}, {-5.6, 1.7}, {9.8, -7.6}, {3.1, 6.2}});

  Portage::Matrix AinvB = A.solve(B, "lapack-sysv");

  Portage::Matrix AinvB_expected({{1.4543699572914695, 0.3263211034432376},
    {-2.077583421788211, 1.5322821892984466},
      {4.470408685674243, -4.354857599247892},
     {-0.8421823492688412, 4.2203641515254}});

  for (int i = 0; i < B.rows(); ++i)
    for (int j = 0; j < B.columns(); ++j)
      ASSERT_NEAR(AinvB_expected[i][j], AinvB[i][j], 1.0e-12);

}

/*!
  @brief Test the matrix solve functionality
 */
TEST(Matrix, SolveWithGESV) {

  /* Mathematica code:
   *
   */

  Portage::Matrix A({{1.5909092362921544, 0.3023197488806776, -0.0008986359725132367,
    0.5719175069965953}, {0.3023197488806776, 3.684725833270652, 0.22450059345997342,
    -0.7267204571061852}, {-0.0008986359725132367, 0.22450059345997342,
    2.410888035652873, 0.6054854449257954}, {0.5719175069965953, -0.7267204571061852,
    0.6054854449257954, 2.3134768947843307}});

  Portage::Matrix B({{1.2, 3.4}, {-5.6, 1.7}, {9.8, -7.6}, {3.1, 6.2}});

  Portage::Matrix AinvB = A.solve(B, "lapack-gesv");

  Portage::Matrix AinvB_expected({{1.4543699572914695, 0.3263211034432376},
    {-2.077583421788211, 1.5322821892984466},
      {4.470408685674243, -4.354857599247892},
     {-0.8421823492688412, 4.2203641515254}});

  for (int i = 0; i < B.rows(); ++i)
    for (int j = 0; j < B.columns(); ++j)
      ASSERT_NEAR(AinvB_expected[i][j], AinvB[i][j], 1.0e-12);

}

/*!
  @brief Test the matrix solve functionality
 */
TEST(Matrix, SolveWithErrorMsg) {

  /* Mathematica code to generate this example:
   *
   * a = RotationMatrix[{{1, 2.5, -1., 3.0}, {2., -3, -2, 1.0}}];
   * eval = DiagonalMatrix[{1., 2., 3., 4.}];
   * aa = a.eval.Transpose[a];
   *
   * The matrix aa is symmetric positive-definite with eigenvalues 1,2,3,4.
   */

  Portage::Matrix A(3,3,0.);
  A[0][0] = 1.; A[1][1] = 1.; A[2][2] = -1.;

  Portage::Matrix B(3,1,0);

  std::string errormsg="ignore";
  Portage::Matrix AinvB = A.solve(B, "lapack-posv", errormsg);
  std::cout << "error message: " << errormsg << std::endl;

  ASSERT_TRUE(errormsg=="ignore");

  errormsg="blahblah";
  AinvB = A.solve(B, "lapack-posv", errormsg);
  std::cout << "error message: " << errormsg << std::endl;

  ASSERT_TRUE(errormsg!="blahblah");

}
