/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <vector>
#include <array>

#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/support/basis.h"
#include "wonton/support/Point.h"

using Portage::Meshfree::Basis::Type;
using Portage::Meshfree::Basis::Unitary;
using Portage::Meshfree::Basis::Linear;
using Portage::Meshfree::Basis::Quadratic;
using Portage::Meshfree::Basis::Traits;
using Portage::Meshfree::Basis::function;
using Portage::Meshfree::Basis::function_size;
using Portage::Meshfree::Basis::jet_size;
using Portage::Meshfree::Basis::shift;
using Portage::Meshfree::Basis::jet;
using Portage::Meshfree::Basis::inverse_jet;
using Portage::Meshfree::Basis::transfactor;
using Wonton::Point;
using std::array;

// Templated test class that takes type of basis (Unitary, Linear,
// Quadratic) and the space dimension

class BasisTest : public ::testing::Test {
 public:
  template<size_t N1, size_t N2, size_t N3>
  array<array<double, N1>, N3> matmatmult(array<array<double, N2>, N1> mat1,
                                          array<array<double, N3>, N2> mat2) {
    array<array<double, N3>, N1> matout{0.0};
    for (size_t i = 0; i < N1; i++)
      for (size_t k = 0; k < N3; k++)
        for (size_t j = 0; j < N2; j++)
          matout[i][k] += mat1[i][j]*mat2[j][k];
    return matout;
  }

  template<size_t N1, size_t N2>
  array<double, N1> matvecmult(array<array<double, N2>, N1> mat,
                               array<double, N2> vec) {
    array<double, N1> vecout{0.0};
    for (size_t i = 0; i < N1; i++)
      for (size_t j = 0; j < N2; j++)
        vecout[i] += mat[i][j]*vec[j];
    return vecout;
  }

  template<size_t N1, size_t N2>
  bool is_identity(array<array<double, N2>, N1> mat, double tol) {
    for (size_t i = 0; i < N1; i++)
      for (size_t j = 0; j < N2; j++)
        if ((i == j && fabs(mat[i][j]-1) > tol) ||
            (i != j && fabs(mat[i][j]) > tol))
          return false;
    return true;
  }

  template<size_t N1, size_t N2>
  bool is_equal(array<array<double, N2>, N1> mat1,
                array<array<double, N2>, N1> mat2, double tol) {
    for (size_t i = 0; i < N1; i++)
      for (size_t j = 0; j < N2; j++)
        if (fabs(mat1[i][j]-mat2[i][j]) > tol)
          return false;
    return true;
  }

  template<size_t N1>
  bool is_equal(array<double, N1> vec1,
                array<double, N1> vec2, double tol) {
    for (size_t i = 0; i < N1; i++)
      if (fabs(vec1[i]-vec2[i]) > tol)
        return false;
    return true;
  }

  template<Type type, int Dim>
  void checkBasis() {
    srand(time(NULL));
    Point<Dim> x, y, c, xc;
    for (size_t d = 0; d < Dim; d++)
      x[d] = ((double)rand())/RAND_MAX;
    for (size_t d = 0; d < Dim; d++)
      y[d] = ((double)rand())/RAND_MAX;
    for (size_t d = 0; d < Dim; d++)
      c[d] = 8*((double)rand())/RAND_MAX;
    for (size_t d = 0; d < Dim; d++)
      xc[d] = x[d]+c[d];

    auto bf_x = function<type, Dim>(x);
    auto bf_y = function<type, Dim>(y);
    auto bf_y_minus_x = function<type, Dim>(Point<Dim>(y-x));
    auto bf_shifted_xy = shift<type, Dim>(x, y);
    auto bj_x = jet<type, Dim>(x);
    auto bj_negx = jet<type, Dim>(-1.0*x);
    auto bjinv_x = inverse_jet<type, Dim>(x);
    auto bf_xc = function<type, Dim>(xc);

    // Check that J(x)*Jinverse(x) is identity
    auto j_jinv = matmatmult(bj_x, bjinv_x);
    ASSERT_TRUE(is_identity(j_jinv, 1.0e-12));
    j_jinv = matmatmult(bjinv_x, bj_x);
    ASSERT_TRUE(is_identity(j_jinv, 1.0e-12));

    // Check that Jinverse(x) = J(-x)
    ASSERT_TRUE(is_equal(bj_negx, bjinv_x, 1.0e-12));

    // Check that b(x) = J(x).e0 or bf_x1 = bj_x1*e0
    array<double, Traits<type, Dim>::function_size>
        e0{1.0};
    array<double, Traits<type, Dim>::function_size>
        vec1 = matvecmult(bj_x, e0);
    ASSERT_TRUE(is_equal(bf_x, vec1, 1.0e-12));

    // Check that b(y-x) = b_shifted(x,y)
    ASSERT_TRUE(is_equal(bf_shifted_xy, bf_y_minus_x, 1.0e-12));

    // Check that b_shifted(x,y) = Jinverse(x)*b(y)
    array<double, Traits<type, Dim>::function_size>
        vec2 = matvecmult(bjinv_x, bf_y);
    ASSERT_TRUE(is_equal(bf_shifted_xy, vec2, 1.0e-8));

    // check external-facing size functions are correct
    size_t fs0, gs0, js0, js1, ks0, ks1;
    fs0=function_size<Dim>(type);
    gs0=Traits<type,Dim>::function_size;
    ASSERT_EQ(fs0, gs0);
    array<size_t,2> ks = Traits<type,Dim>::jet_size;
    js0=jet_size<Dim>(type)[0];
    ASSERT_EQ(js0, ks[0]);
    js1=jet_size<Dim>(type)[1];
    ASSERT_EQ(js1, ks[1]);

    // Check that vector-valued function is correct
    {
      std::vector<double> result(function<Dim>(type, x));
      for (int i=0; i<Traits<type, Dim>::function_size; i++) {
        ASSERT_EQ(bf_x[i], result[i]);
      }
    }

    // Check that vector-valued shift is correct
    {
      std::vector<double> result(shift<Dim>(type, x,y));
      for (int i=0; i<Traits<type, Dim>::function_size; i++) {
        ASSERT_EQ(bf_shifted_xy[i], result[i]);
      }
    }

    // Check that vector<vector>-valued jet is correct
    {
      std::vector<std::vector<double>> result(jet<Dim>(type, x));
      auto jsize = jet_size<Dim>(type);
      for (int i=0; i<jsize[0]; i++) {
        for (int j=0; j<jsize[1]; j++) {
          ASSERT_EQ(bj_x[i][j], result[i][j]);
        }
      }
    }

    // Check that vector<vector>-valued inverse_jet is correct
    {
      auto result(Portage::Meshfree::Basis::inverse_jet<Dim>(type, x));
      auto jsize = jet_size<Dim>(type);
      for (int i=0; i<jsize[0]; i++)
	for (int j=0; j<jsize[1]; j++) {
          ASSERT_EQ(bjinv_x[i][j], result[i][j]);
	}
    }

    // Check that transfactors work correctly
    {
      typename Traits<type,Dim>::matrix_t tf = transfactor<type,Dim>(c);
      ASSERT_EQ(tf.size(), fs0);
      for (int i=0; i<fs0; i++) ASSERT_EQ(tf[i].size(), fs0);
      std::vector<double> tbf_x(fs0,0.);
      for (int i=0; i<fs0; i++) for (int j=0; j<fs0; j++) tbf_x[i] += tf[i][j]*bf_x[j];
      for (int i=0; i<fs0; i++) ASSERT_NEAR(tbf_x[i], bf_xc[i], 1.e-12);

      auto tf2 = transfactor<Dim>(type,c);
      for (int i=0; i<fs0; i++) for (int j=0; j<fs0; j++) ASSERT_EQ(tf[i][j], tf2[i][j]);
    }
  }

};

TEST_F(BasisTest, Unitary_1D) {
  checkBasis<Unitary, 1>();
}

TEST_F(BasisTest, Unitary_2D) {
  checkBasis<Unitary, 2>();
}

TEST_F(BasisTest, Unitary_3D) {
  checkBasis<Unitary, 3>();
}

TEST_F(BasisTest, Linear_1D) {
  checkBasis<Linear, 1>();
}

TEST_F(BasisTest, Linear_2D) {
  checkBasis<Linear, 2>();
}

TEST_F(BasisTest, Linear_3D) {
  checkBasis<Linear, 3>();
}

TEST_F(BasisTest, Quadratic_1D) {
  checkBasis<Quadratic, 1>();
}

TEST_F(BasisTest, Quadratic_2D) {
  checkBasis<Quadratic, 2>();
}

TEST_F(BasisTest, Quadratic_3D) {
  checkBasis<Quadratic, 3>();
}
