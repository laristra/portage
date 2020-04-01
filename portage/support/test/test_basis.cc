/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <vector>
#include <array>
#include <random>
#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/support/basis.h"
#include "wonton/support/Point.h"

using Portage::Meshfree::basis::Type;
using Portage::Meshfree::basis::Unitary;
using Portage::Meshfree::basis::Linear;
using Portage::Meshfree::basis::Quadratic;
using Portage::Meshfree::basis::Traits;
using Portage::Meshfree::basis::function;
using Portage::Meshfree::basis::function_size;
using Portage::Meshfree::basis::jet_size;
using Portage::Meshfree::basis::shift;
using Portage::Meshfree::basis::jet;
using Portage::Meshfree::basis::inverse_jet;
using Portage::Meshfree::basis::transfactor;
using Wonton::Point;

// Templated test class that takes type of basis (Unitary, Linear,
// Quadratic) and the space dimension

class BasisTest : public ::testing::Test {

  template<size_t n>
  using Vector = std::array<double, n>;

  template<size_t n, size_t m>
  using Matrix = std::array<Vector<m>, n>;

public:
  template<size_t n, size_t m, size_t z>
  Matrix<n,z> multiply(Matrix<n,m> A, Matrix<m,z> B) {
    Matrix<n,z> R{0.0};
    for (size_t i = 0; i < n; i++)
      for (size_t k = 0; k < z; k++)
        for (size_t j = 0; j < m; j++)
          R[i][k] += A[i][j] * B[j][k];
    return R;
  }

  template<size_t n, size_t m>
  Vector<n> multiply(Matrix<n,m> A, Vector<m> u) {
    Vector<n> v{0.0};
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < m; j++)
        v[i] += A[i][j] * u[j];
    return v;
  }

  template<size_t n, size_t m>
  bool is_identity(Matrix<n,m> mat) {
    double const eps = 1.0e-12;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < m; j++)
        if ((i == j and std::fabs(mat[i][j] - 1) > eps) or
            (i != j and std::fabs(mat[i][j]) > eps))
          return false;
    return true;
  }

  template<size_t n, size_t m>
  bool is_equal(Matrix<n,m> const& A, Matrix<n,m> const& B, double eps) {
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < m; j++)
        if (std::fabs(A[i][j] - B[i][j]) > eps)
          return false;
    return true;
  }

  template<size_t n>
  bool is_equal(Vector<n> const& u, Vector<n> const& v, double eps) {
    for (size_t i = 0; i < n; i++)
      if (std::fabs(u[i] - v[i]) > eps)
        return false;
    return true;
  }

  template<Type type, int dim>
  void checkBasis() {

    std::random_device device;
    std::mt19937 engine { device() };
    std::uniform_real_distribution<double> generator(0.0, 1.0);

    Point<dim> x, y, c, xc;

    for (int i = 0; i < dim; ++i) {
      x[i] = generator(engine);
      y[i] = generator(engine);
      c[i] = 8.0 * generator(engine);
      xc[i] = x[i] + c[i];
    }

    auto bf_x          = function<type, dim>(x);
    auto bf_y          = function<type, dim>(y);
    auto bf_y_minus_x  = function<type, dim>(Point<dim>(y - x));
    auto bf_shifted_xy = shift<type, dim>(x, y);
    auto bj_x          = jet<type, dim>(x);
    auto bj_negx       = jet<type, dim>(-1.0 * x);
    auto bjinv_x       = inverse_jet<type, dim>(x);
    auto bf_xc         = function<type, dim>(xc);
    int const N        = Traits<type, dim>::function_size;

    // Check that J(x)*Jinverse(x) is identity
    ASSERT_TRUE(is_identity(multiply(bj_x, bjinv_x)));
    ASSERT_TRUE(is_identity(multiply(bjinv_x, bj_x)));

    // Check that Jinverse(x) = J(-x)
    ASSERT_TRUE(is_equal(bj_negx, bjinv_x, 1.E-12));

    // Check that b(x) = J(x).e0 or bf_x1 = bj_x1*e0
    Vector<N> e0 {1.0};
    ASSERT_TRUE(is_equal(bf_x, multiply(bj_x, e0), 1.E-12));

    // Check that b(y-x) = b_shifted(x,y)
    ASSERT_TRUE(is_equal(bf_shifted_xy, bf_y_minus_x, 1.0e-12));

    // Check that b_shifted(x,y) = Jinverse(x)*b(y)
    ASSERT_TRUE(is_equal(bf_shifted_xy, multiply(bjinv_x, bf_y), 1.0e-8));

    // check external-facing size functions are correct
    auto fs0 = function_size<dim>(type);
    ASSERT_EQ(fs0, N);

    auto ks = Traits<type, dim>::jet_size;
    auto js0 = jet_size<dim>(type)[0];
    auto js1 = jet_size<dim>(type)[1];
    ASSERT_EQ(js0, ks[0]);
    ASSERT_EQ(js1, ks[1]);

    // Check that vector-valued function is correct
    {
      std::vector<double> result(function<dim>(type, x));
      for (int i = 0; i < N; i++) {
        ASSERT_EQ(bf_x[i], result[i]);
      }
    }

    // Check that vector-valued shift is correct
    {
      std::vector<double> result(shift<dim>(type, x, y));
      for (int i = 0; i < N; i++) {
        ASSERT_EQ(bf_shifted_xy[i], result[i]);
      }
    }

    // Check that vector<vector>-valued jet is correct
    {
      std::vector<std::vector<double>> result(jet<dim>(type, x));
      auto jsize = jet_size<dim>(type);
      for (int i = 0; i < jsize[0]; i++) {
        for (int j = 0; j < jsize[1]; j++) {
          ASSERT_EQ(bj_x[i][j], result[i][j]);
        }
      }
    }

    // Check that vector<vector>-valued inverse_jet is correct
    {
      auto result(Portage::Meshfree::basis::inverse_jet<dim>(type, x));
      auto jsize = jet_size<dim>(type);
      for (int i = 0; i < jsize[0]; i++)
        for (int j = 0; j < jsize[1]; j++)
          ASSERT_EQ(bjinv_x[i][j], result[i][j]);
    }

    // Check that transfactors work correctly
    {
      auto tf = transfactor<type, dim>(c);
      ASSERT_EQ(tf.size(), fs0);
      for (int i = 0; i < fs0; i++)
        ASSERT_EQ(tf[i].size(), fs0);

      std::vector<double> tbf_x(fs0, 0.);
      for (int i = 0; i < fs0; i++)
        for (int j = 0; j < fs0; j++)
          tbf_x[i] += tf[i][j] * bf_x[j];

      for (int i = 0; i < fs0; i++) {
        ASSERT_NEAR(tbf_x[i], bf_xc[i], 1.e-12);
      }


      auto tf2 = transfactor<dim>(type, c);
      for (int i = 0; i < fs0; i++)
        for (int j = 0; j < fs0; j++)
          ASSERT_EQ(tf[i][j], tf2[i][j]);
    }
  }

};

TEST_F(BasisTest, Unitary_1D) { checkBasis<Unitary, 1>(); }

TEST_F(BasisTest, Unitary_2D) { checkBasis<Unitary, 2>(); }

TEST_F(BasisTest, Unitary_3D) { checkBasis<Unitary, 3>(); }

TEST_F(BasisTest, Linear_1D) { checkBasis<Linear, 1>(); }

TEST_F(BasisTest, Linear_2D) { checkBasis<Linear, 2>(); }

TEST_F(BasisTest, Linear_3D) { checkBasis<Linear, 3>(); }

TEST_F(BasisTest, Quadratic_1D) { checkBasis<Quadratic, 1>(); }

TEST_F(BasisTest, Quadratic_2D) { checkBasis<Quadratic, 2>(); }

TEST_F(BasisTest, Quadratic_3D) { checkBasis<Quadratic, 3>(); }
