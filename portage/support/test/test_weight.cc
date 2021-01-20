/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <cmath>
#include <ctime>

#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "portage/support/weight.h"
#include "portage/support/portage.h"

using Portage::Meshfree::Weight::Geometry;
using Portage::Meshfree::Weight::ELLIPTIC;
using Portage::Meshfree::Weight::TENSOR;
using Portage::Meshfree::Weight::FACETED;

using Portage::Meshfree::Weight::Kernel;
using Portage::Meshfree::Weight::B4;
using Portage::Meshfree::Weight::SQUARE;
using Portage::Meshfree::Weight::EPANECHNIKOV;
using Portage::Meshfree::Weight::POLYRAMP;
using Portage::Meshfree::Weight::INVSQRT;
using Portage::Meshfree::Weight::COULOMB;
using Portage::Meshfree::Weight::STEP;

using Portage::Meshfree::Weight::eval;

using Portage::Meshfree::Weight::FacetData;
using Portage::Meshfree::Weight::faceted;


using Wonton::Point;
using Wonton::Vector;

using std::tuple;
using std::array;
using std::vector;

using ::testing::TestWithParam;
using ::testing::Combine;
using ::testing::Values;

class WeightTest : public TestWithParam<tuple<Geometry, Kernel>> {
 public:
  template<int Dim>
  void checkWeight(const Geometry geo, const Kernel kernel) {

    //    srand(time(NULL));
    Point<Dim> x;
    for (int d = 0; d < Dim; d++)
      x[d] = ((double) rand()) / RAND_MAX;

    static double nominal_h = .2;
    std::array<double, Dim> h {};
    for (int d = 0; d < Dim; d++)
      h[d] = nominal_h;

    size_t nsides = 2 * Dim;
    vector<FacetData<Dim>> facets(nsides);
    if (geo == FACETED) {
      // set up an isothetic box like TENSOR
      size_t side = 0;
      for (size_t j = 0; j < Dim; j++) {
        for (size_t k = 0; k < Dim; k++)
          facets[side].normal[k] = 0.;
        facets[side].normal[j] = -1.;
        facets[side].smoothing = nominal_h;
        side++;
        for (size_t k = 0; k < Dim; k++)
          facets[side].normal[k] = 0.;
        facets[side].normal[j] = 1.;
        facets[side].smoothing = nominal_h;
        side++;
      }
    }

    double weight_at_x;
    if (geo == FACETED) {
      weight_at_x = faceted<Dim>(kernel, x, x, facets, nsides);
    } else {
      weight_at_x = eval<Dim>(geo, kernel, x, x, h);
    }

    np = powl(10, Dim);
    vector<double> weight(np);

    for (size_t i = 0; i < np; i++) {
      // generate a 'y' val in the range (-3h[d],3h[d])  ----
      // rand/RAND_MAX will give us [0,1]; shift it left by
      // 0.5 and multiply by 6 to give us the desired interval

      Point<Dim> y = x;
      for (size_t d = 0; d < Dim; d++)
        y[d] += 6 * h[d] * (((double) rand()) / RAND_MAX - 0.5);

      if (geo == FACETED) {
        // use explicit faceted weight function
        weight[i] = faceted<Dim>(kernel, x, y, facets, nsides);

        // check equivalence to generic interface
        vector<vector<double>> hvec(nsides,vector<double>(Dim+1));
        for (size_t j=0; j<nsides; j++) {
          for (size_t k=0; k<Dim; k++) hvec[j][k] = facets[j].normal[k];
          hvec[j][Dim] = facets[j].smoothing;
        }
        double alt_wt = eval<Dim>(geo, kernel, x, y, hvec);
        EXPECT_NEAR(alt_wt, weight[i], 1.0e-12);
      } else {
        // use explicit weight function
        weight[i] = eval<Dim>(geo, kernel, x, y, h);

        // check equivalence to generic interface
        vector<vector<double>> hvec(1,vector<double>(Dim));
        for (size_t j=0; j<Dim; j++) hvec[0][j] = h[j];
        double alt_wt = eval<Dim>(geo, kernel, x, y, hvec);
        EXPECT_NEAR(alt_wt, weight[i], 1.0e-12);
      }

      // Check that the weight at no other point is greater than that at x
      ASSERT_GE(weight_at_x, weight[i]);
      if (geo==FACETED) {
        ASSERT_NEAR(weight_at_x, 1.0, 1.e-12);
      }

      // Check that the weight is positive inside the support and
      // near zero or zero outside it

      Vector<Dim> v = y - x;
      Vector<Dim> h2;
      for (size_t d = 0; d < Dim; d++)
        h2[d] = h[d];

      bool outside = false;
      if (geo == ELLIPTIC) {
        double s = 0.0;
        for (size_t d = 0; d < Dim; d++)
          s += v[d] * v[d] / (2 * h[d] * 2 * h[d]);
        s = sqrt(s);
        outside = (s > 1.0);
      } else if (geo == TENSOR or geo == FACETED) {
        for (size_t d = 0; d < Dim; d++)
          if (fabs(v[d]) > 2 * h[d])
            outside = true;
      }

      if (outside)
        EXPECT_NEAR(0.0, weight[i], 1.0e-7);
      else
        ASSERT_GT(weight[i], 0.0);
    }
  }

  // Data

  Geometry geo;
  Kernel kernel;
  size_t np;
  vector<double> weights;
};

// Parameterized test for 1D
TEST_P(WeightTest, check_weights_1D) {
  checkWeight<1>(std::get<0>(GetParam()), std::get<1>(GetParam()));
}

// Parameterized test for 2D
TEST_P(WeightTest, check_weights_2D) {
  checkWeight<2>(std::get<0>(GetParam()), std::get<1>(GetParam()));
}

// Parameterized test for 3D
TEST_P(WeightTest, check_weights_3D) {
  checkWeight<3>(std::get<0>(GetParam()), std::get<1>(GetParam()));
}

INSTANTIATE_TEST_CASE_P(
    GeoKernelCombos,
    WeightTest,
    Combine(Values(ELLIPTIC, TENSOR),
            Values(B4, SQUARE, EPANECHNIKOV, INVSQRT, COULOMB)));

INSTANTIATE_TEST_CASE_P(FacetedSupport, WeightTest,
                        Combine(Values(FACETED), Values(POLYRAMP, STEP)));

