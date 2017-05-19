/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <cmath>
#include <ctime>

#include "gtest/gtest.h"

#include "portage/support/weight.h"
#include "portage/support/Point.h"
#include "portage/support/Vector.h"

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

using Portage::Meshfree::Weight::eval;

using Portage::Meshfree::Weight::FacetData;
using Portage::Meshfree::Weight::faceted;

using Portage::Point;
using Portage::Vector;
using std::tuple;
using std::array;
using std::vector;

using ::testing::TestWithParam;
using ::testing::Combine;
using ::testing::Values;

class WeightTest : public TestWithParam<tuple<Geometry, Kernel>> {
 public:
  template <size_t Dim>
      void checkWeight(const Geometry geo, const Kernel kernel) {

    //    srand(time(NULL));
    Point<Dim> x;
    for (int d = 0; d < Dim; d++) x[d] = ((double)rand())/RAND_MAX;
    
    static double nominal_h = .2;
    array<double, Dim> h;
    for (int d = 0; d < Dim; d++) h[d] = nominal_h;

    size_t nsides = 2*Dim;
    vector<FacetData<Dim>> facets(nsides);
    if (geo == FACETED) { 
      // set up an isothetic box like TENSOR
      size_t side = 0;
      for (size_t j=0; j<Dim; j++) {
	for (size_t k=0; k<Dim; k++) facets[side].normal[k] = 0.;
	facets[side].normal[j] = -1.;
	facets[side].smoothing = nominal_h;
	side++;
	for (size_t k=0; k<Dim; k++) facets[side].normal[k] = 0.;
	facets[side].normal[j] = 1.;
	facets[side].smoothing = nominal_h;
	side++;
      }
    }

    double weight_at_x;
    if (geo == FACETED) { 
      weight_at_x = faceted<Dim>(x,x, &facets[0], nsides);
    } else {
      weight_at_x = eval<Dim>(geo, kernel, x, x, h);
    }

    np = powl(10,Dim);
    vector<double> weight(np);

    for (size_t i = 0; i < np; i++) {
      // generate a 'y' val in the range (-3h[d],3h[d])  ----
      // rand/RAND_MAX will give us [0,1]; shift it left by
      // 0.5 and multiply by 6 to give us the desired interval

      Point<Dim> y = x;
      for (size_t d = 0; d < Dim; d++)
        y[d] += 6*h[d]*(((double)rand())/RAND_MAX - 0.5);
      
      if (geo == FACETED) {
	weight[i] = faceted<Dim>(x,y, &facets[0], nsides);
      } else {
	weight[i] = eval<Dim>(geo, kernel, x, y, h);
      }

      // Check that the weight at no other point is greater than that at x
      ASSERT_GE(weight_at_x, weight[i]);

      // Check that the weight is positive inside the support and
      // near zero or zero outside it

      Vector<Dim> v = y-x;
      Vector<Dim> h2;
      for (size_t d = 0; d < Dim; d++)
        h2[d] = h[d];

      bool outside = false;
      if (geo == ELLIPTIC) {
        double s = 0.0;
        for (size_t d = 0; d < Dim; d++)
          s += v[d]*v[d]/(2*h[d]*2*h[d]);
        s = sqrt(s);
        outside = (s > 1.0) ? true : false;
      } else if (geo == TENSOR or geo == FACETED) {
        for (size_t d = 0; d < Dim; d++)
          if (fabs(v[d]) > 2*h[d]) outside = true;
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

INSTANTIATE_TEST_CASE_P(GeoKernelCombos,
                        WeightTest,
                        Combine(Values(ELLIPTIC, TENSOR),
                                Values(B4, SQUARE, EPANECHNIKOV,
                                       INVSQRT, COULOMB)));

INSTANTIATE_TEST_CASE_P(FacetedSupport,
                        WeightTest,
                        Combine(Values(FACETED), Values(POLYRAMP)));
                        
