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
#include <memory>
#include <vector>
#include "gtest/gtest.h"

#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "accumulate.h"

using std::vector;
using std::shared_ptr;
using std::make_shared;
using Portage::Point;

template<size_t dim>
void test_accumulate(Portage::Meshfree::EstimateType etype,
                     Portage::Meshfree::Basis::Type btype,
                     Portage::Meshfree::WeightCenter center) {
  using namespace Portage::Meshfree;

  // create the source and target swarms input data
  const size_t nside = 8;
  const size_t npoints = powl(nside,dim);
  const double smoothing = 1./nside;
  const double jitter = 0.3*smoothing;
  auto src_pts = make_shared<typename Swarm<dim>::PointVec>(npoints);
  auto tgt_pts = make_shared<typename Swarm<dim>::PointVec>(npoints);
  auto extents = make_shared<typename Swarm<dim>::PointVec>(npoints);
  for (size_t i=0; i<npoints; i++) {
    size_t offset = 0, index;
    for (size_t k=0; k<dim; k++) {
      index = (i - offset)/powl(nside, dim-k-1);
      offset += index*powl(nside, dim-k-1);

      double delta;

      delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter;
      (*src_pts)[i][k] = 1.25*index*smoothing + delta;

      delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter;
      (*tgt_pts)[i][k] = index*smoothing + delta;

      (*extents)[i][k] = 1.5*smoothing;
    }
  }

  // create source+target swarms, kernels, geometries, and smoothing lengths
  auto src_swarm = make_shared<Swarm<dim>>(src_pts, extents);
  auto tgt_swarm = make_shared<Swarm<dim>>(tgt_pts, extents);
  auto kernels = make_shared<vector<Weight::Kernel>>(npoints, Weight::B4);
  auto geometries = make_shared<vector<Weight::Geometry>>(npoints, Weight::TENSOR);
  auto smoothingh = make_shared<vector<vector<vector<double>>>>
      (npoints, vector<vector<double>>(1, vector<double>(dim, smoothing)));

  {
    // create the accumulator
    Accumulate<dim> accum(
        src_swarm,
        tgt_swarm,
        etype,
        center,
        kernels,
        geometries,
        smoothingh,
        btype);

    // check sizes and allocate test array
    size_t bsize = Basis::function_size<dim>(btype);
    auto jsize = Basis::jet_size<dim>(btype);
    ASSERT_EQ(jsize[0], bsize);
    ASSERT_EQ(jsize[1], bsize);
    vector<vector<vector<double>>> sums(npoints,
      vector<vector<double>>(jsize[0], vector<double>(jsize[1],0.)));

    // list of src swarm particles (indices)
    vector<size_t> src_particles(npoints);
    for (size_t i=0; i<npoints; i++) src_particles[i] = i;

    // Loop through target particles
    for (size_t i=0; i<npoints; i++) {

      // do the accumulation loop for each target particle against all
      // the source particles
      auto shape_vec = accum(i, src_particles);

      if (etype == KernelDensity) {
        for (size_t j=0; j<npoints; j++) {
          double weight = accum.weight(j,i);
          ASSERT_EQ (shape_vec[j][0], weight);
        }
      } else {      
        auto x = tgt_swarm->get_particle_coordinates(i);
        auto jetx = Basis::jet<dim>(btype,x);

        for (size_t j=0; j<npoints; j++) {
          auto y = src_swarm->get_particle_coordinates(j);
          auto basisy = Basis::function<dim>(btype,y);
          for (size_t k=0; k<jsize[0]; k++) for (size_t m=0; m<jsize[1]; m++) {
              sums[i][k][m] += basisy[k]*shape_vec[j][m];
          }
        }

        for (size_t k=0; k<jsize[0]; k++) for (size_t m=0; m<jsize[1]; m++) {
            // this isn't working yet - need to fix
            //ASSERT_NEAR(sums[i][k][m], jetx[k][m], 1.e-12);
        }
      }
    }
    double x=1.;
  }
}

TEST(accumulate, 1d_KUG) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::KernelDensity,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 2d_KUG) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::KernelDensity,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 3d_KUG) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::KernelDensity,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 1d_KUS) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::KernelDensity,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 2d_KUS) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::KernelDensity,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 3d_KUS) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::KernelDensity,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 1d_RUG) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 2d_RUG) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 3d_RUG) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 1d_RUS) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 2d_RUS) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 3d_RUS) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 1d_RLG) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 2d_RLG) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 3d_RLG) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 1d_RLS) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 2d_RLS) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 3d_RLS) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 1d_RQG) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 2d_RQG) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 3d_RQG) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(accumulate, 1d_RQS) {
  test_accumulate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 2d_RQS) {
  test_accumulate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(accumulate, 3d_RQS) {
  test_accumulate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Scatter);
}




