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
#include <ostream>
#include "gtest/gtest.h"

#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"
#include "estimate.h"

using std::vector;
using std::shared_ptr;
using std::make_shared;
using std::cout;
using Portage::Point;

template<size_t dim>
void test_estimate(Portage::Meshfree::EstimateType etype,
                   Portage::Meshfree::Basis::Type btype,
                   Portage::Meshfree::WeightCenter center) {
  using namespace Portage::Meshfree;

  // create the source and target swarms input data
  const size_t nside = 3;
  const double smoothing = 1./nside;
  const double jitter = 0.;  // 0.25*smoothing;
  const size_t nsrc = powl(nside+1,dim);
  const size_t ntgt = powl(nside+3,dim);
  const double tsmoothing = 1./(nside+2);
  auto src_pts = make_shared<typename Swarm<dim>::PointVec>(nsrc);
  auto tgt_pts = make_shared<typename Swarm<dim>::PointVec>(ntgt);
  auto sextents = make_shared<typename Swarm<dim>::PointVec>(nsrc);
  auto textents = make_shared<typename Swarm<dim>::PointVec>(ntgt);
  auto extents = make_shared<typename Swarm<dim>::PointVec>(ntgt);
  for (size_t i=0; i<nsrc; i++) {
    size_t offset = 0, index;
    for (size_t k=0; k<dim; k++) {
      index = (i - offset)/powl(nside+1, dim-k-1);
      offset += index*powl(nside+1, dim-k-1);

      double delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter*1.5;
      (*src_pts)[i][k] = -.25 + 1.5*index*smoothing + delta;
      (*sextents)[i][k] = 2.0*smoothing;
      if (index==0 or index == nside) (*sextents)[i][k]*=2;
    }
    extents = sextents;
  }
  for (size_t i=0; i<ntgt; i++) {
    size_t offset = 0, index;
    for (size_t k=0; k<dim; k++) {
      index = (i - offset)/powl(nside+3, dim-k-1);
      offset += index*powl(nside+3, dim-k-1);

      double delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter;
      (*tgt_pts)[i][k] = index*tsmoothing + delta;

      (*textents)[i][k] = 2.0*tsmoothing;
      if (index==0 or index == nside+2) (*textents)[i][k]*=2;
    }
    extents = textents;
  }

  // create source+target swarms, kernels, geometries, and smoothing lengths
  auto src_swarm = Swarm<dim>(src_pts);
  auto src_state = SwarmState<dim>(src_swarm);
  auto tgt_swarm = Swarm<dim>(tgt_pts);
  auto tgt_state = SwarmState<dim>(tgt_swarm);
  size_t nkern;
  if (center == Gather) nkern = ntgt; else if (center == Scatter) nkern = nsrc;
  auto kernels = vector<Weight::Kernel>(nkern, Weight::B4);
  auto geometries = vector<Weight::Geometry>(nkern, Weight::TENSOR);
  auto smoothingh = vector<vector<vector<double>>>
      (nkern, vector<vector<double>>(1, vector<double>(dim)));
  for (size_t i=0; i<nkern; i++) for (size_t j=0; j<dim; j++) smoothingh[i][0][j] = (*extents)[i][j];

  // create the accumulator
  Accumulate<dim, Swarm<dim>, Swarm<dim>> accum(
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

  // list of src swarm particles (indices)
  vector<unsigned int> src_particles(nsrc);
  for (unsigned int i=0; i<nsrc; i++) src_particles[i] = i;

  // add fields
  string cnums[15]={"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14"};
  string fnames[15];
  size_t nbasis=jsize[0]; 
  for (size_t i=0; i<nbasis; i++) fnames[i] = string("field")+cnums[i];
  for (size_t i=0; i<nbasis; i++) {
    typename SwarmState<dim>::DblVecPtr sfieldp=make_shared<typename SwarmState<dim>::DblVec>(nsrc,0.);
    src_state.add_field(fnames[i], sfieldp);
    typename SwarmState<dim>::DblVecPtr tfieldp=make_shared<typename SwarmState<dim>::DblVec>(ntgt,0.);
    tgt_state.add_field(fnames[i], tfieldp);
  }

  // fill in source state
  for (size_t j=0; j<nbasis; j++) {
    typename SwarmState<dim>::DblVecPtr fieldp; src_state.get_field(fnames[j], fieldp);
    auto &field(*fieldp);
    for (size_t i=0; i<nsrc; i++) {
      auto y = src_swarm.get_particle_coordinates(i);
      auto basisy = Basis::function<dim>(btype,y);
      field[i] = basisy[j];
    }
  }

  // make the estimator
  Estimate<dim, SwarmState<dim>> estimate(src_state);
    
  // Loop through target particles
  for (size_t i=0; i<ntgt; i++) {
    auto x = tgt_swarm.get_particle_coordinates(i);
    auto jetx = Basis::jet<dim>(btype,x);
    auto sources_and_mults = accum(i, src_particles);

    // count actual neighbors
    size_t nnbr=0;
    for (size_t k=0; k<nsrc; k++) if (accum.weight(k,i)!=0.) nnbr++;
    if (nnbr < nbasis) {
      std::cout << "number of neighbors "<< nnbr << " is too small at " << i << "\n";
    }

    // Loop through fields
    for (size_t j=0; j<nbasis; j++) {

      // get the target field
      typename SwarmState<dim>::DblVecPtr fieldp; tgt_state.get_field(fnames[j], fieldp);
      auto &field(*fieldp);

      // Loop through derivatives
      for (size_t k=0; k<nbasis; k++) {
	// tell estimate what to do
	estimate.set_variable(fnames[j], k);

	// do the estimate
	double result = estimate(i, sources_and_mults);

	// save value
	if (k==0) field[i] = result;

	// check the estimate
	ASSERT_NEAR(result, jetx[j][k], 1.e-12);
      }
    }
  }
}

TEST(estimate, 1d_RUG) {
  test_estimate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 2d_RUG) {
  test_estimate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 3d_RUG) {
  test_estimate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 1d_RUS) {
  test_estimate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 2d_RUS) {
  test_estimate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 3d_RUS) {
  test_estimate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Unitary,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 1d_RLG) {
  test_estimate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 2d_RLG) {
  test_estimate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 3d_RLG) {
  test_estimate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 1d_RLS) {
  test_estimate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 2d_RLS) {
  test_estimate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 3d_RLS) {
  test_estimate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Linear,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 1d_RQG) {
  test_estimate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 2d_RQG) {
  test_estimate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 3d_RQG) {
  test_estimate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Gather);
}

TEST(estimate, 1d_RQS) {
  test_estimate<1>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 2d_RQS) {
  test_estimate<2>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Scatter);
}

TEST(estimate, 3d_RQS) {
  test_estimate<3>(Portage::Meshfree::EstimateType::LocalRegression,
                     Portage::Meshfree::Basis::Type::Quadratic,
                     Portage::Meshfree::WeightCenter::Scatter);
}




