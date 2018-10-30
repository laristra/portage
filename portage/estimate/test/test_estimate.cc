/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <cmath>
#include <memory>
#include <vector>
#include <ostream>
#include "gtest/gtest.h"

#include "portage/estimate/estimate.h"

#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"
#include "portage/support/portage.h"

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
      {Point<dim> pt=(*src_pts)[i]; pt[k] = -.25 + 1.5*index*smoothing + delta; (*src_pts)[i] = pt;}
      {Point<dim> pt=(*sextents)[i]; pt[k] = 2.0*smoothing;
       if (index==0 or index == nside) pt[k]*=2; (*sextents)[i]=pt;}
    }
    extents = sextents;
  }
  for (size_t i=0; i<ntgt; i++) {
    size_t offset = 0, index;
    for (size_t k=0; k<dim; k++) {
      index = (i - offset)/powl(nside+3, dim-k-1);
      offset += index*powl(nside+3, dim-k-1);

      double delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter;
      {Point<dim> pt=(*tgt_pts)[i]; pt[k] = index*tsmoothing + delta; (*tgt_pts)[i]=pt;}

      {Point<dim> pt=(*textents)[i]; pt[k] = 2.0*tsmoothing;
       if (index==0 or index == nside+2) pt[k]*=2; (*textents)[i]=pt;}
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
  auto kernels = Portage::vector<Weight::Kernel>(nkern, Weight::B4);
  auto geometries = Portage::vector<Weight::Geometry>(nkern, Weight::TENSOR);
  auto smoothingh = Portage::vector<std::vector<std::vector<double>>>
      (nkern, std::vector<std::vector<double>>(1, std::vector<double>(dim)));
  for (size_t i=0; i<nkern; i++) for (size_t j=0; j<dim; j++) {
      std::vector<std::vector<double>> vv=smoothingh[i];
      Point<dim> pt = (*extents)[i];
      vv[0][j] = pt[j];
      smoothingh[i] = vv;
    }

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
  std::vector<unsigned int> src_particles(nsrc);
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

	// check the estimate, verifying reproducing property
	ASSERT_NEAR(result, jetx[j][k], 1.e-11);
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




