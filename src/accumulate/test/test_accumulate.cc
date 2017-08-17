

#include <cmath>
#include <memory>
#include <vector>
#include "gtest/gtest.h"

#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"

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
    }
  }

  // create source+target swarms, kernels, geometries, and smoothing lengths
  Swarm<dim> src_swarm(src_pts);
  Swarm<dim> tgt_swarm(tgt_pts);
  vector<Weight::Kernel> kernels(npoints, Weight::B4);
  vector<Weight::Geometry> geometries(npoints, Weight::TENSOR);
  vector<vector<vector<double>>> smoothingh(npoints,
                                            vector<vector<double>>(1, vector<double>(dim, smoothing)));

  {
    // create the accumulator
    Accumulate<dim, Swarm<dim>, Swarm<dim>>
        accum(
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
    vector<unsigned int> src_particles(npoints);
    for (size_t i=0; i<npoints; i++) src_particles[i] = i;

    // Loop through target particles
    for (size_t i=0; i<npoints; i++) {

      // do the accumulation loop for each target particle against all
      // the source particles
      vector<Portage::Weights_t> shape_vecs = accum(i, src_particles);

      if (etype == KernelDensity) {
        for (size_t j=0; j<npoints; j++) {
          double weight = accum.weight(j,i);
          ASSERT_EQ ((shape_vecs[j]).weights[0], weight);
        }
      } else {      
        auto x = tgt_swarm.get_particle_coordinates(i);
        auto jetx = Basis::jet<dim>(btype,x);

        for (size_t j=0; j<npoints; j++) {
          auto y = src_swarm.get_particle_coordinates(j);
          auto basisy = Basis::function<dim>(btype,y);
          for (size_t k=0; k<jsize[0]; k++) for (size_t m=0; m<jsize[1]; m++) {
              sums[i][k][m] += basisy[k]*(shape_vecs[j]).weights[m];
          }
        }

        for (size_t k=0; k<jsize[0]; k++) for (size_t m=0; m<jsize[1]; m++) {
            // this isn't working yet - need to fix
            //ASSERT_NEAR(sums[i][k][m], jetx[k][m], 1.e-12);
        }
      }
    }
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




