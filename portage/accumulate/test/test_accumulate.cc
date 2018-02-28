/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <cmath>
#include <memory>
#include <vector>
#include "gtest/gtest.h"

#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"
#include "portage/support/test_operator_data.h"

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
  const size_t nside = 6;
  const size_t npoints = powl(nside,dim);
  const double deltax = 1./nside;
  const double smoothing = 2.5*deltax;
  const double jitter = 0.2;
  auto src_pts = make_shared<typename Swarm<dim>::PointVec>(npoints);
  auto tgt_pts = make_shared<typename Swarm<dim>::PointVec>(npoints);
  for (size_t i=0; i<npoints; i++) {
    size_t offset = 0, index;
    for (size_t k=0; k<dim; k++) {
      index = (i - offset)/powl(nside, dim-k-1);
      offset += index*powl(nside, dim-k-1);

      double delta;

      delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter*deltax;
      (*src_pts)[i][k] = index*deltax + delta;

      delta = 2.*(((double)rand())/RAND_MAX -.5)*jitter*deltax;
      (*tgt_pts)[i][k] = index*deltax + delta;
    }
  }

  // create source+target swarms, kernels, geometries, and smoothing lengths
  Swarm<dim> src_swarm(src_pts);
  Swarm<dim> tgt_swarm(tgt_pts);
  vector<Weight::Kernel> kernels(npoints, Weight::B4);
  vector<Weight::Geometry> geometries(npoints, Weight::TENSOR);
  vector<vector<vector<double>>> 
    smoothingh(npoints, vector<vector<double>>(1, vector<double>(dim, smoothing)));

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

    // list of src swarm particles (indices)
    vector<unsigned int> src_particles(npoints);
    for (size_t i=0; i<npoints; i++) src_particles[i] = i;

    // Loop through target particles
    for (size_t i=0; i<npoints; i++) {
      vector<vector<double>> sums(jsize[0], vector<double>(jsize[1],0.));

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
              sums[k][m] += basisy[k]*(shape_vecs[j]).weights[m];
          }
        }

        for (size_t k=0; k<jsize[0]; k++) for (size_t m=0; m<jsize[1]; m++) {
            ASSERT_NEAR(sums[k][m], jetx[k][m], 1.e-11);
        }
      }
    }
  }
}


// Test the reproducing property of basis integration operators
template<
Portage::Meshfree::Basis::Type btype, 
Portage::Meshfree::Operator::Type opertype, 
Portage::Meshfree::Operator::Domain domain
>
void test_operator(Portage::Meshfree::WeightCenter center) {

  using namespace Portage::Meshfree;

  // create the source swarm input geometry data
  constexpr size_t dim=Operator::dimension(domain);
  const size_t nside = 6;
  const size_t npoints = powl(nside,dim);
  const double deltax = 1./nside;
  const double jitter = 0.3;
  const double smoothing = 2.5*(1.+jitter)*deltax;
  auto src_pts = make_shared<typename Swarm<dim>::PointVec>(npoints);
  for (size_t i=0; i<npoints; i++) {
    size_t offset = 0, index;
    for (size_t k=0; k<dim; k++) {
      index = (i - offset)/powl(nside, dim-k-1);
      offset += index*powl(nside, dim-k-1);

      double delta;

      delta = (((double)rand())/RAND_MAX)*jitter*deltax;
      (*src_pts)[i][k] = index*deltax + delta;
    }
  } 

  // create the target swarm input geometry data based on integration domains.
  auto tgt_pts = make_shared<typename Swarm<dim>::PointVec>(1);
  vector<vector<Point<dim>>> domain_points(1);
  domain_points[0] = reference_points<domain>();
  Point<dim> centroid(vector<double>(dim,0.));
  for (int i=0; i<domain_points[0].size(); i++) {
    for (int j=0; j<dim; j++) centroid[j] += domain_points[0][i][j];
  }
  for (int k=0; k<dim; k++) {
    centroid[k] /= domain_points[0].size();
    (*tgt_pts)[0][k] = centroid[k];
  }

  // create source+target swarms, kernels, geometries, and smoothing lengths
  Swarm<dim> src_swarm(src_pts);
  Swarm<dim> tgt_swarm(tgt_pts);
  vector<Weight::Kernel> kernels(npoints, Weight::B4);
  vector<Weight::Geometry> geometries(npoints, Weight::TENSOR);
  vector<vector<vector<double>>> 
    smoothingh(npoints, vector<vector<double>>(1, vector<double>(dim, smoothing)));

  // create the accumulator
  Accumulate<dim, Swarm<dim>, Swarm<dim>>
    accum(
	  src_swarm,
	  tgt_swarm,
	  OperatorRegression,
	  center,
	  kernels,
	  geometries,
	  smoothingh,
	  btype, 
	  opertype,
	  domain, 
	  domain_points);

  // check sizes and allocate test array
  size_t bsize = Basis::function_size<dim>(btype);
  size_t jsize = Operator::Operator<opertype, btype, domain>::operator_size;

  // list of src swarm particles (indices)
  vector<unsigned int> src_particles(npoints);
  for (size_t i=0; i<npoints; i++) src_particles[i] = i;

  vector<vector<double>> sums(bsize, vector<double>(jsize,0.));

  // do the accumulation loop for each target particle against all
  // the source particles
  vector<Portage::Weights_t> shape_vecs = accum(0, src_particles);

  for (size_t j=0; j<npoints; j++) {
    auto y = src_swarm.get_particle_coordinates(j);
    auto basisy = Basis::function<dim>(btype,y);
    for (size_t k=0; k<bsize; k++) 
      for (size_t m=0; m<jsize; m++) {
	sums[k][m] += basisy[k]*(shape_vecs[j]).weights[m];
      }
  }

  for (size_t k=0; k<bsize; k++) { 
    for (size_t m=0; m<jsize; m++) {
      if (opertype == Operator::VolumeIntegral) {
	ASSERT_EQ(jsize, 1);
	switch (btype) {
	case Basis::Unitary: {
	  switch (domain) {
	  case Operator::Interval:      {ASSERT_NEAR(sums[k][m],exactUnitaryInterval[k], 1.e-12);     break;}
	  case Operator::Quadrilateral: {ASSERT_NEAR(sums[k][m],exactUnitaryQuadrilateral[k], 1.e-12);break;}
	  case Operator::Triangle:      {ASSERT_NEAR(sums[k][m],exactUnitaryTriangle[k], 1.e-12);     break;}
	  case Operator::Hexahedron:    {ASSERT_NEAR(sums[k][m],exactUnitaryHexahedron[k], 1.e-12);   break;}
	  case Operator::Wedge:         {ASSERT_NEAR(sums[k][m],exactUnitaryWedge[k], 1.e-12);        break;}
	  case Operator::Tetrahedron:   {ASSERT_NEAR(sums[k][m],exactUnitaryTetrahedron[k], 1.e-12);  break;}
	  default:
	    break;
	  }
	}
	case Basis::Linear: {
	  switch (domain) {
	  case Operator::Interval:      {ASSERT_NEAR(sums[k][m],exactLinearInterval[k], 1.e-12);     break;}
	  case Operator::Quadrilateral: {ASSERT_NEAR(sums[k][m],exactLinearQuadrilateral[k], 1.e-12);break;}
	  case Operator::Triangle:      {ASSERT_NEAR(sums[k][m],exactLinearTriangle[k], 1.e-12);     break;}
	  case Operator::Hexahedron:    {ASSERT_NEAR(sums[k][m],exactLinearHexahedron[k], 1.e-12);   break;}
	  case Operator::Wedge:         {ASSERT_NEAR(sums[k][m],exactLinearWedge[k], 1.e-12);        break;}
	  case Operator::Tetrahedron:   {ASSERT_NEAR(sums[k][m],exactLinearTetrahedron[k], 1.e-12);  break;}
	  default:
	    break;
	  }
	}
	case Basis::Quadratic: {
	  switch (domain) {
	  case Operator::Interval:      {ASSERT_NEAR(sums[k][m],exactQuadraticInterval[k], 1.e-12);     break;}
	  case Operator::Quadrilateral: {ASSERT_NEAR(sums[k][m],exactQuadraticQuadrilateral[k], 1.e-12);break;}
	  case Operator::Triangle:      {ASSERT_NEAR(sums[k][m],exactQuadraticTriangle[k], 1.e-12);     break;}
	  case Operator::Tetrahedron:   {ASSERT_NEAR(sums[k][m],exactQuadraticTetrahedron[k], 1.e-12);  break;}
	  default:
	    break;
	  }
	}
	default:
	  break;
	}
      }
    }
  }
}

// test the pointwise estimation capability

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

// test the operator capability

TEST(operator, UnitaryInterval) {
  test_operator<Portage::Meshfree::Basis::Type::Unitary,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Interval>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, UnitaryQuadrilateral) {
  test_operator<Portage::Meshfree::Basis::Type::Unitary,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Quadrilateral>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, UnitaryTriangle) {
  test_operator<Portage::Meshfree::Basis::Type::Unitary,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Triangle>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, UnitaryHexahedron) {
  test_operator<Portage::Meshfree::Basis::Type::Unitary,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Hexahedron>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, UnitaryWedge) {
  test_operator<Portage::Meshfree::Basis::Type::Unitary,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Wedge>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, UnitaryTetrahedron) {
  test_operator<Portage::Meshfree::Basis::Type::Unitary,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Tetrahedron>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, LinearInterval) {
  test_operator<Portage::Meshfree::Basis::Type::Linear,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Interval>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, LinearQuadrilateral) {
  test_operator<Portage::Meshfree::Basis::Type::Linear,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Quadrilateral>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, LinearTriangle) {
  test_operator<Portage::Meshfree::Basis::Type::Linear,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Triangle>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, LinearHexahedron) {
  test_operator<Portage::Meshfree::Basis::Type::Linear,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Hexahedron>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, LinearWedge) {
  test_operator<Portage::Meshfree::Basis::Type::Linear,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Wedge>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, LinearTetrahedron) {
  test_operator<Portage::Meshfree::Basis::Type::Linear,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Tetrahedron>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, QuadraticInterval) {
  test_operator<Portage::Meshfree::Basis::Type::Quadratic,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Interval>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, QuadraticQuadrilateral) {
  test_operator<Portage::Meshfree::Basis::Type::Quadratic,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Quadrilateral>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, QuadraticTriangle) {
  test_operator<Portage::Meshfree::Basis::Type::Quadratic,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Triangle>
    (Portage::Meshfree::WeightCenter::Scatter);
}

TEST(operator, QuadraticTetrahedron) {
  test_operator<Portage::Meshfree::Basis::Type::Quadratic,
                Portage::Meshfree::Operator::VolumeIntegral, 
		Portage::Meshfree::Operator::Tetrahedron>
    (Portage::Meshfree::WeightCenter::Scatter);
}




