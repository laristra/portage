/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <cassert>
#include <cmath>
#include <string>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/driver_swarm.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/operator.h"
#include "portage/support/faceted_setup.h"

#include "portage/support/portage.h"
#include "portage/search/search_points_by_cells.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

namespace {

using std::shared_ptr;
using std::make_shared;
using Wonton::Point;

using namespace Portage::Meshfree;

double TOL = 1e-6;

// This is a set of integration tests for the swarm-swarm remap driver.
// There will be at least one test corresponding to each case found in
// main.cc. This is a test fixture and must be derived from the
// ::testing::Test class. Specializations of this class, such as 2D/3D
// coincident and non-coincident remaps should be derived from this.

template<size_t dim>
class BaseTest : public ::testing::Test {
protected:
  // swarms and states
  Swarm<dim> source_swarm;
  Swarm<dim> target_swarm;
  SwarmState<dim> source_state;
  SwarmState<dim> target_state;

  // smoothing lengths
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths_;

  // kernel and geometry specifications
  Portage::vector<Weight::Kernel> kernels_;
  Portage::vector<Weight::Geometry> geometries_;

  // operator info
  WeightCenter center_;
  Operator::Type operator_;
  Portage::vector<Operator::Domain> domains_;
  Portage::vector<std::vector<Point<dim>>> operator_data_;

  // Constructor for Driver test
  BaseTest(int nb_source, int nb_target, int distrib,
           double x_min = 0.0, double x_max = 0.0,
           double y_min = 0.0, double y_max = 0.0,
           double z_min = 0.0, double z_max = 0.0,
           Operator::Type op = Operator::LastOperator)
    : source_swarm(nb_source, distrib, 0,
                   x_min, x_max, y_min, y_max, z_min, z_max),
      target_swarm(nb_target, distrib, 0,
                   x_min, x_max, y_min, y_max, z_min, z_max),
      source_state(source_swarm),
      target_state(target_swarm),
      center_(Gather),
      operator_(op)
  {

    if (op != Operator::LastOperator) {

      int const num_owned_target = target_swarm.num_owned_particles();

      domains_ = Portage::vector<Operator::Domain>(num_owned_target);
      size_t npoints[3]={2,4,8};
      operator_data_.resize(num_owned_target, std::vector<Point<dim>>(npoints[dim - 1]));


      size_t npdim = static_cast<size_t>(pow(1.001*operator_data_.size(),1./dim));
      size_t nptot = 1; for (size_t m=0; m<dim; m++) nptot *= npdim;
      assert(nptot == target_swarm->num_owned_particles());
      assert(npdim>1);
      size_t n=0;
      size_t ij[npoints[dim-1]]; 
      size_t i=0,j=0,k=0;
      Point<dim> pt;
      std::vector<Point<dim>> points(npoints[dim-1]);

      // Create points determining the integration volume.
      // Assumes target swarm is created from SwarmFactory and represents a perfect cubic array of points.
      for (size_t n=0; n<nptot; n++) {
        switch(dim) {
        case 1:
          i = n;
          ij[0]=i; ij[1]=(i+1);
          break;
        case 2:
          i = n/npdim;
          j = n - i*npdim;
          ij[0]=i*npdim+j; ij[1]=(i+1)*npdim+j; ij[2]=(i+1)*npdim+j+1; ij[3]=i*npdim+j+1;
          break;
        case 3:
          i = n/(npdim*npdim);
          j = (n - i*npdim*npdim)/npdim;
          k = n - i*npdim*npdim - j*npdim;
          ij[0]=npdim*(i*npdim+j)+k;         ij[1]=npdim*((i+1)*npdim+j)+k;
          ij[2]=npdim*((i+1)*npdim+j+1)+k;   ij[3]=npdim*(i*npdim+j+1)+k;
          ij[4]=npdim*(i*npdim+j)+k+1;       ij[5]=npdim*((i+1)*npdim+j)+k+1;
          ij[6]=npdim*((i+1)*npdim+j+1)+k+1; ij[7]=npdim*(i*npdim+j+1)+k+1;
          break;
          default: break;
        }

        for (int m=0; m<npoints[dim-1]; m++) {
          if (i==npdim-1 or j==npdim-1 or k==npdim-1) {
            // particles on the upper boundaries get an integration volume of zero
            pt = target_swarm.get_particle_coordinates(ij[0]);
          } else {
            pt = target_swarm.get_particle_coordinates(ij[m]);
          }
          points[m] = pt;
        }
        operator_data_[n] = points;
        domains_[n] = Operator::domain_from_points<dim>(points);
      }
    }
  }

  /*void set_smoothing_lengths(shared_ptr<Portage::vector<std::vector<std::vector<double>>>> smoothing_lengths,
			     Portage::Meshfree::WeightCenter center=Portage::Meshfree::Gather) {
    smoothing_lengths_ = smoothing_lengths;
    center_ = center;
  }*/

  void set_smoothing_lengths(const int* n, double h, WeightCenter center = Gather) {
    assert(n != nullptr);
    std::vector<std::vector<double>> const default_lengths(n[1], std::vector<double>(n[2], h));
    smoothing_lengths_.resize(n[0], default_lengths);
    center_ = center;
  }

  // This is the basic test method to be called for each unit test.
  // It will work for 1, 2-D and 3-D swarms
  //
  template <template<int, class, class> class Search,
            Basis::Type basis>
  void unitTest(double compute_initial_field(Wonton::Point<dim> coord),
                double expected_answer) {

    // Fill the source state data with the specified profile
    const int nsrcpts = source_swarm.num_owned_particles();
    const int ntarpts = target_swarm.num_owned_particles();
//    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr sourceData =
//        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(nsrcpts);

    Portage::vector<double> sourceData(nsrcpts);
    Portage::vector<double> targetData(nsrcpts, 0.0);

    // Create the source data for given function
    for (unsigned int p = 0; p < nsrcpts; ++p) {
      auto coord = source_swarm.get_particle_coordinates(p);
      sourceData[p] = compute_initial_field(coord);
    }

    source_state.add_field("particledata", sourceData);
    target_state.add_field("particledata", targetData);

    // Build the main driver data for this mesh type
    SwarmDriver<Search, Accumulate, Estimate, dim,
                Swarm<dim>, SwarmState<dim>,
                Swarm<dim>, SwarmState<dim>>
        d(source_swarm, source_state,
          target_swarm, target_state,
          smoothing_lengths_,
          Weight::B4, Weight::ELLIPTIC, center_);

    EstimateType estimator=LocalRegression;
    if (operator_ != Portage::Meshfree::Operator::LastOperator) 
      estimator = Portage::Meshfree::OperatorRegression;

    // Register the variable name with the driver
    std::vector<std::string> remap_fields;
    remap_fields.emplace_back("particledata");
    d.set_remap_var_names(remap_fields, remap_fields,
                          estimator, basis, 
                          operator_, domains_, operator_data_);

    // run on one processor (no argument implies serial run)
    d.run();

    // Check the answer
    double toterr=0.;
    //typename Portage::Meshfree::SwarmState<dim>::DblVecPtr vecout;
    auto vecout = target_state.get_field("particledata");

    ASSERT_NE(nullptr, vecout);
    if (operator_ == Portage::Meshfree::Operator::LastOperator) {
      for (int p = 0; p < ntarpts; ++p) {
        Wonton::Point<dim> coord = target_swarm.get_particle_coordinates(p);
        double error;
        error = compute_initial_field(coord) - vecout[p];
        // dump diagnostics for each particle
        if (dim == 1)
          std::printf("Particle=% 4d Coord = (% 5.3lf)", p, coord[0]);
        else if (dim == 2)
          std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf)", p, coord[0],
                      coord[1]);
        else if (dim == 3)
          std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf,% 5.3lf)", p,
                      coord[0], coord[1], coord[2]);
	{double val=vecout[p];
	  std::printf("  Value = % 10.6lf  Err = % lf\n", val, error);}
        toterr += error*error;
      }
    
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
      ASSERT_NEAR(expected_answer, sqrt(toterr), TOL);
    } else if (operator_ == Portage::Meshfree::Operator::VolumeIntegral) {
      double total = 0.;
      for (int p = 0; p < ntarpts; ++p) {
        total += vecout[p];
      }
      ASSERT_NEAR(expected_answer, total, TOL);
    }
  }



  // This unit test exercises the alternate more detailed constructor.
  // It will work for 1, 2-D and 3-D swarms
  //
  template <template<int, class, class> class Search,
            Portage::Meshfree::Basis::Type basis>
  void unitTestAlt(double compute_initial_field(Wonton::Point<dim> coord),
                double expected_answer) {

    // Fill the source state data with the specified profile
    const int nsrcpts = source_swarm->num_owned_particles();
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr sourceData = 
        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(nsrcpts);

    // Create the source data for given function
    for (unsigned int p = 0; p < nsrcpts; ++p) {
      Wonton::Point<dim> coord =
          source_swarm->get_particle_coordinates(p);
      (*sourceData)[p] = compute_initial_field(coord);
    }
    source_state->add_field("particledata", sourceData);

    // Build the target state storage
    const int ntarpts = target_swarm->num_owned_particles();
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr targetData = 
        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(ntarpts, 0.0);
    target_state->add_field("particledata", targetData);

    // Set kernels and geometries
    if (center_ == Portage::Meshfree::Gather) {
      kernels_.resize(ntarpts, Portage::Meshfree::Weight::B4);
      geometries_.resize(ntarpts, Portage::Meshfree::Weight::ELLIPTIC); 
    } else if (center_ == Portage::Meshfree::Scatter) {
      kernels_.resize(nsrcpts, Portage::Meshfree::Weight::B4);
      geometries_.resize(nsrcpts, Portage::Meshfree::Weight::ELLIPTIC); 
    }

    // Build the main driver data for this mesh type
    Portage::Meshfree::SwarmDriver<Search,
                                   Portage::Meshfree::Accumulate,
                                   Portage::Meshfree::Estimate,
                                   dim,
                                   Portage::Meshfree::Swarm<dim>,
                                   Portage::Meshfree::SwarmState<dim>,
                                   Portage::Meshfree::Swarm<dim>,
                                   Portage::Meshfree::SwarmState<dim>>
        d(*source_swarm, *source_state, *target_swarm, *target_state, *smoothing_lengths_,
          kernels_, geometries_, center_);

    Portage::Meshfree::EstimateType estimator=Portage::Meshfree::LocalRegression;
    if (operator_ != Portage::Meshfree::Operator::LastOperator) 
      estimator = Portage::Meshfree::OperatorRegression;

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("particledata");
    d.set_remap_var_names(remap_fields, remap_fields,
                          estimator, basis, 
                          operator_, domains_, operator_data_);

    // run on one processor (no argument implies serial run)
    d.run();

    // Check the answer
    double toterr=0.;
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr vecout;
    target_state->copy_field("particledata", vecout);
    ASSERT_NE(nullptr, vecout);
    if (operator_ == Portage::Meshfree::Operator::LastOperator) {
      for (int p = 0; p < ntarpts; ++p) {
        Wonton::Point<dim> coord = target_swarm->get_particle_coordinates(p);
        double error;
        error = compute_initial_field(coord) - (*vecout)[p];
        // dump diagnostics for each particle
        if (dim == 1)
          std::printf("Particle=% 4d Coord = (%10.3e)", p, coord[0]);
        else if (dim == 2)
          std::printf("Particle=% 4d Coord = (%10.3e,%10.3e)", p, coord[0],
                      coord[1]);
        else if (dim == 3)
          std::printf("Particle=% 4d Coord = (%10.3e,%10.3e,%10.3e)", p,
                      coord[0], coord[1], coord[2]);
	{double val=(*vecout)[p]; 
	  std::printf("  Value = %10.3e  Err = %10.3e\n", val, error);}
        toterr += error*error;
      }
    
      std::printf("\n\nL2 NORM OF ERROR = %e10.3\n\n", sqrt(toterr));
      ASSERT_NEAR(expected_answer, sqrt(toterr/ntarpts), TOL);
    } else if (operator_ == Portage::Meshfree::Operator::VolumeIntegral) {
      double total = 0.;
      for (int p = 0; p < ntarpts; ++p) {
        total += (*vecout)[p];
      }
      ASSERT_NEAR(expected_answer, total, TOL);
    }
  }

};

//  std::shared_ptr<Swarm<2>> SwarmFactory(double xmin, double ymin,
//                                         double xmax, double ymax,
//                                         unsigned int nparticles,
//                                         unsigned int distribution,
//                                         unsigned int rand_seed=0)

//  Swarm(int num_particles, int distribution,
//        unsigned user_seed = 0,
//        double x_min = 0.0, double x_max = 0.0,
//        double y_min = 0.0, double y_max = 0.0,
//        double z_min = 0.0, double z_max = 0.0);

// Class which constructs a pair of 1-D swarms (random distribution) for remaps
class DriverTest1DGather : public BaseTest<1> {
public:
  DriverTest1DGather() : BaseTest<1>(7, 5, 2, 0.0, 1.0) {
    int const dim[] = { 5, 1, 1 };
    BaseTest<1>::set_smoothing_lengths(dim, 0.5, Gather);
  }
};

class DriverTest2DGather : public BaseTest<2> {
public:
  DriverTest2DGather() : BaseTest<2>(7*7, 5*5, 2, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 5*5, 1, 2 };
    BaseTest<2>::set_smoothing_lengths(dim, 0.5, Gather);
  }
};

class DriverTest3DGather : public BaseTest<3> {
public:
  DriverTest3DGather() : BaseTest<3>(7*7*7, 5*5*5, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 5*5*5, 1, 3 };
    BaseTest<3>::set_smoothing_lengths(dim, 0.5, Gather);
  }
};

class DriverTest1DScatter : public BaseTest<1> {
public:
  DriverTest1DScatter() : BaseTest<1>(7, 5, 2, 0.0, 1.0) {
    int const dim[] = { 5, 1, 1 };
    BaseTest<1>::set_smoothing_lengths(dim, 0.5, Scatter);
  }
};

class DriverTest2DScatter : public BaseTest<2> {
public:
  DriverTest2DScatter() : BaseTest<2>(7*7, 5*5, 2, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 5*5, 1, 2 };
    BaseTest<2>::set_smoothing_lengths(dim, 0.5, Scatter);
  }
};

class DriverTest3DScatter : public BaseTest<3> {
public:
  DriverTest3DScatter() : BaseTest<3>(7*7*7, 5*5*5, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 5*5*5, 1, 3 };
    BaseTest<3>::set_smoothing_lengths(dim, 0.5, Scatter);
  }
};

// Class which constructs a pair of 1-D swarms (ordered distribution) for remaps, and integrates
class IntegrateDriverTest1D : public BaseTest<1> {
public:
  IntegrateDriverTest1D() : BaseTest<1>(7, 5, 1, 0.0, 1.0) {
    int const dim[] = { 5, 1, 1 };
    BaseTest<1>::set_smoothing_lengths(dim, 0.5);
  }
};

class IntegrateDriverTest2D : public BaseTest<2> {
public:
  IntegrateDriverTest2D() : BaseTest<2>(7*7, 5*5, 1, 0.0, 1.0) {
    int const dim[] = { 5*5, 1, 2 };
    BaseTest<2>::set_smoothing_lengths(dim, 0.5);
  }
};

class IntegrateDriverTest3D : public BaseTest<3> {
public:
  IntegrateDriverTest3D() : BaseTest<3>(7*7*7, 5*5*5, 1, 0.0, 1.0) {
    int const dim[] = { 5*5*5, 1, 3 };
    BaseTest<3>::set_smoothing_lengths(dim, 0.5);
  }
};

template<int dim>
double compute_constant_field(Wonton::Point<dim> coord) {
  return 25.0;
}

// Methods for computing initial field values
template<int dim>
double compute_linear_field(Wonton::Point<dim> coord) {
  double val = 0.0;
  for (int i = 0; i < dim; i++)
    val += coord[i];
  return val;
}

template<int dim>
double compute_quadratic_field(Wonton::Point<dim> coord) {
  double val = 0.0;
  for (int i = 0; i < dim; i++)
    for (int j = i; j < dim; j++)
      val += coord[i] * coord[j];
  return val;
}

template<size_t dim>
double compute_cubic_field(Wonton::Point<dim> coord) {
  double val = 0.0;
  for (int i = 0; i < dim; i++)
    for (int j = i; j < dim; j++)
      for (int k = j; k < dim; k++)
      val += coord[i] * coord[j] * coord[k];
  return val;
}


// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(DriverTest1DGather, 1D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<1>, 0.0);
}

TEST_F(DriverTest1DGather, 1D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<1>, 0.0);
}

TEST_F(DriverTest1DGather, 1D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

TEST_F(DriverTest1DScatter, 1D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_ConstantFieldUnitaryBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_LinearFieldLinearBasisAlt) {
  unitTestAlt<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatterAlt) {
  unitTestAlt<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_LinearFieldLinearBasisAlt) {
  unitTestAlt<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatterAlt) {
  unitTestAlt<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}


TEST_F(IntegrateDriverTest1D, 1D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<1>, 1./2.);
}

TEST_F(IntegrateDriverTest1D, 1D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 1./3.);
}

TEST_F(IntegrateDriverTest2D, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 1.0);
}

TEST_F(IntegrateDriverTest2D, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
    (compute_quadratic_field<2>, 2./3. + 1./4.); // 0.9166666666666666
}

TEST_F(IntegrateDriverTest3D, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
    (compute_linear_field<3>, 3./2.);
}

  TEST(Part, 2D) {
    const size_t NCELLS = 4;
    std::shared_ptr<Wonton::Simple_Mesh> mesh_ptr = 
      std::make_shared<Wonton::Simple_Mesh>(-1., -1., 1., 1., NCELLS, NCELLS);
    Wonton::Simple_Mesh &mesh = *mesh_ptr;
    Wonton::Simple_Mesh_Wrapper mwrapper(mesh);
    double factor = 1.25, bfactor=0.5, dx=0.5;

    size_t ncellsmesh = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
    assert(ncellsmesh == NCELLS*NCELLS);
    std::vector<double> values(ncellsmesh, 1.0);
#define DEBUG_HERE 0
#ifdef DEBUG_HERE 
    std::cout << "src={"<<std::endl;
#endif
    for (int i=0; i<ncellsmesh; i++) {
      Wonton::Point<2> pnt;
      mwrapper.cell_centroid(i, &pnt);
      if      (pnt[0]<0. and pnt[1]<0.) values[i] = 2.;
      else if (pnt[0]>0. and pnt[1]>0.) values[i] = 2.;
#ifdef DEBUG_HERE
      std::cout<<"{"<<pnt[0]<<","<<pnt[1]<<","<<values[i]<<"}";
      if (i<ncellsmesh-1) std::cout <<","<<std::endl;
#endif
    }
    double *valptr = values.data();
    std::cout <<"};"<<std::endl;

    Wonton::Simple_State state(mesh_ptr); 
    std::vector<double> &added = state.add("indicate", Wonton::CELL, valptr);
    Wonton::Simple_State_Wrapper swrapper(state);

    std::shared_ptr<Portage::Meshfree::Swarm<2>> src_swarm_ptr = 
      Portage::Meshfree::SwarmFactory<2, Wonton::Simple_Mesh_Wrapper>(mwrapper, Wonton::CELL);

    std::shared_ptr<Portage::Meshfree::SwarmState<2>> src_state_ptr = 
      Portage::Meshfree::SwarmStateFactory<2, Wonton::Simple_State_Wrapper>(swrapper, Wonton::CELL);
    
    int ntarget = (2*NCELLS+2)*(2*NCELLS+2); // keep all target swarm points from overlying any source points
    auto tgt_swarm_ptr = Portage::Meshfree::SwarmFactory(-1.,-1.,1.,1.,ntarget,1);
    auto tgt_state_ptr = std::make_shared<Portage::Meshfree::SwarmState<2>>(*tgt_swarm_ptr);

    auto tvalues_ptr = std::make_shared<Portage::vector<double>>(ntarget);
    tgt_state_ptr->add_field("indicate",  tvalues_ptr);

    Portage::vector<std::vector<std::vector<double>>> smoothing;
    Portage::vector<Wonton::Point<2>> extents;
    Portage::vector<Wonton::Point<2>> dummy;
    Portage::Meshfree::Weight::faceted_setup_cell<2,Wonton::Simple_Mesh_Wrapper>
                                           (mwrapper, smoothing, extents, 1.0, 2.0);

    Portage::Meshfree::SwarmDriver<Portage::SearchPointsByCells,
                                   Portage::Meshfree::Accumulate,
                                   Portage::Meshfree::Estimate,
                                   2,
                                   Portage::Meshfree::Swarm<2>,
                                   Portage::Meshfree::SwarmState<2>,
                                   Portage::Meshfree::Swarm<2>,
                                   Portage::Meshfree::SwarmState<2>>
      driver(*src_swarm_ptr, *src_state_ptr, *tgt_swarm_ptr, *tgt_state_ptr, smoothing,
               extents, dummy, Portage::Meshfree::Scatter);

    std::vector<std::string> svars={"indicate"};
    std::vector<std::string> tvars={"indicate"};
    Portage::vector<std::vector<std::vector<double>>> psmoothing;
    Portage::Meshfree::Weight::faceted_setup_cell<2,Wonton::Simple_Mesh_Wrapper>
                                           (mwrapper, psmoothing, dummy, 0.25, 1.0);
    driver.set_remap_var_names(svars, tvars, 
                               Portage::Meshfree::LocalRegression, 
                               Portage::Meshfree::Basis::Unitary,
                               Portage::Meshfree::Operator::LastOperator,
                               Portage::vector<Portage::Meshfree::Operator::Domain>(0),
                               Portage::vector<std::vector<Point<2>>>(0,std::vector<Point<2>>(0)),
                               "indicate", 0.25, psmoothing);

    driver.run();

    Portage::Meshfree::SwarmState<2>::DblVecPtr indicator_ptr;
    tgt_state_ptr->copy_field("indicate", indicator_ptr);
    Portage::vector<double> &indicator=*indicator_ptr;

#ifdef DEBUG_HERE
    std::cout << "dat={"<<std::endl;
    for (int i=0; i<ntarget; i++) {
      Wonton::Point<2> p=tgt_swarm_ptr->get_particle_coordinates(i);
      double value=0.;
      if      (p[0]<0. and p[1]<0.) value = 2.;
      else if (p[0]>0. and p[1]>0.) value = 2.;
      else value = 1.0;
      std::cout<<"{"<<p[0]<<","<<p[1]<<","<<indicator[i]<<","<<value<<"}";
      if (i<ntarget-1) std::cout <<","<<std::endl;
    }
    std::cout <<"};"<<std::endl;
#endif
    
    for (size_t i=0; i<ntarget; i++) {
      Wonton::Point<2> pnt=tgt_swarm_ptr->get_particle_coordinates(i);
      if      (pnt[0]<0. and pnt[1]<0.) {ASSERT_NEAR(2.0, indicator[i], 1.e-12);}
      else if (pnt[0]>0. and pnt[1]>0.) {ASSERT_NEAR(2.0, indicator[i], 1.e-12);}
      else                              {ASSERT_NEAR(1.0, indicator[i], 1.e-12);}
    }

  }


  TEST(Part, 3D) {
    const size_t NCELLS = 4;
    std::shared_ptr<Wonton::Simple_Mesh> mesh_ptr = 
      std::make_shared<Wonton::Simple_Mesh>(-1., -1., -1., 1., 1., 1., NCELLS, NCELLS, NCELLS);
    Wonton::Simple_Mesh &mesh = *mesh_ptr;
    Wonton::Simple_Mesh_Wrapper mwrapper(mesh);
    double factor = 1.25, bfactor=0.5, dx=0.5;

    size_t ncellsmesh = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
    assert(ncellsmesh == NCELLS*NCELLS*NCELLS);
    std::vector<double> values(ncellsmesh, 1.0);
#ifdef DEBUG_HERE 
    std::cout << "src={"<<std::endl;
#endif
    for (int i=0; i<ncellsmesh; i++) {
      Wonton::Point<3> pnt;
      mwrapper.cell_centroid(i, &pnt);
      if      (pnt[0]<0. and pnt[1]<0. and pnt[2]<0.) values[i] = 2.;
      else if (pnt[0]>0. and pnt[1]>0. and pnt[2]>0.) values[i] = 2.;
#ifdef DEBUG_HERE
      std::cout<<"{"<<pnt[0]<<","<<pnt[1]<<","<<pnt[2]<<","<<values[i]<<"}";
      if (i<ncellsmesh-1) std::cout <<","<<std::endl;
#endif
    }
    double *valptr = values.data();
    std::cout <<"};"<<std::endl;

    Wonton::Simple_State state(mesh_ptr); 
    std::vector<double> &added = state.add("indicate", Wonton::CELL, valptr);
    Wonton::Simple_State_Wrapper swrapper(state);

    std::shared_ptr<Portage::Meshfree::Swarm<3>> src_swarm_ptr = 
      Portage::Meshfree::SwarmFactory<3, Wonton::Simple_Mesh_Wrapper>(mwrapper, Wonton::CELL);

    std::shared_ptr<Portage::Meshfree::SwarmState<3>> src_state_ptr = 
      Portage::Meshfree::SwarmStateFactory<3, Wonton::Simple_State_Wrapper>(swrapper, Wonton::CELL);
    
    int ntarget = (2*NCELLS+2)*(2*NCELLS+2)*(2*NCELLS+2); // keep all target swarm points from overlying any source points
    auto tgt_swarm_ptr = Portage::Meshfree::SwarmFactory(-1.,-1.,-1.,1.,1.,1.,ntarget,1);
    auto tgt_state_ptr = std::make_shared<Portage::Meshfree::SwarmState<3>>(*tgt_swarm_ptr);

    auto tvalues_ptr = std::make_shared<Portage::vector<double>>(ntarget);
    tgt_state_ptr->add_field("indicate",  tvalues_ptr);

    Portage::vector<std::vector<std::vector<double>>> smoothing;
    Portage::vector<Wonton::Point<3>> extents;
    Portage::vector<Wonton::Point<3>> dummy;
    Portage::Meshfree::Weight::faceted_setup_cell<3,Wonton::Simple_Mesh_Wrapper>
                                           (mwrapper, smoothing, extents, 1.0, 2.0);

    Portage::Meshfree::SwarmDriver<Portage::SearchPointsByCells,
                                   Portage::Meshfree::Accumulate,
                                   Portage::Meshfree::Estimate,
                                   3,
                                   Portage::Meshfree::Swarm<3>,
                                   Portage::Meshfree::SwarmState<3>,
                                   Portage::Meshfree::Swarm<3>,
                                   Portage::Meshfree::SwarmState<3>>
      driver(*src_swarm_ptr, *src_state_ptr, *tgt_swarm_ptr, *tgt_state_ptr, smoothing,
               extents, dummy, Portage::Meshfree::Scatter);

    std::vector<std::string> svars={"indicate"};
    std::vector<std::string> tvars={"indicate"};
    Portage::vector<std::vector<std::vector<double>>> psmoothing;
    Portage::Meshfree::Weight::faceted_setup_cell<3,Wonton::Simple_Mesh_Wrapper>
                                           (mwrapper, psmoothing, dummy, 0.25, 1.0);
    driver.set_remap_var_names(svars, tvars, 
                               Portage::Meshfree::LocalRegression, 
                               Portage::Meshfree::Basis::Unitary,
                               Portage::Meshfree::Operator::LastOperator,
                               std::vector<Portage::Meshfree::Operator::Domain>(0),
                               std::vector<std::vector<Point<3>>>(0,std::vector<Point<3>>(0)),
                               "indicate", 0.25, psmoothing);

    driver.run();

    Portage::Meshfree::SwarmState<3>::DblVecPtr indicator_ptr;
    tgt_state_ptr->copy_field("indicate", indicator_ptr);
    Portage::vector<double> &indicator=*indicator_ptr;

#ifdef DEBUG_HERE
    std::cout << "dat={"<<std::endl;
    for (int i=0; i<ntarget; i++) {
      Wonton::Point<3> p=tgt_swarm_ptr->get_particle_coordinates(i);
      double value=0.;
      if      (p[0]<0. and p[1]<0. and p[2]<0.) value = 2.;
      else if (p[0]>0. and p[1]>0. and p[2]>0.) value = 2.;
      else value = 1.0;
      std::cout<<"{"<<p[0]<<","<<p[1]<<","<<p[2]<<","<<indicator[i]<<","<<value<<"}";
      if (i<ntarget-1) std::cout <<","<<std::endl;
    }
    std::cout <<"};"<<std::endl;
#endif
#undef DEBUG_HERE
    
    for (size_t i=0; i<ntarget; i++) {
      Wonton::Point<3> pnt=tgt_swarm_ptr->get_particle_coordinates(i);
      if      (pnt[0]<0. and pnt[1]<0. and pnt[2]<0.) {ASSERT_NEAR(2.0, indicator[i], 1.e-12);}
      else if (pnt[0]>0. and pnt[1]>0. and pnt[2]>0.) {ASSERT_NEAR(2.0, indicator[i], 1.e-12);}
      else                                            {ASSERT_NEAR(1.0, indicator[i], 1.e-12);}
    }

  }

}  // end namespace
