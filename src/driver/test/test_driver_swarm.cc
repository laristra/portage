/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/driver_swarm.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

namespace {

using std::vector;
using std::shared_ptr;
using std::make_shared;
using Portage::Meshfree::SwarmFactory;


double TOL = 1e-6;

// This is a set of integration tests for the swarm-swarm remap driver.
// There will be at least one test corresponding to each case found in
// main.cc. This is a test fixture and must be derived from the
// ::testing::Test class. Specializations of this class, such as 2D/3D
// coincident and non-coincident remaps should be derived from this.

template<size_t dim>
class DriverTest : public ::testing::Test {
 protected:
  // Source and target swarms
  shared_ptr<Portage::Meshfree::Swarm<dim>> sourceSwarm;
  shared_ptr<Portage::Meshfree::Swarm<dim>> targetSwarm;

  // Source and target mesh state
  shared_ptr<Portage::Meshfree::SwarmState<dim>> sourceState;
  shared_ptr<Portage::Meshfree::SwarmState<dim>> targetState;

  shared_ptr<vector<vector<vector<double>>>> smoothing_lengths_;

  Portage::Meshfree::WeightCenter center_;

  // Constructor for Driver test
  DriverTest(shared_ptr<Portage::Meshfree::Swarm<dim>> s,
             shared_ptr<Portage::Meshfree::Swarm<dim>> t) :
    sourceSwarm(s), targetSwarm(t), smoothing_lengths_(nullptr),
    center_(Portage::Meshfree::Gather) {
    sourceState = make_shared<Portage::Meshfree::SwarmState<dim>>(*sourceSwarm);
    targetState = make_shared<Portage::Meshfree::SwarmState<dim>>(*targetSwarm);
  }

  void set_smoothing_lengths(shared_ptr<vector<vector<vector<double>>>>
                             smoothing_lengths,
			     Portage::Meshfree::WeightCenter center=Portage::Meshfree::Gather) {
    smoothing_lengths_ = smoothing_lengths;
    center_ = center;
  }



  // This is the basic test method to be called for each unit test.
  // It will work for 1, 2-D and 3-D swarms
  //
  template <template<int, class, class> class Search,
            Portage::Meshfree::Basis::Type basis>
  void unitTest(double compute_initial_field(Portage::Point<dim> coord),
                double expected_answer) {

    // Fill the source state data with the specified profile
    const int nsrcpts = sourceSwarm->num_owned_particles();
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr sourceData = 
        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(nsrcpts);

    // Create the source data for given function
    for (unsigned int p = 0; p < nsrcpts; ++p) {
      Portage::Point<dim> coord =
          sourceSwarm->get_particle_coordinates(p);
      (*sourceData)[p] = compute_initial_field(coord);
    }
    sourceState->add_field("particledata", sourceData);

    // Build the target state storage
    const int ntarpts = targetSwarm->num_owned_particles();
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr targetData = 
        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(ntarpts, 0.0);
    targetState->add_field("particledata", targetData);

    // Build the main driver data for this mesh type
    // Register the variable name and interpolation order with the driver
    vector<std::string> remap_fields;
    remap_fields.push_back("particledata");

    Portage::Meshfree::SwarmDriver<Search,
                                   Portage::Meshfree::Accumulate,
                                   Portage::Meshfree::Estimate,
                                   dim,
                                   Portage::Meshfree::Swarm<dim>,
                                   Portage::Meshfree::SwarmState<dim>,
                                   Portage::Meshfree::Swarm<dim>,
                                   Portage::Meshfree::SwarmState<dim>>
        d(*sourceSwarm, *sourceState, *targetSwarm, *targetState, *smoothing_lengths_,
	  Portage::Meshfree::Weight::B4, Portage::Meshfree::Weight::ELLIPTIC, center_);
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          basis);
    // run on one processor
    d.run(false);

    // Check the answer
    Portage::Point<dim> nodexy;
    double stdval, err;
    double toterr=0.;

    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr vecout;
    targetState->get_field("particledata", vecout);
    ASSERT_NE(nullptr, vecout);

    for (int p = 0; p < ntarpts; ++p) {
      Portage::Point<dim> coord = targetSwarm->get_particle_coordinates(p);
      double error;
      error = compute_initial_field(coord) - (*vecout)[p];
      // dump diagnostics for each particle
      if (dim == 1)
        std::printf("Particle=% 4d Coord = (% 5.3lf)", p, coord[0]);
      else if (dim == 2)
        std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf)", p, coord[0],
                    coord[1]);
      else if (dim == 3)
        std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf,% 5.3lf)", p,
                    coord[0], coord[1], coord[2]);
      std::printf("  Value = % 10.6lf  Err = % lf\n", (*vecout)[p], error);
      toterr += error*error;
    }
    
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    ASSERT_NEAR(expected_answer, sqrt(toterr), TOL);
  }

};



// Class which constructs a pair of 1-D swarms (random distribution) for remaps
struct DriverTest1D : DriverTest<1> {
  DriverTest1D() : DriverTest(SwarmFactory(0.0, 1.0, 7, 2),
                              SwarmFactory(0.0, 1.0, 5, 2))
  {
    auto smoothing_lengths = make_shared<vector<vector<vector<double>>>>(5,
                   vector<vector<double>>(1, vector<double>(1, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};


// Class which constructs a pair of 2-D swarms (random distribution) for remaps
struct DriverTest2D : DriverTest<2> {
  DriverTest2D() : DriverTest(SwarmFactory(0.0, 0.0, 1.0, 1.0, 7*7, 2),
                              SwarmFactory(0.0, 0.0, 1.0, 1.0, 5*5, 2)) {
    auto smoothing_lengths = make_shared<vector<vector<vector<double>>>>(5*5,
                   vector<vector<double>>(1, vector<double>(2, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};


// Class which constructs a pair of 3-D swarms (random distribution) for remaps
struct DriverTest3D : DriverTest<3> {
  DriverTest3D(): DriverTest(SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          7*7*7, 2),
                             SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          5*5*5, 2)) {
    auto smoothing_lengths = make_shared<vector<vector<vector<double>>>>(5*5*5,
                   vector<vector<double>>(1, vector<double>(3, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};



// Class which constructs a pair of 1-D swarms (random distribution) for remaps
struct DriverTest1DScatter : DriverTest<1> {
  DriverTest1DScatter() : DriverTest(SwarmFactory(0.0, 1.0, 7, 2),
                              SwarmFactory(0.0, 1.0, 5, 2))
  {
    auto smoothing_lengths = make_shared<vector<vector<vector<double>>>>(7,
                   vector<vector<double>>(1, vector<double>(1, 2.0/6)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};


// Class which constructs a pair of 2-D swarms (random distribution) for remaps
struct DriverTest2DScatter : DriverTest<2> {
  DriverTest2DScatter() : DriverTest(SwarmFactory(0.0, 0.0, 1.0, 1.0, 7*7, 2),
                              SwarmFactory(0.0, 0.0, 1.0, 1.0, 5*5, 2)) {
    auto smoothing_lengths = make_shared<vector<vector<vector<double>>>>(7*7,
                   vector<vector<double>>(1, vector<double>(2, 2.0/6)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};


// Class which constructs a pair of 3-D swarms (random distribution) for remaps
struct DriverTest3DScatter : DriverTest<3> {
  DriverTest3DScatter(): DriverTest(SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          7*7*7, 2),
                             SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          5*5*5, 2)) {
    auto smoothing_lengths = make_shared<vector<vector<vector<double>>>>(7*7*7,
                   vector<vector<double>>(1, vector<double>(3, 2.0/6)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};


template<size_t Dimension>
    double compute_constant_field(Portage::Point<Dimension> coord) {
  return 25.0;
}

// Methods for computing initial field values
template<size_t Dimension>
    double compute_linear_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++) val += coord[i];
  return val;
}

template<size_t Dimension>
    double compute_quadratic_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++) val += coord[i]*coord[i];
  return val;
}

template<size_t Dimension>
    double compute_cubic_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++) val += coord[i]*coord[i]*coord[i];
  return val;
}


// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(DriverTest1D, 1D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<1>, 0.0);
}

//TEST_F(DriverTest1D, 1D_LinearFieldUnitaryBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Unitary>
//      (compute_linear_field<1>, 0.0095429796560267122);
//}

TEST_F(DriverTest1D, 1D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<1>, 0.0);
}

//TEST_F(DriverTest1D, 1D_QuadraticFieldLinearBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Linear>
//      (compute_quadratic_field<1>, 0.0095429796560267122);
//}

TEST_F(DriverTest1D, 1D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

  /*
TEST_F(DriverTest1DScatter, 1D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}
  */

//TEST_F(DriverTest1D, 1D_CubicFieldQuadraticBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
//      (compute_quadratic_field<1>, 0.0);
//}



TEST_F(DriverTest2D, 2D_ConstantFieldUnitaryBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

//TEST_F(DriverTest2D, 2D_LinearFieldUnitaryBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Unitary>
//      (compute_linear_field<2>, 0.0095429796560267122);
//}

TEST_F(DriverTest2D, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

//TEST_F(DriverTest2D, 2D_QuadraticFieldLinearBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Linear>
//      (compute_quadratic_field<2>, 0.0095429796560267122);
//}

TEST_F(DriverTest2D, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

  /*
TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}
  */

//TEST_F(DriverTest2D, 2D_CubicFieldQuadraticBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
//      (compute_quadratic_field<2>, 0.0);
//}


TEST_F(DriverTest3D, 3D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<3>, 0.0);
}

//TEST_F(DriverTest3D, 3D_LinearFieldUnitaryBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Unitary>
//      (compute_linear_field<3>, 0.0095429796560267122);
//}

TEST_F(DriverTest3D, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

//TEST_F(DriverTest3D, 3D_QuadraticFieldLinearBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Linear>
//      (compute_quadratic_field<3>, 0.0095429796560267122);
//}

TEST_F(DriverTest3D, 3D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

  /*
TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}
  */

//TEST_F(DriverTest1D, 3D_CubicFieldQuadraticBasis) {
//  unitTest<Portage::SearchSimplePoints, Portage::Meshfree::Basis::Quadratic>
//      (compute_quadratic_field<3>, 0.0);
//}

}  // end namespace
