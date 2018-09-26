/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <cassert>
#include <cmath>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/driver_swarm.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/operator.h"
#include "portage/support/Point.h"
#include "portage/support/portage.h"
#include "portage/search/search_points_by_cells.h"

namespace {

using std::shared_ptr;
using std::make_shared;
using Portage::Meshfree::SwarmFactory;
using Portage::Point;


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

  shared_ptr<Portage::vector<std::vector<std::vector<double>>>> smoothing_lengths_;

  // Operator info
  Portage::Meshfree::Operator::Type operator_;
  Portage::vector<Portage::Meshfree::Operator::Domain> domains_;
  Portage::vector<std::vector<Point<dim>>> operator_data_;

  Portage::Meshfree::WeightCenter center_;

  // Constructor for Driver test
  DriverTest(shared_ptr<Portage::Meshfree::Swarm<dim>> s,
             shared_ptr<Portage::Meshfree::Swarm<dim>> t, 
             Portage::Meshfree::Operator::Type op=Portage::Meshfree::Operator::LastOperator):
    sourceSwarm(s), targetSwarm(t), smoothing_lengths_(nullptr),
    center_(Portage::Meshfree::Gather), operator_(op)
  {
    sourceState = make_shared<Portage::Meshfree::SwarmState<dim>>(*sourceSwarm);
    targetState = make_shared<Portage::Meshfree::SwarmState<dim>>(*targetSwarm); 

    if (op != Portage::Meshfree::Operator::LastOperator) {
      domains_ = Portage::vector<Portage::Meshfree::Operator::Domain>(targetSwarm->num_owned_particles());
      size_t npoints[3]={2,4,8};
      operator_data_ = Portage::vector<std::vector<Point<dim>>>
	(targetSwarm->num_owned_particles(),std::vector<Point<dim>>(npoints[dim-1]));
      size_t npdim = static_cast<size_t>(pow(1.001*operator_data_.size(),1./dim));
      size_t nptot = 1; for (size_t m=0; m<dim; m++) nptot *= npdim;
      assert(nptot == targetSwarm->num_owned_particles());
      assert(npdim>1);
      size_t n=0;
      size_t ij[npoints[dim-1]]; 
      size_t i=0,j=0,k=0;
      Point<dim> pt;
      std::vector<Point<dim>> points(npoints[dim-1]);
      // Create points determining the integration volume.
      // Assumes target swarm is created from SwarmFactory and represents a perfect cubic array of points.
      for (size_t n=0; n<nptot; n++) {
	switch(dim){
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
	}
	for (int m=0; m<npoints[dim-1]; m++) {
	  if (i==npdim-1 or j==npdim-1 or k==npdim-1) {
	    // particles on the upper boundaries get an integration volume of zero
	    pt = targetSwarm->get_particle_coordinates(ij[0]);
	  } else {
	    pt = targetSwarm->get_particle_coordinates(ij[m]);
	  }
	  points[m] = pt;
	}
	operator_data_[n] = points;
	domains_[n] = Portage::Meshfree::Operator::domain_from_points<dim>(points);
      }
    }
  }

  void set_smoothing_lengths(shared_ptr<Portage::vector<std::vector<std::vector<double>>>> smoothing_lengths,
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
    std::vector<std::string> remap_fields;
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

    Portage::Meshfree::EstimateType estimator=Portage::Meshfree::LocalRegression;
    if (operator_ != Portage::Meshfree::Operator::LastOperator) 
      estimator = Portage::Meshfree::OperatorRegression;

    d.set_remap_var_names(remap_fields, remap_fields,
                          estimator, basis, 
                          operator_, domains_, operator_data_);
    // run on one processor
    d.run(false);

    // Check the answer
    double toterr=0.;
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr vecout;
    targetState->get_field("particledata", vecout);
    ASSERT_NE(nullptr, vecout);
    if (operator_ == Portage::Meshfree::Operator::LastOperator) {
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
	{double val=(*vecout)[p]; 
	  std::printf("  Value = % 10.6lf  Err = % lf\n", val, error);}
        toterr += error*error;
      }
    
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
      ASSERT_NEAR(expected_answer, sqrt(toterr), TOL);
    } else if (operator_ == Portage::Meshfree::Operator::VolumeIntegral) {
      double total = 0.;
      for (int p = 0; p < ntarpts; ++p) {
        total += (*vecout)[p];
      }
      ASSERT_NEAR(expected_answer, total, TOL);
    }
  }

};



// Class which constructs a pair of 1-D swarms (random distribution) for remaps
struct DriverTest1D : DriverTest<1> {
  DriverTest1D() : DriverTest(SwarmFactory(0.0, 1.0, 7, 2),
                              SwarmFactory(0.0, 1.0, 5, 2))
  {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(5,
                   std::vector<std::vector<double>>(1, std::vector<double>(1, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};


// Class which constructs a pair of 2-D swarms (random distribution) for remaps
struct DriverTest2D : DriverTest<2> {
  DriverTest2D() : DriverTest(SwarmFactory(0.0, 0.0, 1.0, 1.0, 7*7, 2),
                              SwarmFactory(0.0, 0.0, 1.0, 1.0, 5*5, 2)) {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(5*5,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};


// Class which constructs a pair of 3-D swarms (random distribution) for remaps
struct DriverTest3D : DriverTest<3> {
  DriverTest3D(): DriverTest(SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          7*7*7, 2),
                             SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          5*5*5, 2)) {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(5*5*5,
                   std::vector<std::vector<double>>(1, std::vector<double>(3, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};



// Class which constructs a pair of 1-D swarms (random distribution) for remaps
struct DriverTest1DScatter : DriverTest<1> {
  DriverTest1DScatter() : DriverTest(SwarmFactory(0.0, 1.0, 7, 2),
                              SwarmFactory(0.0, 1.0, 5, 2))
  {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(7,
                   std::vector<std::vector<double>>(1, std::vector<double>(1, 2.0/6)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};


// Class which constructs a pair of 2-D swarms (random distribution) for remaps
struct DriverTest2DScatter : DriverTest<2> {
  DriverTest2DScatter() : DriverTest(SwarmFactory(0.0, 0.0, 1.0, 1.0, 7*7, 2),
                              SwarmFactory(0.0, 0.0, 1.0, 1.0, 5*5, 2)) {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(7*7,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 2.0/6)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};


// Class which constructs a pair of 3-D swarms (random distribution) for remaps
struct DriverTest3DScatter : DriverTest<3> {
  DriverTest3DScatter(): DriverTest(SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          7*7*7, 2),
                             SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          5*5*5, 2)) {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(7*7*7,
                   std::vector<std::vector<double>>(1, std::vector<double>(3, 2.0/6)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};



// Class which constructs a pair of 1-D swarms (ordered distribution) for remaps, and integrates
struct IntegrateDriverTest1D : DriverTest<1> {
  IntegrateDriverTest1D() : DriverTest(SwarmFactory(0.0, 1.0, 7, 1),
                                       SwarmFactory(0.0, 1.0, 5, 1),
                                       Portage::Meshfree::Operator::VolumeIntegral)
  {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(5,
                   std::vector<std::vector<double>>(1, std::vector<double>(1, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};


// Class which constructs a pair of 2-D swarms (ordered distribution) for remaps, and integrates
struct IntegrateDriverTest2D : DriverTest<2> {
  IntegrateDriverTest2D() : DriverTest(SwarmFactory(0.0, 0.0, 1.0, 1.0, 7*7, 1),
                                       SwarmFactory(0.0, 0.0, 1.0, 1.0, 5*5, 1),
                                       Portage::Meshfree::Operator::VolumeIntegral) {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(5*5,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};


// Class which constructs a pair of 3-D swarms (ordered distribution) for remaps, and integrates
struct IntegrateDriverTest3D : DriverTest<3> {
  IntegrateDriverTest3D(): DriverTest(SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          7*7*7, 1),
                                      SwarmFactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          5*5*5, 1),
                                      Portage::Meshfree::Operator::VolumeIntegral) {
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(5*5*5,
                   std::vector<std::vector<double>>(1, std::vector<double>(3, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
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
  for (size_t i = 0; i < Dimension; i++) 
    for (size_t j = i; j < Dimension; j++) 
      val += coord[i]*coord[j];
  return val;
}

template<size_t Dimension>
double compute_cubic_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++)
    for (size_t j = i; j < Dimension; j++) 
      for (size_t k = j; k < Dimension; k++) 
      val += coord[i]*coord[j]*coord[k];
  return val;
}


// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(DriverTest1D, 1D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<1>, 0.0);
}

TEST_F(DriverTest1D, 1D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<1>, 0.0);
}

TEST_F(DriverTest1D, 1D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

TEST_F(DriverTest1DScatter, 1D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

TEST_F(DriverTest2D, 2D_ConstantFieldUnitaryBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2D, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2D, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest3D, 3D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3D, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3D, 3D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
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

}  // end namespace
