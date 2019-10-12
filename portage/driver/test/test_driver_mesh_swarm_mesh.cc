/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <string>
#include <limits>

#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

// portage includes
#include "portage/driver/mmdriver.h"
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/operator.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_mm_wrapper.h"
#include "wonton/state/state_vector_uni.h"
#include "wonton/mesh/flat/flat_mesh_wrapper.h"

namespace {

double TOL = 1e-6;

// This is a set of integration tests based off of main.cc.  There
// will be at least one test corresponding to each case found in
// main.cc.  This is a test fixture and must be derived from the
// ::testing::Test class.  Specializations of this class, such as
// 2D/3D coincident and non-coincident remaps should be derived from
// this.
class MSMDriverTest : public ::testing::Test {
 protected:
  // Source and target meshes
  std::shared_ptr<Wonton::Simple_Mesh> sourceMesh;
  std::shared_ptr<Wonton::Simple_Mesh> targetMesh;
  
  // Wrappers for interfacing with the underlying mesh data structures
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper;
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> sourceStateWrapper;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> targetStateWrapper;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> targetStateWrapper2;
  
  // Operator domains and data
  Portage::vector<Portage::Meshfree::Operator::Domain> domains_;

  //  This is the basic test method to be called for each unit test. It will work
  //  for 2-D and 3-D, coincident and non-coincident cell-centered remaps.
  template <
    template<Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class,
    class, class> class Intersect,
    template<int, Portage::Entity_kind, class, class, class, class,
    template<class, int, class, class> class,
    class, class, class> class Interpolate,
    template <int, class, class> class SwarmSearch,
    int Dimension = 3
  >
  void unitTest(double compute_initial_field(Portage::Point<Dimension> centroid),
                double smoothing_factor, Portage::Meshfree::Basis::Type basis, 
		Portage::Meshfree::WeightCenter center=Portage::Meshfree::Gather, 
		Portage::Meshfree::Operator::Type oper8or=Portage::Meshfree::Operator::LastOperator,
		bool faceted=false)
  {
    // check dimension - no 1D
    assert(Dimension > 1);

    //  Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells();
    std::vector<double> sourceData(nsrccells);
    const int nsrcnodes = sourceMeshWrapper.num_owned_nodes();
    std::vector<double> sourceDataNode(nsrcnodes);

    // Create the source data for given function
    Wonton::Flat_Mesh_Wrapper<double> sourceFlatMesh;
    sourceFlatMesh.initialize(sourceMeshWrapper);
  
    for (unsigned int c = 0; c < nsrccells; ++c) {
      Portage::Point<Dimension> cen;
      sourceFlatMesh.cell_centroid(c, &cen);
      sourceData[c] = compute_initial_field(cen);
    }
    
    sourceStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"celldata", Portage::Entity_kind::CELL, sourceData)
    );

    for (unsigned int c = 0; c < nsrcnodes; ++c) {
      Portage::Point<Dimension> cen;
      sourceFlatMesh.node_get_coordinates(c, &cen);
      sourceDataNode[c] = compute_initial_field(cen);
    }
    sourceStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"nodedata", Portage::Entity_kind::NODE, sourceDataNode));

    // Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    std::vector<double> targetData(ntarcells), targetData2(ntarcells);
    std::vector<double> targetDataNode(ntarnodes), targetData2Node(ntarnodes);
    
    targetStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"celldata", Portage::Entity_kind::CELL, targetData));
    targetStateWrapper2.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"celldata", Portage::Entity_kind::CELL, targetData2));
    targetStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"nodedata", Portage::Entity_kind::NODE, targetDataNode));
    targetStateWrapper2.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"nodedata", Portage::Entity_kind::NODE, targetData2Node));

    // Make list of field names
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    if (not faceted) {
      remap_fields.push_back("nodedata");
    }

    //  Build the mesh-mesh driver data for this mesh type
    Portage::MMDriver<Portage::SearchKDTree,
    Intersect,
    Interpolate,
    Dimension,
    Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>,
    Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>>
    mmdriver(sourceMeshWrapper, sourceStateWrapper,
             targetMeshWrapper, targetStateWrapper);
    mmdriver.set_remap_var_names(remap_fields);
    mmdriver.run();  // run in serial (executor argument defaults to nullptr)

    // Set up the operator information if needed
    Portage::vector<std::vector<Portage::Point<Dimension>>> data;
    std::vector<double> exact;
    Portage::Meshfree::Operator::Domain domain_types[3]=
      {Portage::Meshfree::Operator::Interval,
       Portage::Meshfree::Operator::Quadrilateral,
       Portage::Meshfree::Operator::Hexahedron};
    if (oper8or == Portage::Meshfree::Operator::VolumeIntegral) {
      int numcells = targetMesh->num_entities(Portage::Entity_kind::CELL, 
                                              Portage::Entity_type::ALL);
      domains_.resize(numcells);
      data.resize(numcells);
      exact.resize(numcells);  // each element is resized automatically
      for (int c=0; c<numcells; c++) {
	// get integration domains
        domains_[c] = domain_types[Dimension-1];
        std::vector<Portage::Point<Dimension>> cellverts;
        targetMesh->cell_get_coordinates(c, &cellverts);
        data[c] = cellverts;

	// get exact value of integral
	std::vector<std::vector<double>> result;
	Portage::Meshfree::Operator::apply<Dimension>
          ( Portage::Meshfree::Operator::VolumeIntegral, 
            basis, domain_types[Dimension-1], data[c], result);
	if (Dimension==2) {
	  if (basis==Portage::Meshfree::Basis::Linear) 
	    exact[c] = result[1][0]+result[2][0];
	  if (basis==Portage::Meshfree::Basis::Quadratic) 
	    exact[c] = 2.*result[3][0]+result[4][0]+2.*result[5][0];
	}
	if (Dimension==3) {
	  if (basis==Portage::Meshfree::Basis::Linear) 
	    exact[c] = result[1][0]+result[2][0]+result[3][0];
	  if (basis==Portage::Meshfree::Basis::Quadratic) 
	    exact[c] = 2.*result[4][0]+result[5][0]+result[6][0]+
	               2.*result[7][0]+result[8][0]+2.*result[9][0];
	}
      }
    }

    //  Build the mesh-swarm-mesh driver data for this mesh type
    using MSM_Driver_Type =
    Portage::MSM_Driver<
      SwarmSearch,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      Dimension,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>
      >;

    MSM_Driver_Type *msmdriver_ptr;

    if (faceted) {
      msmdriver_ptr = new MSM_Driver_Type
        (sourceMeshWrapper, sourceStateWrapper,
         targetMeshWrapper, targetStateWrapper2,
         smoothing_factor, 
         smoothing_factor, 
         Portage::Meshfree::Weight::FACETED, 
         Portage::Meshfree::Weight::POLYRAMP, 
         center, 
         std::string("NONE"), 
         std::numeric_limits<double>::infinity());
    } else {
      msmdriver_ptr = new MSM_Driver_Type
        (sourceMeshWrapper, sourceStateWrapper,
         targetMeshWrapper, targetStateWrapper2,
         smoothing_factor, 
         smoothing_factor, 
         Portage::Meshfree::Weight::TENSOR, 
         Portage::Meshfree::Weight::B4, 
         center);
    }
    MSM_Driver_Type &msmdriver(*msmdriver_ptr);
    
    Portage::Meshfree::EstimateType estimator=
      Portage::Meshfree::LocalRegression;
    if (oper8or == Portage::Meshfree::Operator::VolumeIntegral) 
      estimator = Portage::Meshfree::OperatorRegression;

    msmdriver.set_remap_var_names(remap_fields, remap_fields, 
                                  estimator, 
                                  basis,
                                  oper8or,
                                  domains_,
                                  data);

    msmdriver.run();  // run in serial (executor argument defaults to nullptr)

    // Check the answer
    double stdval, err;
    double toterr = 0.;

    auto& cellvecout = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper.get("celldata"))->get_data();
    auto& cellvecout2 = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper2.get("celldata"))->get_data();
    auto& nodevecout = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper.get("nodedata"))->get_data();
    auto& nodevecout2 = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper2.get("nodedata"))->get_data();

    Wonton::Flat_Mesh_Wrapper<double> targetFlatMesh;
    targetFlatMesh.initialize(targetMeshWrapper);
  
    if (oper8or != Portage::Meshfree::Operator::VolumeIntegral) {

      for (int c = 0; c < ntarcells; ++c) {
	Portage::Point<Dimension> ccen;
	targetFlatMesh.cell_centroid(c, &ccen);
	double value = compute_initial_field(ccen);
	double merror, serror;
	merror = cellvecout[c] - value;
	serror = cellvecout2[c] - value;
	//  dump diagnostics for each cell
	if (Dimension == 2) {
	  std::printf("Cell=% 4d Centroid=(% 5.3lf,% 5.3lf)", c,
		      ccen[0], ccen[1]);
	} else if (Dimension == 3) {
	  std::printf("Cell=% 4d Centroid=(% 5.3lf,% 5.3lf,% 5.3lf)", c,
		      ccen[0], ccen[1], ccen[2]);
	}
	std::printf(" Val=% 10.6lf MM=% 10.6lf Err=% 10.6lf MSM=% 10.6lf Err=% lf\n",
		    value, cellvecout[c], merror, cellvecout2[c], serror);
	toterr = std::max(toterr, std::fabs(serror));
      }

      std::printf("\n\nLinf NORM OF MSM CELL ERROR = %lf\n\n", toterr);
      ASSERT_LT(toterr, TOL);

      if (not faceted) {
        toterr = 0.;
        for (int n = 0; n < ntarnodes; ++n) {
          Portage::Point<Dimension> node;
          targetFlatMesh.node_get_coordinates(n, &node);
          double value = compute_initial_field(node);
          double merror, serror;
          merror = nodevecout[n] - value;
          serror = nodevecout2[n] - value;
          //  dump diagnostics for each node
          if (Dimension == 2) {
            std::printf("Node=% 4d Coords=(% 5.3lf,% 5.3lf)", n,
                        node[0], node[1]);
          } else if (Dimension == 3) {
            std::printf("Node=% 4d Coords=(% 5.3lf,% 5.3lf,% 5.3lf)", n,
                        node[0], node[1], node[2]);
          }
          std::printf(" Val=% 10.6lf MM=% 10.6lf Err=% 10.6lf MSM=% 10.6lf Err=% lf\n",
                      value, nodevecout[n], merror, nodevecout2[n], serror);
          toterr = std::max(toterr, std::fabs(serror));
        }

        std::printf("\n\nLinf NORM OF MSM NODE ERROR = %lf\n\n", toterr);
        ASSERT_LT(toterr, TOL);
      }

    } else {

      for (int c = 0; c < ntarcells; ++c) {
	Portage::Point<Dimension> ccen;
	targetFlatMesh.cell_centroid(c, &ccen);
	double serror;
	serror = cellvecout2[c] - exact[c];
	if (Dimension == 2) {
	  std::printf("Centroid=% 4d Coords=(% 5.3lf,% 5.3lf)", c,
		      ccen[0], ccen[1]);
	} else if (Dimension == 3) {
	  std::printf("Centroid=% 4d Coords=(% 5.3lf,% 5.3lf,% 5.3lf)", c,
		      ccen[0], ccen[1], ccen[2]);
	}
	std::printf(" Exact=% 10.6lf MSM=% 10.6lf Err=% 10.6lf\n",
		    exact[c], cellvecout2[c], serror);
	toterr = std::max(toterr, std::fabs(serror));
      }

      std::printf("\n\nLinf NORM OF MSM OPERATOR ERROR = %lf\n\n", toterr);
      ASSERT_LT(toterr, TOL);

    }

    delete msmdriver_ptr;
  }

  // Constructor for Driver test
  MSMDriverTest(std::shared_ptr<Wonton::Simple_Mesh> s,
                std::shared_ptr<Wonton::Simple_Mesh> t) :
    sourceMesh(s), targetMesh(t),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceMeshWrapper), targetStateWrapper(targetMeshWrapper),
    targetStateWrapper2(targetMeshWrapper)
  {}

};

// Class which constructs a pair of simple 2-D meshes, target
// contained in source
struct MSMDriverTest2D : MSMDriverTest {
  MSMDriverTest2D(): MSMDriverTest(
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 10, 10),
    std::make_shared<Wonton::Simple_Mesh>(0.3, 0.3, 0.7, 0.7, 4,  4)) {}
};


// Class which constructs a pair of simple 3-D meshes, target
// contained in source
struct MSMDriverTest3D : MSMDriverTest {
  MSMDriverTest3D(): MSMDriverTest(
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10),
    std::make_shared<Wonton::Simple_Mesh>(0.3, 0.3, 0.3, 0.7, 0.7, 0.7,  4,  4,  4)) {}
};

// Methods for computing initial field values

double compute_linear_field_2d(Portage::Point<2> centroid) {
  return centroid[0]+centroid[1];
}
double compute_quadratic_field_2d(Portage::Point<2> centroid) {
  return centroid[0]*centroid[0] + centroid[1]*centroid[1] +
      centroid[0]*centroid[1];
}

// Methods for computing initial field values
double compute_linear_field_3d(Portage::Point<3> centroid) {
  return centroid[0]+centroid[1]+centroid[2];
}
double compute_quadratic_field_3d(Portage::Point<3> centroid) {
  return centroid[0]*centroid[0] + centroid[1]*centroid[1] +
      centroid[2]*centroid[2] + centroid[0]*centroid[1] +
      centroid[1]*centroid[2] + centroid[2]*centroid[0];
}

// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(MSMDriverTest2D, 2D1stOrderLinear) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_linear_field_2d, .75, Portage::Meshfree::Basis::Linear);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadratic) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, Portage::Meshfree::Basis::Quadratic);
}

TEST_F(MSMDriverTest2D, 2D1stOrderLinearScatter) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_linear_field_2d, .75, Portage::Meshfree::Basis::Linear, 
     Portage::Meshfree::Scatter);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticScatter) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Scatter);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticGatherFaceted) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Gather, Portage::Meshfree::Operator::LastOperator, true);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticScatterFaceted) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Scatter, Portage::Meshfree::Operator::LastOperator, true);
}

TEST_F(MSMDriverTest2D, 2D1stOrderLinearIntegrate) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_linear_field_2d, .75, Portage::Meshfree::Basis::Linear, 
     Portage::Meshfree::Gather, Portage::Meshfree::Operator::VolumeIntegral);
}


TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticIntegrate) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Gather, Portage::Meshfree::Operator::VolumeIntegral);
}


TEST_F(MSMDriverTest3D, 3D1stOrderLinear) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, Portage::Meshfree::Basis::Linear);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadratic) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, Portage::Meshfree::Basis::Quadratic);
}

TEST_F(MSMDriverTest3D, 3D1stOrderLinearScatter) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, Portage::Meshfree::Basis::Linear, 
     Portage::Meshfree::Scatter);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticScatter) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Scatter);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticGatherFaceted) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Gather, Portage::Meshfree::Operator::LastOperator, true);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticScatterFaceted) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Scatter, Portage::Meshfree::Operator::LastOperator, true);
}

TEST_F(MSMDriverTest3D, 3D1stOrderLinearIntegrate) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, Portage::Meshfree::Basis::Linear, 
     Portage::Meshfree::Gather, Portage::Meshfree::Operator::VolumeIntegral);
}


}  // namespace
