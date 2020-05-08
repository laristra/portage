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
#ifdef WONTON_ENABLE_MPI
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
// avoid long namespaces
using namespace Portage::Meshfree;
// default numerical tolerance
double const epsilon = 1e-6;

// This is a set of integration tests based off of main.cc.  There
// will be at least one test corresponding to each case found in
// main.cc.  This is a test fixture and must be derived from the
// ::testing::Test class.  Specializations of this class, such as
// 2D/3D coincident and non-coincident remaps should be derived from
// this.
class MSMDriverTest : public ::testing::Test {
public:
  /**
   * @brief Driver test constructor
   *
   * @param source: a shared pointer to the source mesh
   * @param target: a shared pointer to the target mesh
   */
  MSMDriverTest(std::shared_ptr<Wonton::Simple_Mesh> source,
                std::shared_ptr<Wonton::Simple_Mesh> target)
    : source_mesh(source), target_mesh(target),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state(source_mesh_wrapper),
      target_state_one(target_mesh_wrapper),
      target_state_two(target_mesh_wrapper) {}

  /**
   * @brief Test main method.
   *
   * @tparam Intersect
   * @tparam Interpolate
   * @tparam SwarmSearch
   * @tparam dim
   * @param compute_initial_field
   * @param smoothing_factor
   * @param basis
   * @param center
   * @param op
   * @param faceted
   */
  template <
    template<Wonton::Entity_kind,
      class, class, class,
      template<class, int, class, class> class,
      class, class> class Intersect,
    template<int, Wonton::Entity_kind,
      class, class, class, class, class,
      template<class, int, class, class> class,
      class, class, class> class Interpolate,
    template <int, class, class> class SwarmSearch,
    int dim = 3
  >
  void unitTest(double compute_initial_field(Wonton::Point<dim> const& centroid),
                double smoothing_factor,
                basis::Type basis,
                WeightCenter center = Gather,
                oper::Type op = oper::LastOperator,
                bool faceted = false) {
    // check dimension - no 1D
    assert(dim > 1);

    using Field = Wonton::StateVectorUni<>;

    using MeshRemap = Portage::MMDriver<Portage::SearchKDTree,
                                       Intersect, Interpolate, dim,
                                       Wonton::Simple_Mesh_Wrapper,
                                       Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>>;

    using SwarmRemap = Portage::MSM_Driver<SwarmSearch, Accumulate, Estimate, dim,
                                           Wonton::Simple_Mesh_Wrapper,
                                           Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>>;

    int const nb_source_cells = source_mesh_wrapper.num_owned_cells();
    int const nb_source_nodes = source_mesh_wrapper.num_owned_nodes();
    int const nb_target_cells = target_mesh_wrapper.num_owned_cells();
    int const nb_target_nodes = target_mesh_wrapper.num_owned_nodes();

    //  Fill the source state data with the specified profile
    std::vector<double> cell_data(nb_source_cells);
    std::vector<double> node_data(nb_source_nodes);

    // Create the source data for given function
    Wonton::Flat_Mesh_Wrapper<double> source_flat_mesh;
    source_flat_mesh.initialize(source_mesh_wrapper);
  
    for (int i = 0; i < nb_source_cells; ++i) {
      Wonton::Point<dim> cen;
      source_flat_mesh.cell_centroid(i, &cen);
      cell_data[i] = compute_initial_field(cen);
    }

    for (int i = 0; i < nb_source_nodes; ++i) {
      Wonton::Point<dim> cen;
      source_flat_mesh.node_get_coordinates(i, &cen);
      node_data[i] = compute_initial_field(cen);
    }

    source_state.add(std::make_shared<Field>("celldata", Wonton::CELL, cell_data));
    source_state.add(std::make_shared<Field>("nodedata", Wonton::NODE, node_data));

    // Build the target state storage
    std::vector<double> target_data[4];
    target_data[0].resize(nb_target_cells);
    target_data[1].resize(nb_target_cells);
    target_data[2].resize(nb_target_nodes);
    target_data[3].resize(nb_target_nodes);

    target_state_one.add(std::make_shared<Field>("celldata", Wonton::CELL, target_data[0]));
    target_state_two.add(std::make_shared<Field>("celldata", Wonton::CELL, target_data[1]));
    target_state_one.add(std::make_shared<Field>("nodedata", Wonton::NODE, target_data[2]));
    target_state_two.add(std::make_shared<Field>("nodedata", Wonton::NODE, target_data[3]));

    // Make list of field names
    std::vector<std::string> remap_fields = { "celldata" };
    if (not faceted)
      remap_fields.emplace_back("nodedata");

    //  Build the mesh-mesh driver data for this mesh type
    MeshRemap mesh_remap(source_mesh_wrapper, source_state,
                         target_mesh_wrapper, target_state_one);
    mesh_remap.set_remap_var_names(remap_fields);
    mesh_remap.run();  // run in serial (executor argument defaults to nullptr)

    // Set up the operator information if needed
    Wonton::vector<std::vector<Wonton::Point<dim>>> data;
    std::vector<double> exact;
    oper::Domain domain_types[3] = {oper::Interval,
                                    oper::Quadrilateral,
                                    oper::Hexahedron };

    if (op == oper::VolumeIntegral) {
      int ncells = target_mesh->num_entities(Wonton::CELL,Wonton::ALL);
      domains_.resize(ncells);
      data.resize(ncells);
      exact.resize(ncells);  // each element is resized automatically

      for (int i = 0; i < ncells; i++) {

        std::vector<Wonton::Point<dim>> vertices;
        std::vector<std::vector<double>> result;

	      // get integration domains
        domains_[i] = domain_types[dim - 1];
        target_mesh->cell_get_coordinates(i, &vertices);
        data[i] = vertices;

        // get exact value of integral
        oper::apply<dim>(oper::VolumeIntegral, basis,
                         domain_types[dim-1], data[i], result);
        if (dim == 2) {
          switch (basis) {
            case basis::Linear:
              exact[i] = result[1][0] + result[2][0]; break;
            case basis::Quadratic:
              exact[i] = 2. * result[3][0] + result[4][0] + 2. * result[5][0]; break;
            default:
              break;
          }
        } else if (dim == 3) {
          switch (basis) {
            case basis::Linear:
              exact[i] = result[1][0] + result[2][0] + result[3][0]; break;
            case basis::Quadratic:
              exact[i] = 2. * result[4][0] + result[5][0] + result[6][0]
                       + 2. * result[7][0] + result[8][0] + 2. * result[9][0]; break;
            default:
              break;
          }
        }
      }
    }

    //  Build the mesh-swarm-mesh driver data for this mesh type
    SwarmRemap swarm_remap(source_mesh_wrapper, source_state,
                           target_mesh_wrapper, target_state_two,
                           smoothing_factor, smoothing_factor,
                           faceted ? Weight::FACETED : Weight::TENSOR,
                           faceted ? Weight::POLYRAMP : Weight::B4, center);

    auto estimator = (op == oper::VolumeIntegral ? OperatorRegression
                                                 : LocalRegression);

    swarm_remap.set_remap_var_names(remap_fields, remap_fields, estimator,
                                    basis, op, domains_, data);
    swarm_remap.run();

    // Check the answer
    double toterr = 0.;
#ifdef ENABLE_DEBUG
    auto& cell_field1 = target_state_one.get<Field>("celldata")->get_data();
    auto& node_field1 = target_state_one.get<Field>("nodedata")->get_data();
#endif
    auto& cell_field2 = target_state_two.get<Field>("celldata")->get_data();
    auto& node_field2 = target_state_two.get<Field>("nodedata")->get_data();

    Wonton::Flat_Mesh_Wrapper<double> target_flat_mesh;
    target_flat_mesh.initialize(target_mesh_wrapper);

    if (op != oper::VolumeIntegral) {
      for (int i = 0; i < nb_target_cells; ++i) {
        Wonton::Point<dim> c;
        target_flat_mesh.cell_centroid(i, &c);
        double const value = compute_initial_field(c);
        double const swarm_error = cell_field2[i] - value;

        #if ENABLE_DEBUG
          double const mesh_error  = cell_field1[i] - value;
          //  dump diagnostics for each cell
          switch (dim) {
            case 2: std::printf("cell: %4d coord: (%5.3lf, %5.3lf)", i, c[0], c[1]); break;
            case 3: std::printf("cell: %4d coord: (%5.3lf, %5.3lf, %5.3lf)", i, c[0], c[1], c[2]); break;
            default: break;
          }

          std::printf(" value: %10.6lf,", value);
          std::printf(" mesh: %10.6lf error: %lf", cell_field1[i], mesh_error);
          std::printf(" swarm: %10.6lf error: %lf\n", cell_field2[i], swarm_error);
        #endif
        toterr = std::max(toterr, std::abs(swarm_error));
      }

      #if ENABLE_DEBUG
        std::printf("\n\nL^inf NORM OF MSM CELL ERROR = %lf\n\n", toterr);
      #endif
      ASSERT_LT(toterr, epsilon);

      if (not faceted) {
        toterr = 0.;
        for (int i = 0; i < nb_target_nodes; ++i) {
          Wonton::Point<dim> p;
          target_flat_mesh.node_get_coordinates(i, &p);
          double const value = compute_initial_field(p);
          double const swarm_error = node_field2[i] - value;

          #if ENABLE_DEBUG
            double const mesh_error  = node_field1[i] - value;

            //  dump diagnostics for each node
            switch (dim) {
              case 2: std::printf("node: %4d coord: (%5.3lf, %5.3lf)", i, p[0], p[1]); break;
              case 3: std::printf("node: %4d coord: (%5.3lf, %5.3lf, %5.3lf)", i, p[0], p[1], p[2]); break;
              default: break;
            }

            std::printf(" value: %10.6lf,", value);
            std::printf(" mesh: %10.6lf error: %lf", node_field1[i], mesh_error);
            std::printf(" swarm: %10.6lf error: %lf\n", node_field2[i], swarm_error);
          #endif
          toterr = std::max(toterr, std::abs(swarm_error));
        }

        #if ENABLE_DEBUG
          std::printf("\n\nL^inf NORM OF MSM NODE ERROR = %lf\n\n", toterr);
        #endif
        ASSERT_LT(toterr, epsilon);
      }
    } else {

      for (int i = 0; i < nb_target_cells; ++i) {
        Wonton::Point<dim> c;
        target_flat_mesh.cell_centroid(i, &c);
        double const error = cell_field2[i] - exact[i];

        #if ENABLE_DEBUG
          switch (dim) {
            case 2: std::printf("cell: %4d coord: (%5.3lf, %5.3lf)", i, c[0], c[1]); break;
            case 3: std::printf("cell: %4d coord: (%5.3lf, %5.3lf, %5.3lf)", i, c[0], c[1], c[2]); break;
            default: break;
          }
          std::printf(" exact: %10.6lf,", exact[i]);
          std::printf(" swarm: %10.6lf error: %10.6lf\n", cell_field2[i], error);
        #endif
        toterr = std::max(toterr, std::abs(error));
      }

      #if ENABLE_DEBUG
        std::printf("\n\nL^inf NORM OF MSM OPERATOR ERROR = %lf\n\n", toterr);
      #endif
      ASSERT_LT(toterr, epsilon);
    }
  }

protected:
  // Source and target meshes
  std::shared_ptr<Wonton::Simple_Mesh> source_mesh;
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh;

  // Wrappers for interfacing with the underlying mesh data structures
  Wonton::Simple_Mesh_Wrapper source_mesh_wrapper;
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> source_state;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> target_state_one;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> target_state_two;

  // Operator domains and data
  Wonton::vector<oper::Domain> domains_;
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

double compute_linear_field_2d(Wonton::Point<2> const& centroid) {
  return centroid[0] + centroid[1];
}

double compute_quadratic_field_2d(Wonton::Point<2> const& centroid) {
  return centroid[0] * centroid[0] +
         centroid[1] * centroid[1] +
         centroid[0] * centroid[1];
}

// Methods for computing initial field values
double compute_linear_field_3d(Wonton::Point<3> const& centroid) {
  return centroid[0] + centroid[1] + centroid[2];
}
double compute_quadratic_field_3d(Wonton::Point<3> const& centroid) {
  return centroid[0] * centroid[0] + centroid[1] * centroid[1] +
         centroid[2] * centroid[2] + centroid[0] * centroid[1] +
         centroid[1] * centroid[2] + centroid[2] * centroid[0];
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
    (&compute_linear_field_2d, .75, basis::Linear);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadratic) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, basis::Quadratic);
}

TEST_F(MSMDriverTest2D, 2D1stOrderLinearScatter) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_linear_field_2d, .75, basis::Linear,
     Scatter);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticScatter) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, basis::Quadratic,
     Scatter);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticGatherFaceted) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, basis::Quadratic,
     Gather, oper::LastOperator, true);
}

TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticScatterFaceted) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, basis::Quadratic,
     Scatter, oper::LastOperator, true);
}

TEST_F(MSMDriverTest2D, 2D1stOrderLinearIntegrate) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_linear_field_2d, .75, basis::Linear,
     Gather, oper::VolumeIntegral);
}


TEST_F(MSMDriverTest2D, 2D2ndOrderQuadraticIntegrate) {
  unitTest<Portage::IntersectR2D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 2>
    (&compute_quadratic_field_2d, .75, basis::Quadratic,
     Gather, oper::VolumeIntegral);
}


TEST_F(MSMDriverTest3D, 3D1stOrderLinear) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, basis::Linear);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadratic) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, basis::Quadratic);
}

TEST_F(MSMDriverTest3D, 3D1stOrderLinearScatter) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, basis::Linear,
     Scatter);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticScatter) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, basis::Quadratic,
     Scatter);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticGatherFaceted) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, basis::Quadratic,
     Gather, oper::LastOperator, true);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticScatterFaceted) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, basis::Quadratic,
     Scatter, oper::LastOperator, true);
}

TEST_F(MSMDriverTest3D, 3D1stOrderLinearIntegrate) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, basis::Linear,
     Gather, oper::VolumeIntegral);
}


}  // namespace
