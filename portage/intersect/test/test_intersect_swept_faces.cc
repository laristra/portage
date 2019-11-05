/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "gtest/gtest.h"

#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_swept_faces.h"
#include "portage/driver/coredriver.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"

/**
 * @brief Fixture class for swept face intersection tests.
 *
 * Here, we consider a planar cartesian grid which is advected
 * by a unique displacement vector and we aim to check the volume
 * of each swept face for a given source cell. Notice that mesh
 * topology remains unchanged, hence source and target cells
 * and faces indices are kept.
 */
class IntersectSweptBase : public testing::Test {

protected:
  using Intersector = Portage::IntersectSweptFace<2, Wonton::Entity_kind::CELL,
                                                  Wonton::Jali_Mesh_Wrapper,
                                                  Wonton::Jali_State_Wrapper,
                                                  Wonton::Jali_Mesh_Wrapper>;
public:
  /**
   * @brief Disabled default constructor
   */
  IntersectSweptBase() = delete;

  /**
   * @brief Setup each test-case.
   *
   * It initializes both source and target meshes and states,
   * then computes and assigns a density field on source mesh.
   */
  IntersectSweptBase(double x0, double y0, double x1, double y1)
    : source_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 6.0, 6.0, 3, 3)),
      target_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(x0, y0, x1, y1, 3, 3)),
      source_state(Jali::State::create(source_mesh)),
      target_state(Jali::State::create(target_mesh)),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state_wrapper(*source_state),
      target_state_wrapper(*target_state)
  {
    num_tols.use_default();
  }

  /**
   * @brief Retrieve all cells incident to the faces of a given cell.
   *
   * @param cell: the current cell.
   * @return a list of cell neighbors incident to the faces of that cell.
   */
  std::vector<int> search(int current) const {

    std::vector<int> list, faces, dirs, cells;
    source_mesh_wrapper.cell_get_faces_and_dirs(current, &faces, &dirs);

    list.emplace_back(current);
    for (auto &&f : faces) {
      cells.clear();
      source_mesh_wrapper.face_get_cells(f, Wonton::Entity_type::ALL, &cells);
      if (cells.size() == 2) {
        int const index = (cells[0] == current ? 1 : 0);
        list.emplace_back(cells[index]);
      }
    }
    return list;
  }


  /**
   *
   * @param moments
   * @return
   */
  double compute_swept_area(std::vector<Wonton::Weights_t> const& moments) const {
    return std::accumulate(moments.begin(), moments.end(), 0.0,
      [](double prec, auto const& moment) {
        return prec + moment.weights[0];
      });
  }

protected:
  // useful constants and aliases
  static constexpr double upper_bound = std::numeric_limits<double>::max();
  static constexpr double lower_bound = std::numeric_limits<double>::min();
  static constexpr double epsilon = 1.E-10;
  static constexpr auto CELL = Wonton::Entity_kind::CELL;

  Portage::NumericTolerances_t num_tols;

  // source and target meshes and states
  std::shared_ptr<Jali::Mesh> source_mesh;
  std::shared_ptr<Jali::Mesh> target_mesh;
  std::shared_ptr<Jali::State> source_state;
  std::shared_ptr<Jali::State> target_state;

  // wrappers for interfacing with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper;
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper;
  Wonton::Jali_State_Wrapper source_state_wrapper;
  Wonton::Jali_State_Wrapper target_state_wrapper;

  /* since target grid is advected using a unique displacement vector
   * we expect to have a constant area for each swept face.
   */
  double const unit_face_area = 2.;
};

/**
 * @brief Fixture class for intersection moment computation tests
 *        when target cells are swept forward like below.
 *
 *    .............      displacement vector: (1,1)
 *   _:___:___:_  :      source mesh: plain
 *  | :.|.:.|.:.|.:      target mesh: dotted
 *  |_:_|_:_|_:_| :
 *  | :.|.:.|.:.|.:
 *  |_:_|_:_|_:_| :
 *  | :.|.:.|.:.|.:
 *  |___|___|___|
 *  0   2   4   6
 */
class IntersectSweptForward : public IntersectSweptBase {
protected:
  IntersectSweptForward() : IntersectSweptBase(1, 1, 7, 7) {}
};

/**
 * @brief Fixture class for intersection moment computation tests
 *        when target cells are swept backward like below.
 *
 *     ___________
 *  ..|...|...|.. |       displacement vector: (-1,-1)
 *  : |_:_|_:_|_:_|       source mesh: plain
 *  :.|...|...|.. |       target mesh: dotted
 *  : |_:_|_:_|_:_|
 *  :.|...|...|.. |
 *  : |_:_|_:_|_:_|
 *  :...:...:...:
 */
class IntersectSweptBackward : public IntersectSweptBase {
protected:
  IntersectSweptBackward() : IntersectSweptBase(-1, -1, 5, 5) {}
};

TEST_F(IntersectSweptForward, AreaCheck) {

  Intersector intersector(source_mesh_wrapper,
                          source_state_wrapper,
                          target_mesh_wrapper,
                          num_tols);

  /* pick internal, boundary and corner source cells.
   *    .............
   *   _:___:___:_  :     cell index ordering:
   *  | :.|.:.|.:.|.:     2 5 8
   *  |_:_|_:_|_:_| :     1 4 7
   *  | :.|.:.|.:.|.:     0 3 6
   *  |_:_|_:_|_:_| :
   *  | :.|.:.|.:.|.:
   *  |___|___|___|
   */
  int const internal_cell = 4;
  int const boundary_cell = 7;
  int const corner_cell   = 8;

  // search for candidate cells and compute moments of intersection
  auto weights_internal = intersector(internal_cell, search(internal_cell));
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));
  auto weights_corner   = intersector(corner_cell, search(corner_cell));

  /* check interior cell case:
   * - we have two contributing cells with non-zero weights.
   * - the target cell area is perfectly preserved after advection.
   * - the source cell self-contribution is zero.
   */
  ASSERT_EQ(weights_internal.size(), unsigned(2));

  // check area preservation
  double source_area = source_mesh_wrapper.cell_volume(internal_cell);
  double target_area = compute_swept_area(weights_internal);
  ASSERT_NEAR(source_area, target_area, epsilon);

  for (auto const& moments : weights_internal) {
    ASSERT_NEAR(moments.weights[0], unit_face_area, epsilon);
  }

  /* check boundary cell case:
   * - we are not conservative anymore since swept faces lying outside
   *   the source mesh are not taken into account.
   * - we have a unique contributing neighbor.
   * - the source cell self-contribution is still zero.
   */
  ASSERT_EQ(weights_boundary.size(), unsigned(1));

  source_area = source_mesh_wrapper.cell_volume(boundary_cell);
  target_area = compute_swept_area(weights_boundary);
  ASSERT_NEAR(source_area, 2 * target_area, epsilon);

  for (auto const& moments : weights_boundary) {
    ASSERT_NEAR(moments.weights[0], unit_face_area, epsilon);
  }

  /* check corner cell case:
   * - the source cell self-contribution is still zero.
   * - we have no more contributing neighbor since all swept faces are
   *   lying outside the source mesh and their volume are not extrapolated.
   */
  ASSERT_TRUE(weights_corner.empty());
}


TEST_F(IntersectSweptBackward, AreaCheck) {

  Intersector intersector(source_mesh_wrapper,
                          source_state_wrapper,
                          target_mesh_wrapper,
                          num_tols);

  /* pick internal, boundary and corner source cells.
   *
   *     ___________      cell indices ordering:
   *  ..|...|...|.. |     2 5 8       
   *  : |_:_|_:_|_:_|     1 4 7                                 
   *  :.|...|...|.. |     0 3 6          
   *  : |_:_|_:_|_:_|
   *  :.|...|...|.. |
   *  : |_:_|_:_|_:_|
   *  :...:...:...:
   */
  int const internal_cell = 4;
  int const boundary_cell = 7;
  int const corner_cell   = 8;

  // search for candidate cells and compute moments of intersection
  auto weights_internal = intersector(internal_cell, search(internal_cell));
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));
  auto weights_corner   = intersector(corner_cell, search(corner_cell));

  /* we should have the same results as the previous test for interior
   * cell case, that is:
   * - two contributing cells with non-zero weights.
   * - a target cell area is perfectly preserved after advection.
   * - a source cell self-contribution reduced to zero.
   */
  ASSERT_EQ(weights_internal.size(), unsigned(2));

  // check area preservation
  double source_area = source_mesh_wrapper.cell_volume(internal_cell);
  double target_area = compute_swept_area(weights_internal);
  ASSERT_NEAR(source_area, target_area, epsilon);

  for (auto const& moments : weights_internal) {
    ASSERT_NEAR(moments.weights[0], unit_face_area, epsilon);
  }

  /* for the boundary cell case:
   * - all swept faces are covered by the source mesh this time,
   *   so we are in the very same situation as for internal cell.
   * - we have two contributing neighbors.
   * - we have one source cell self-contribution this time.
   */
  ASSERT_EQ(weights_boundary.size(), weights_internal.size());

  source_area = source_mesh_wrapper.cell_volume(boundary_cell);
  target_area = compute_swept_area(weights_boundary);
  ASSERT_NEAR(source_area, target_area, epsilon);

  for (auto const& moments : weights_boundary) {
    ASSERT_NEAR(moments.weights[0], unit_face_area, epsilon);
    switch(moments.entityID) {
      case 4:
        ASSERT_NEAR(moments.weights[1], 3.5, epsilon);
        ASSERT_NEAR(moments.weights[2], 2.5, epsilon);
        break;
      case 6:
        ASSERT_NEAR(moments.weights[1], 4.5, epsilon);
        ASSERT_NEAR(moments.weights[2], 1.5, epsilon);
        break;
      default: FAIL() << "Unexpected moment entity index";
    }
//    std::cout << "boundary_swept_centroid["<< moments.entityID <<"]: ";
//    std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
  }

  /* check corner cell case:
   * - the source cell self-contribution is still zero.
   * - we have no more contributing neighbor since all swept faces are
   *   lying outside the source mesh and their volume are not extrapolated.
   */
  //ASSERT_TRUE(weights_corner.empty());
}
