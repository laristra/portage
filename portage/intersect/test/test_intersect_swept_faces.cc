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
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"

/**
 * @brief Fixture class for swept face moments computation tests.
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
                           [](double previous, auto const& moment) {
                             return previous + moment.weights[0];
                           });
  }

  /**
   *
   * @param id
   * @param moments
   * @return
   */
  double compute_contribution(int id, std::vector<Wonton::Weights_t> const& moments) const {
    double contrib = source_mesh_wrapper.cell_volume(id);
    for (auto const& moment : moments) {
      if (moment.entityID == id) {
        contrib += moment.weights[0];
      }
    }
    return contrib;
  }

protected:
  // useful constants and aliases
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
 *     ___________       displacement vector: (-1,-1) 
 *  ..|...|...|.. |      source mesh: plain            
 *  : |_:_|_:_|_:_|      target mesh: dotted           
 *  :.|...|...|.. |
 *  : |_:_|_:_|_:_|
 *  :.|...|...|.. |
 *  : |_:_|_:_|_:_|
 *  :...:...:...:
 */
class IntersectSweptBackward : public IntersectSweptBase {
protected:
  IntersectSweptBackward() : IntersectSweptBase(-1, -1, 5, 5) {}
};

TEST_F(IntersectSweptForward, MomentsCheck) {

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

  /* for interior cell case, we have:
   * - two contributing neighbor cells with non-zero weights (5:top, 7:right)
   * - a reconstructed cell area perfectly preserved after advection.
   * - self-contribution is reduced to zero.
   */
  double source_area = source_mesh_wrapper.cell_volume(internal_cell);
  double target_area = source_area + compute_swept_area(weights_internal);
  double self_contrib = compute_contribution(internal_cell, weights_internal);
  bool first_source_weight = true;

  ASSERT_EQ(weights_internal.size(), unsigned(4));
  ASSERT_DOUBLE_EQ(target_area, source_area);
  ASSERT_DOUBLE_EQ(self_contrib, 0.0);

  for (auto const& moments : weights_internal) {
    #if DEBUG
      std::cout << "forward::internal_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
    #endif
    ASSERT_DOUBLE_EQ(std::abs(moments.weights[0]), unit_face_area);
    switch(moments.entityID) {
      case internal_cell:
        if (first_source_weight) {
          ASSERT_DOUBLE_EQ(moments.weights[1], 3.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 2.5);
        } else {
          ASSERT_DOUBLE_EQ(moments.weights[1], 2.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 3.5);
        }
        first_source_weight = false;
        break;
      case 5:
        ASSERT_DOUBLE_EQ(moments.weights[1], 3.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 4.5);
        break;
      case 7:
        ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 3.5);
        break;
      default: FAIL() << "forward::internal: unexpected moment entity index";
    }
  }

  /* check boundary cell case:
   * - we are not conservative anymore since swept faces lying outside
   *   the source mesh are not taken into account.
   * - we have a unique contributing neighbor (8:top).
   * - the source cell self-contribution is still zero.
   */
  source_area = source_mesh_wrapper.cell_volume(boundary_cell);
  target_area = source_area + compute_swept_area(weights_boundary);
  self_contrib = compute_contribution(boundary_cell, weights_boundary);
  first_source_weight = true;

  ASSERT_EQ(weights_boundary.size(), unsigned(3));
  ASSERT_DOUBLE_EQ(target_area, 0.5 * source_area);
  ASSERT_DOUBLE_EQ(self_contrib, 0.0);

  for (auto const& moments : weights_boundary) {
    #if DEBUG
      std::cout << "forward::boundary_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
    #endif
    ASSERT_DOUBLE_EQ(std::abs(moments.weights[0]), unit_face_area);
    switch(moments.entityID) {
      case boundary_cell:
        if (first_source_weight) {
          ASSERT_DOUBLE_EQ(moments.weights[1], 5.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 2.5);
        } else {
          ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 3.5);
        }
        first_source_weight = false;
        break;
      case 8:
        ASSERT_DOUBLE_EQ(moments.weights[1], 5.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 4.5);
        break;
      default: FAIL() << "forward::boundary: unexpected moment entity index";
    }
  }

  /* check corner cell case:
   * - the source cell self-contribution is still zero.
   * - we have no more contributing neighbor since all swept faces are
   *   lying outside the source mesh and their volume are not extrapolated.
   */
  source_area = source_mesh_wrapper.cell_volume(corner_cell);
  target_area = source_area + compute_swept_area(weights_corner);
  self_contrib = compute_contribution(corner_cell, weights_corner);
  first_source_weight = true;

  ASSERT_EQ(weights_corner.size(), unsigned(2));
  ASSERT_DOUBLE_EQ(target_area, 0.0);
  ASSERT_DOUBLE_EQ(self_contrib, 0.0);

  for (auto const& moments : weights_corner) {
    #if DEBUG
      std::cout << "forward::corner_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
    #endif
    ASSERT_DOUBLE_EQ(std::abs(moments.weights[0]), unit_face_area);
    switch(moments.entityID) {
      case corner_cell:
        if (first_source_weight) {
          ASSERT_DOUBLE_EQ(moments.weights[1], 5.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 4.5);
        } else {
          ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 5.5);
        }
        first_source_weight = false;
        break;
      default: FAIL() << "forward::boundary: unexpected moment entity index";
    }
  }
}


TEST_F(IntersectSweptBackward, MomentsCheck) {

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

  /* we should have the same results as the previous test for interior case:
   * - two contributing neighbors cells with non-zero weights.
   * - reconstructed cell area is perfectly preserved after advection.
   * - no cell self-contribution.
   */
  double source_area = source_mesh_wrapper.cell_volume(internal_cell);
  double target_area = source_area + compute_swept_area(weights_internal);
  double self_contrib = compute_contribution(internal_cell, weights_internal);
  bool first_source_weight = true;

  ASSERT_EQ(weights_internal.size(), unsigned(4));
  ASSERT_DOUBLE_EQ(source_area, target_area);
  ASSERT_DOUBLE_EQ(self_contrib, 0.0);

  for (auto const& moments : weights_internal) {
    #if DEBUG
      std::cout << "backward::internal_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
    #endif
    ASSERT_DOUBLE_EQ(std::abs(moments.weights[0]), unit_face_area);
    switch(moments.entityID) {
      case internal_cell:
        if (first_source_weight) {
          ASSERT_DOUBLE_EQ(moments.weights[1], 3.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 2.5);
        } else {
          ASSERT_DOUBLE_EQ(moments.weights[1], 2.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 3.5);
        }
        first_source_weight = false;
        break;
      case 1:
        ASSERT_DOUBLE_EQ(moments.weights[1], 1.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 2.5);
        break;
      case 3:
        ASSERT_DOUBLE_EQ(moments.weights[1], 2.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 1.5);
        break;
      default: FAIL() << "backward::internal: unexpected moment entity index";
    }
  }

  /* for the boundary cell case:
   * - all swept faces are covered by the source mesh this time,
   *   so we are in the very same situation as for internal cell.
   * - we have two contributing neighbors (4:left, 6:bottom).
   * - no cell self-contribution.
   */
  source_area = source_mesh_wrapper.cell_volume(boundary_cell);
  target_area = source_area + compute_swept_area(weights_boundary);
  self_contrib = compute_contribution(boundary_cell, weights_boundary);
  first_source_weight = true;

  ASSERT_EQ(weights_boundary.size(), weights_internal.size());
  ASSERT_DOUBLE_EQ(source_area, target_area);
  ASSERT_DOUBLE_EQ(self_contrib, 0.0);

  for (auto const& moments : weights_boundary) {
    #if DEBUG
      std::cout << "backward::boundary_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
    #endif
    ASSERT_DOUBLE_EQ(std::abs(moments.weights[0]), unit_face_area);
    switch(moments.entityID) {
      case boundary_cell:
        if (first_source_weight) {
          ASSERT_DOUBLE_EQ(moments.weights[1], 5.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 2.5);
        } else {
          ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 3.5);
        }
        first_source_weight = false;
        break;
      case 4:
        ASSERT_DOUBLE_EQ(moments.weights[1], 3.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 2.5);
        break;
      case 6:
        ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 1.5);
        break;
      default: FAIL() << "backward::boundary: unexpected moment entity index";
    }
  }

  /* for the corner cell case (as opposed to the forward case):
   * - all swept faces are covered by the source mesh so we get
   *   the same moments as for an internal cell.
   * - the reconstructed cell area is preserved after advection.
   * - we have two contributing neighbors (5:left, 7:bottom).
   * - the source cell self-contribution is still zero.
   */
  source_area = source_mesh_wrapper.cell_volume(corner_cell);
  target_area = source_area + compute_swept_area(weights_corner);
  self_contrib = compute_contribution(corner_cell, weights_corner);
  first_source_weight = true;

  ASSERT_EQ(weights_corner.size(), weights_internal.size());
  ASSERT_DOUBLE_EQ(source_area, target_area);
  ASSERT_DOUBLE_EQ(self_contrib, 0.0);

  for (auto const& moments : weights_corner) {
    #if DEBUG
      std::cout << "backward::corner_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << moments.weights[1] <<", "<< moments.weights[2] << std::endl;
    #endif
    ASSERT_DOUBLE_EQ(std::abs(moments.weights[0]), unit_face_area);
    switch(moments.entityID) {
      case corner_cell:
        if (first_source_weight) {
          ASSERT_DOUBLE_EQ(moments.weights[1], 5.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 4.5);
        } else {
          ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
          ASSERT_DOUBLE_EQ(moments.weights[2], 5.5);
        }
        first_source_weight = false;
        break;
      case 5:
        ASSERT_DOUBLE_EQ(moments.weights[1], 3.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 4.5);
        break;
      case 7:
        ASSERT_DOUBLE_EQ(moments.weights[1], 4.5);
        ASSERT_DOUBLE_EQ(moments.weights[2], 3.5);
        break;
      default: FAIL() << "backward::corner: unexpected moment entity index";
    }
  }
}
