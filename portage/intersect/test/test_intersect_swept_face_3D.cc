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
#include "portage/intersect/intersect_swept_face.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"
#ifdef HAVE_TANGRAM
  #include "tangram/intersect/split_r3d.h"
  #include "tangram/reconstruct/MOF.h"
#endif
/**
 * @brief Fixture class for swept volume moments computation tests.
 *
 * Here, we consider a 3D cartesian grid which is advected
 * by a unique displacement vector and we aim to check the volume
 * of each swept region for a given source cell. Notice that mesh
 * topology remains unchanged, hence source and target cells
 * and faces indices are kept.
 */
class IntersectSweptBase3D : public testing::Test {

protected:
#ifdef HAVE_TANGRAM
  using Intersector = Portage::IntersectSweptFace3D<Wonton::Entity_kind::CELL,
                                                    Wonton::Jali_Mesh_Wrapper,
                                                    Wonton::Jali_State_Wrapper,
                                                    Wonton::Jali_Mesh_Wrapper,
                                                    Tangram::MOF,
                                                    Tangram::SplitR3D,
                                                    Tangram::ClipR3D>;
#else
  using Intersector = Portage::IntersectSweptFace3D<Wonton::Entity_kind::CELL,
                                                    Wonton::Jali_Mesh_Wrapper,
                                                    Wonton::Jali_State_Wrapper,
                                                    Wonton::Jali_Mesh_Wrapper>;
#endif
public:
  /**
   * @brief Disabled default constructor
   */
  IntersectSweptBase3D() = delete;

  /**
   * @brief Setup each test-case.
   *
   * It initializes both source and target meshes and states,
   * then computes and assigns a density field on source mesh.
   */
  IntersectSweptBase3D(double x0, double y0, double z0, double x1, double y1, double z1)
    : source_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 6.0, 6.0, 6.0, 3, 3, 3)),
      target_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(x0, y0, z0, x1, y1, z1, 3, 3, 3)),
      source_state(Jali::State::create(source_mesh)),
      target_state(Jali::State::create(target_mesh)),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state_wrapper(*source_state),
      target_state_wrapper(*target_state)
  {
    num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<3>;
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
   * @brief Compute the total swept region volume.
   *
   * @param moments: list of previously computed moments.
   * @return: the total area of all swept faces.
   */
  double compute_swept_volume(std::vector<Wonton::Weights_t> const& moments) const {
    double value = 0.;
    for (auto const& moment : moments) {
      value += moment.weights[0];
    }
    return (std::abs(value) < num_tols.min_absolute_volume ? 0.0 : value);
  }

  /**
   * @brief Compute the given cell volume contribution.
   *
   * @param id: cell index.
   * @param moments: list of previously computed moments.
   * @return the area contribution of the given cell.
   */
  static double compute_contribution(int id, std::vector<Wonton::Weights_t> const& moments) {
    double contrib = 0.;
    for (auto const& moment : moments) {
      if (moment.entityID == id) {
        contrib += moment.weights[0];
      }
    }
    return contrib;
  }

  /**
   * @brief Deduce the centroid coordinates from current cell moments.
   *
   * @param moment: current cell moments
   * @return its centroid coordinates
   */
  static Wonton::Point<3> deduce_centroid(Wonton::Weights_t const& moment) {
    return Wonton::createP3(moment.weights[1] / moment.weights[0],
                            moment.weights[2] / moment.weights[0],
                            moment.weights[3] / moment.weights[0]);
  }

protected:
  // numerical tolerances, and threshold for floating-point comparison
  Portage::NumericTolerances_t num_tols;
  double const epsilon = 1.E-14;

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
   * we expect to have a constant volume for each swept region.
   */
  double const unit_region_volume = 4.0;
  int const nb_hex_faces  = 6;

  // enable or disable debug prints
  bool verbose = false;
};


/**
 * @brief Fixture class for intersection moment computation tests
 *        when target cells are swept forward.
 *           .......
 *          . :   .:
 *         . :   . :          displacement vector: (1,1,1)
 *       ....:...  :          source cell: plain
 *     __:___:..:..:          target cell: dotted
 *    /| :  /|  : .
 *   / | : / |  :.
 *  /__|_:/..|..:       z
 * |   |__|__|          |  y
 * |  /   |  /          | /
 * | /    | /           |/___ x
 * |/_____|/
 *
 */
class IntersectSweptForward3D : public IntersectSweptBase3D {
protected:
  IntersectSweptForward3D() : IntersectSweptBase3D(1, 1, 1, 7, 7, 7) {}
};

/**
 * @brief Fixture class for intersection moment computation tests
 *        when target cells are swept backward.
 *           ______
 *          /|    /|
 *         / |   / |          displacement vector: (-1,-1,-1)
 *        /__|__/  |          source cell: plain
 *     ..|...:__|__|          target cell: dotted
 *    .: |  .:  |  /
 *   . : | . :  | /
 *  ...:.|.__:__|/       z
 * :   :..:..:          |  y
 * :  .   :  .          | /
 * : .    : .           |/___ x
 * :......:.
 *
 */
class IntersectSweptBackward3D : public IntersectSweptBase3D {
protected:
  IntersectSweptBackward3D() : IntersectSweptBase3D(-1, -1, -1, 5, 5, 5) {}
};

/**
 * @brief Fixture class for intersection moment computation tests
 *        when target cells are swept only along x-axis.
 *                    
 *                             displacement vector: (1,0,0)
 *     ____......              source cell: plain
 *    /|  .:/| .:              target cell: dotted
 *   / | . / |. :
 *  /__|../:.|  :
 * |   :__|:_:..:       z
 * |  /:  |  :  .      |  y
 * | / : .| /: .       | /
 * |/__:..|..:.        |/___ x
 *
 */
class IntersectSweptOneAxis3D : public IntersectSweptBase3D {
protected:
  IntersectSweptOneAxis3D() : IntersectSweptBase3D(1, 0, 0, 7, 6, 6) {}
};

TEST_F(IntersectSweptForward3D, MomentsCheck) {

  Intersector intersector(source_mesh_wrapper,
                          source_state_wrapper,
                          target_mesh_wrapper,
                          num_tols);

  /* pick internal, boundary and corner source cells.
   * swept region configuration for forward displacement:
   *
   *   source hex        target hex           frontal face
   *                     7'......6'          swept polyhedron:
   *                      .:    .:             4'......5'
   *    7______6         . :   . :             /:    /:
   *    /|    /|      4'...:..5' :            / :   / :
   *   / |   / |       :   :..:..2'          /  :  /  :
   * 4 __|__5  |       :  .   :  .         4____:_5...:
   * |   3__|__2       : .    : .          |   /0'|  /1'
   * |  /   |  /       :......:            |  /   | /
   * | /    | /        0'     1'           | /    |/
   * |/_____|/                             |/_____/
   * 0      1                              0      1
   */
  int const internal_cell = 13;
  int const boundary_cell = 25;
  int const corner_cell   = 26;

  // search for candidate cells and compute moments of intersection
  auto weights_internal = intersector(internal_cell, search(internal_cell));
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));
  auto weights_corner   = intersector(corner_cell, search(corner_cell));

  /* for interior hexahedral cell, we have:
   * - 3 positive swept volumes attached to cells 14:top|16:right|22:forward (V=-12).
   * - 3 negative swept volume attached to the source cell itself (V=12).
   * - the absolute volume of the source cell (V=8).
   * hence the cell self-contribution is -4.
   */
  double source_volume = source_mesh_wrapper.cell_volume(internal_cell);
  double target_volume = compute_swept_volume(weights_internal);
  double self_contrib  = compute_contribution(internal_cell, weights_internal);
  int nb_self_weights = 0;

  ASSERT_EQ(weights_internal.size(), unsigned(nb_hex_faces + 1));
  ASSERT_DOUBLE_EQ(source_volume, target_volume);
  ASSERT_DOUBLE_EQ(self_contrib, -unit_region_volume);

  for (auto const& moments : weights_internal) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
    #if DEBUG
    if (verbose) {
      std::cout << "forward::internal_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << centroid[0] <<", "<< centroid[1] << ", " << centroid[2];
      std::cout << std::endl;
    }
    #endif
    switch (moments.entityID) {
      case internal_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 3.0);
            ASSERT_DOUBLE_EQ(centroid[1], 3.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 3.5);
            ASSERT_DOUBLE_EQ(centroid[1], 3.5);
            ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
          case 3:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 3.5);
            ASSERT_DOUBLE_EQ(centroid[1], 2.5);
            ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
          case 4:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 2.5);
            ASSERT_DOUBLE_EQ(centroid[1], 3.5);
            ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
          default: FAIL() << "forward::internal: invalid self weights count";
        } break;
      case 14:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 3.5);
        ASSERT_DOUBLE_EQ(centroid[1], 3.5);
        ASSERT_DOUBLE_EQ(centroid[2], 4.5); break;
      case 16:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 3.5);
        ASSERT_DOUBLE_EQ(centroid[1], 4.5);
        ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
      case 22:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 4.5);
        ASSERT_DOUBLE_EQ(centroid[1], 3.5);
        ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
      default: FAIL() << "forward::internal: unexpected moment entity index";
    }
  }

  /* check boundary hex cell case:
   * - the total swept volume is not conserved anymore since swept region lying
   *   outside the source mesh are not taken into account.
   * - we have a unique contributing neighbor (26:top).
   * - the source cell self-contribution is still -4.
   */
  source_volume = source_mesh_wrapper.cell_volume(boundary_cell);
  target_volume = compute_swept_volume(weights_boundary);
  self_contrib  = compute_contribution(boundary_cell, weights_boundary);
  nb_self_weights = 0;

  ASSERT_EQ(weights_boundary.size(), unsigned(nb_hex_faces - 1));
  ASSERT_NEAR(target_volume, 0.0, 1.0e-15);
  ASSERT_DOUBLE_EQ(self_contrib, -unit_region_volume);

  for (auto const& moments : weights_boundary) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "forward::boundary_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << centroid[0] <<", "<< centroid[1] << ", " << centroid[2];
      std::cout << std::endl;
    }
#endif
    switch (moments.entityID) {
      case boundary_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.0);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.5);
            ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
          case 3:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.5);
            ASSERT_DOUBLE_EQ(centroid[1], 4.5);
            ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
          case 4:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 4.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.5);
            ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
          default: FAIL() << "forward::boundary: invalid self weights count";
        } break;
      case 26:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 5.5);
        ASSERT_DOUBLE_EQ(centroid[1], 5.5);
        ASSERT_DOUBLE_EQ(centroid[2], 4.5); break;
      default: FAIL() << "forward::boundary: unexpected moment entity index";
    }
  }

  /* check corner hex cell case:
   * - the source cell self-contribution is still -4.
   * - we have no more contributing neighbor since all swept faces are
   *   lying outside the source mesh and their volume are not extrapolated.
   * - the total swept volume became negative since it corresponds to the self-
   *   contribution in this case.
   */
  target_volume = compute_swept_volume(weights_corner);
  self_contrib  = compute_contribution(corner_cell, weights_corner);
  nb_self_weights = 0;

  ASSERT_EQ(weights_corner.size(), unsigned(nb_hex_faces - 2));
  ASSERT_DOUBLE_EQ(self_contrib, -unit_region_volume);
  ASSERT_DOUBLE_EQ(target_volume, self_contrib);

  for (auto const& moments : weights_corner) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "forward::corner_swept_centroid[" << moments.entityID << "]: ";
      std::cout << centroid[0] << ", " << centroid[1] << ", " << centroid[2];
      std::cout << std::endl;
    }
#endif
    switch(moments.entityID) {
      case corner_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.0);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 5.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.5);
            ASSERT_DOUBLE_EQ(centroid[2], 4.5); break;
          case 3:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.5);
            ASSERT_DOUBLE_EQ(centroid[1], 4.5);
            ASSERT_DOUBLE_EQ(centroid[2], 5.5); break;
          case 4:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 4.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.5);
            ASSERT_DOUBLE_EQ(centroid[2], 5.5); break;
          default: FAIL() << "forward::corner: invalid self weights count";
        } break;
      default: FAIL() << "forward::corner: unexpected moment entity index";
    }
  }
}

TEST_F(IntersectSweptBackward3D, MomentsCheck) {

  Intersector intersector(source_mesh_wrapper,
                          source_state_wrapper,
                          target_mesh_wrapper,
                          num_tols);

  /* pick internal, boundary and corner source cells.
   * swept region configuration for forward displacement:
   *
   *     target hex      source hex        frontal face
   *                       7______6        swept polyhedron:
   *                       /|    /|          4'______5'
   *     7'......6'       / |   / |           .|    .|
   *      .:    .:      4 __|__5  |          . |   . |
   *     . :   . :      |   3__|__2         .  |  .  |
   *  4'...:..5' :      |  /   |  /       4....|_5___|
   *   :   :..:..2'     | /    | /        :   .0':  .1'
   *   :  .   :  .      |/_____|/         :  .   : .
   *   : .    : .       0      1          : .    :.
   *   :......:                           :......:
   *   0'     1'                          0      1
   */
  int const internal_cell = 13;
  int const boundary_cell = 25;
  int const corner_cell   = 26;

  // search for candidate cells and compute moments of intersection
  auto weights_internal = intersector(internal_cell, search(internal_cell));
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));
  auto weights_corner   = intersector(corner_cell, search(corner_cell));

  /* we should have the same configuration as the previous tests for interior
   * hexahedral cell, that is:
   * - 3 positive swept volumes attached to cells 4:backward|10:left|12:bottom (V=-12).
   * - 3 negative swept volume attached to the source cell itself (V=12).
   * - the absolute volume of the source cell (V=8).
   * - the self-contribution is still -4.
   */
  double source_volume = source_mesh_wrapper.cell_volume(internal_cell);
  double target_volume = compute_swept_volume(weights_internal);
  double self_contrib  = compute_contribution(internal_cell, weights_internal);
  int nb_self_weights = 0;

  ASSERT_EQ(weights_internal.size(), unsigned(nb_hex_faces + 1));
  ASSERT_NEAR(source_volume, target_volume, epsilon);
  ASSERT_DOUBLE_EQ(self_contrib, -unit_region_volume);

  for (auto const& moments : weights_internal) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "backward::internal_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << centroid[0] <<", "<< centroid[1] << ", " << centroid[2];
      std::cout << std::endl;
    }
#endif
    switch (moments.entityID) {
      case internal_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 3.0);
            ASSERT_DOUBLE_EQ(centroid[1], 3.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 2.5);
            ASSERT_DOUBLE_EQ(centroid[1], 2.5);
            ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
          case 3:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 3.5);
            ASSERT_DOUBLE_EQ(centroid[1], 2.5);
            ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
          case 4:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 2.5);
            ASSERT_DOUBLE_EQ(centroid[1], 3.5);
            ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
          default: FAIL() << "backward::internal: invalid self weights count";
        } break;
      case 4:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 1.5);
        ASSERT_DOUBLE_EQ(centroid[1], 2.5);
        ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
      case 10:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 2.5);
        ASSERT_DOUBLE_EQ(centroid[1], 1.5);
        ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
      case 12:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 2.5);
        ASSERT_DOUBLE_EQ(centroid[1], 2.5);
        ASSERT_DOUBLE_EQ(centroid[2], 1.5); break;
      default: FAIL() << "backward::internal: unexpected moment entity index";
    }
  }

  /* for the boundary hex cell case:
   * - all swept volumes are covered by the source mesh this time,
   *   so we are in the very same situation as for internal cell.
   * - we have three contributing neighbors (16:backward|22:left|24:bottom).
   * - the self-contribution is still -4.
   */
  source_volume = source_mesh_wrapper.cell_volume(boundary_cell);
  target_volume = compute_swept_volume(weights_boundary);
  self_contrib  = compute_contribution(boundary_cell, weights_boundary);
  nb_self_weights = 0;

  ASSERT_EQ(weights_boundary.size(), unsigned(nb_hex_faces + 1));
  ASSERT_NEAR(source_volume, target_volume, epsilon);
  ASSERT_NEAR(self_contrib, -unit_region_volume, epsilon);

  for (auto const& moments : weights_boundary) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "backward::boundary_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << centroid[0] <<", "<< centroid[1] << ", " << centroid[2];
      std::cout << std::endl;
    }
#endif
    switch (moments.entityID) {
      case boundary_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.0);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 4.5);
            ASSERT_DOUBLE_EQ(centroid[1], 4.5);
            ASSERT_DOUBLE_EQ(centroid[2], 3.5); break;
          case 3:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.5);
            ASSERT_DOUBLE_EQ(centroid[1], 4.5);
            ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
          case 4:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 4.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.5);
            ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
          default: FAIL() << "backward::boundary: invalid self weights count";
        } break;
      case 16:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 3.5);
        ASSERT_DOUBLE_EQ(centroid[1], 4.5);
        ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
      case 22:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 4.5);
        ASSERT_DOUBLE_EQ(centroid[1], 3.5);
        ASSERT_DOUBLE_EQ(centroid[2], 2.5); break;
      case 24:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 4.5);
        ASSERT_DOUBLE_EQ(centroid[1], 4.5);
        ASSERT_DOUBLE_EQ(centroid[2], 1.5); break;
      default: FAIL() << "backward::boundary: unexpected moment entity index";
    }
  }

  /* for the corner cell case (as opposed to the forward case):
   * - all swept regions are covered by the source mesh so we get
   *   the same moments as for an internal cell.
   * - the reconstructed cell volume is preserved after displacement.
   * - we have three contributing neighbors (17:backward|23:left|25:bottom).
   * - the source cell self-contribution is still -4.
   */
  source_volume = source_mesh_wrapper.cell_volume(corner_cell);
  target_volume = compute_swept_volume(weights_corner);
  self_contrib  = compute_contribution(corner_cell, weights_corner);
  nb_self_weights = 0;

  ASSERT_EQ(weights_corner.size(), weights_internal.size());
  ASSERT_NEAR(source_volume, target_volume, epsilon);
  ASSERT_NEAR(self_contrib, -unit_region_volume, epsilon);

  for (auto const& moments : weights_corner) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "backward::corner_swept_centroid["<< moments.entityID <<"]: ";
      std::cout << centroid[0] <<", "<< centroid[1] << ", " << centroid[2];
      std::cout << std::endl;
    }
#endif
    switch(moments.entityID) {
      case corner_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_NEAR(volume, 2 * unit_region_volume, epsilon);
            ASSERT_NEAR(centroid[0], 5.0, epsilon);
            ASSERT_NEAR(centroid[1], 5.0, epsilon);
            ASSERT_NEAR(centroid[2], 5.0, epsilon); break;
          case 2:
            ASSERT_NEAR(volume, unit_region_volume, epsilon);
            ASSERT_NEAR(centroid[0], 4.5, epsilon);
            ASSERT_NEAR(centroid[1], 4.5, epsilon);
            ASSERT_NEAR(centroid[2], 5.5, epsilon); break;
          case 3:
            ASSERT_NEAR(volume, unit_region_volume, epsilon);
            ASSERT_NEAR(centroid[0], 5.5, epsilon);
            ASSERT_NEAR(centroid[1], 4.5, epsilon);
            ASSERT_NEAR(centroid[2], 4.5, epsilon); break;
          case 4:
            ASSERT_NEAR(volume, unit_region_volume, epsilon);
            ASSERT_NEAR(centroid[0], 4.5, epsilon);
            ASSERT_NEAR(centroid[1], 5.5, epsilon);
            ASSERT_NEAR(centroid[2], 4.5, epsilon); break;
          default: FAIL() << "backward::corner: invalid self weights count";
        } break;
      case 17:
        ASSERT_NEAR(volume, unit_region_volume, epsilon);
        ASSERT_NEAR(centroid[0], 3.5, epsilon);
        ASSERT_NEAR(centroid[1], 4.5, epsilon);
        ASSERT_NEAR(centroid[2], 4.5, epsilon); break;
      case 23:
        ASSERT_NEAR(volume, unit_region_volume, epsilon);
        ASSERT_NEAR(centroid[0], 4.5, epsilon);
        ASSERT_NEAR(centroid[1], 3.5, epsilon);
        ASSERT_NEAR(centroid[2], 4.5, epsilon); break;
      case 25:
        ASSERT_NEAR(volume, unit_region_volume, epsilon);
        ASSERT_NEAR(centroid[0], 4.5, epsilon);
        ASSERT_NEAR(centroid[1], 4.5, epsilon);
        ASSERT_NEAR(centroid[2], 3.5, epsilon); break;
      default: FAIL() << "backward::corner: unexpected moment entity index";
    }
  }
}

TEST_F(IntersectSweptOneAxis3D, MomentsCheck) {

  Intersector intersector(source_mesh_wrapper,
                          source_state_wrapper,
                          target_mesh_wrapper,
                          num_tols);

  /* pick internal, boundary and corner source cells.
   * swept region configuration for forward displacement:
   *
   *   source hex        target hex        frontal face
   *                                       swept polyhedron:
   *
   *    7______6         7'......6'
   *    /|    /|          .:    .:
   *   / |   / |         . :   . :
   * 4 __|__5  |      4'...:..5' :         4___4'..5..5'
   * |   3__|__2       :   :..:..2'        |   :  |   :
   * |  /   |  /       :  .   :  .         |   :  |   :
   * | /    | /        : .    : .          |   :  |   :
   * |/_____|/         :......:            |___:..|...:
   * 0      1          0'     1'           0   0' 1   1'
   */
  int const internal_cell = 13;
  int const boundary_cell = 25;
  int const corner_cell   = 26;

  // search for candidate cells and compute moments of intersection
  auto weights_internal = intersector(internal_cell, search(internal_cell));
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));
  auto weights_corner   = intersector(corner_cell, search(corner_cell));

  /* since the target cell was swept along the x-axis only then we should have:
   * - two flat swept regions with zero volume which will be excluded.
   * - one positive swept volume associated with cell 22 (right).
   * - one negative swept volume associated with the source cell itself.
   * - the unsigned volume of the source cell itself.
   * hence the source cell self-contribution is half of its volume.
   */
  double source_volume = source_mesh_wrapper.cell_volume(internal_cell);
  double target_volume = compute_swept_volume(weights_internal);
  double self_contrib  = compute_contribution(internal_cell, weights_internal);
  int nb_self_weights = 0;

  ASSERT_EQ(weights_internal.size(), unsigned(3));
  ASSERT_DOUBLE_EQ(source_volume, target_volume);
  ASSERT_DOUBLE_EQ(self_contrib, unit_region_volume);

  for (auto const& moments : weights_internal) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "one-axis::internal_swept_centroid[" << moments.entityID << "]: ";
      std::cout << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << std::endl;
    }
#endif
    switch(moments.entityID) {
      case internal_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 3.0);
            ASSERT_DOUBLE_EQ(centroid[1], 3.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 2.5);
            ASSERT_DOUBLE_EQ(centroid[1], 3.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          default: FAIL() << "one-axis::internal: invalid self weights count";
        } break;
      case 22:
        ASSERT_DOUBLE_EQ(volume, unit_region_volume);
        ASSERT_DOUBLE_EQ(centroid[0], 4.5);
        ASSERT_DOUBLE_EQ(centroid[1], 3.0);
        ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
      default: FAIL() << "one-axis::internal: unexpected moment entity index";
    }
  }

  /* for the boundary cell case:
   * - one positive swept volume lying outside the source domain and hence excluded.
   * - one negative swept volume associated with the source cell itself.
   * - the unsigned volume of the source cell itself.
   * hence the source cell self-contribution is again half of its volume.
   */
  source_volume = source_mesh_wrapper.cell_volume(boundary_cell);
  target_volume = compute_swept_volume(weights_boundary);
  self_contrib  = compute_contribution(boundary_cell, weights_boundary);
  nb_self_weights = 0;

  ASSERT_EQ(weights_boundary.size(), unsigned(2));
  ASSERT_DOUBLE_EQ(target_volume, 0.5 * source_volume);
  ASSERT_DOUBLE_EQ(self_contrib, unit_region_volume);

  for (auto const& moments : weights_boundary) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "one-axis::boundary_swept_centroid[" << moments.entityID << "]: ";
      std::cout << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << std::endl;
    }
#endif
    switch(moments.entityID) {
      case boundary_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.0);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 4.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 3.0); break;
          default: FAIL() << "one-axis::boundary: invalid self weights count";
        } break;
      default: FAIL() << "one-axis::boundary: unexpected moment entity index";
    }
  }

  /* the corner cell case is no different from the boundary one, since the
   * target cell is swept along the x-axis only. Hence we have:
   * - one positive swept volume lying outside the source domain and hence excluded.
   * - one negative swept volume associated with the source cell itself.
   * - the unsigned volume of the source cell itself.
   * hence the source cell self-contribution is exactly the same as above.
   */
  source_volume = source_mesh_wrapper.cell_volume(corner_cell);
  target_volume = compute_swept_volume(weights_corner);
  self_contrib  = compute_contribution(corner_cell, weights_corner);
  nb_self_weights = 0;

  ASSERT_EQ(weights_corner.size(), weights_boundary.size());
  ASSERT_DOUBLE_EQ(target_volume, 0.5 * source_volume);
  ASSERT_DOUBLE_EQ(self_contrib, unit_region_volume);

  for (auto const& moments : weights_corner) {
    auto const volume = std::abs(moments.weights[0]);
    auto const centroid = deduce_centroid(moments);
#if DEBUG
    if (verbose) {
      std::cout << "one-axis::corner_swept_centroid[" << moments.entityID << "]: ";
      std::cout << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << std::endl;
    }
#endif
    switch(moments.entityID) {
      case corner_cell:
        switch (++nb_self_weights) {
          case 1:
            ASSERT_DOUBLE_EQ(volume, 2 * unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 5.0);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 5.0); break;
          case 2:
            ASSERT_DOUBLE_EQ(volume, unit_region_volume);
            ASSERT_DOUBLE_EQ(centroid[0], 4.5);
            ASSERT_DOUBLE_EQ(centroid[1], 5.0);
            ASSERT_DOUBLE_EQ(centroid[2], 5.0); break;
          default: FAIL() << "one-axis::corner: invalid self weights count";
        } break;
      default: FAIL() << "one-axis::corner: unexpected moment entity index";
    }
  }
}
