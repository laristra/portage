/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "gtest/gtest.h"

#ifdef PORTAGE_ENABLE_MPI
  #include "mpi.h"
#endif

#include <numeric>

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_swept_face.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"


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
  using Intersector = Portage::IntersectSweptFace3D<Wonton::Entity_kind::CELL,
                                                    Wonton::Jali_Mesh_Wrapper,
                                                    Wonton::Jali_State_Wrapper,
                                                    Wonton::Jali_Mesh_Wrapper>;
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
    return (std::abs(value) < num_tols.min_relative_volume ? 0.0 : value);
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
   * we expect to have a constant volume for each swept region.
   */
  double const unit_region_volume = 4.0;

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
//  int const nb_cells = source_mesh_wrapper.num_entities(Wonton::Entity_kind::CELL,
//                                                        Wonton::Entity_type::ALL);
//  for (int i = 0; i < nb_cells; ++i) {
//    auto centroid = source_mesh->cell_centroid(i);
//    std::cout << "cell["<< i <<"]: "<< centroid << std::endl;
//  }
////
////  source_mesh->write_to_exodus_file("source_mesh.exo");

  int const internal_cell = 13;
  int const boundary_cell = 25;
  int const corner_cell   = 26;
  int const nb_hex_faces  = 6;

  auto centroid_internal = source_mesh->cell_centroid(internal_cell);
  auto centroid_boundary = source_mesh->cell_centroid(boundary_cell);
  auto centroid_corner   = source_mesh->cell_centroid(corner_cell);

  // search for candidate cells and compute moments of intersection
  auto weights_internal = intersector(internal_cell, search(internal_cell));
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));
  auto weights_corner   = intersector(corner_cell, search(corner_cell));

  /* for interior hexahedral cell, we have:
   * - three positive swept volumes attached to cells 14:top|16:right|22:forward (V=-12).
   * - three negative swept volume attached to the source cell itself (V=12).
   * - the unsigned volume of the source cell itself (V=8).
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
   *
   */
  source_volume = source_mesh_wrapper.cell_volume(boundary_cell);
  target_volume = compute_swept_volume(weights_boundary);
  self_contrib  = compute_contribution(boundary_cell, weights_boundary);
  nb_self_weights = 0;

  ASSERT_EQ(weights_boundary.size(), unsigned(5));
  std::cout << "target_volume: " << target_volume << std::endl;
  std::cout << "self_contrib: " << self_contrib << std::endl;

}