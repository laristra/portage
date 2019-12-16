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
  static double compute_swept_volume(std::vector<Wonton::Weights_t> const& moments) {
    return std::accumulate(moments.begin(), moments.end(), 0.0,
                           [](double previous, auto const& moment) {
                             return previous + moment.weights[0];
                           });
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
  double const unit_region_volume = 2.0;

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

  auto centroid_internal = source_mesh->cell_centroid(internal_cell);
  auto centroid_boundary = source_mesh->cell_centroid(boundary_cell);
  auto centroid_corner   = source_mesh->cell_centroid(corner_cell);

  // search for candidate cells and compute moments of intersection
  std::cout << "= internal =" << std::endl;
  auto weights_internal = intersector(internal_cell, search(internal_cell));

  std::cout << "= boundary =" << std::endl;
  auto weights_boundary = intersector(boundary_cell, search(boundary_cell));

  std::cout << "= corner =" << std::endl;
  auto weights_corner = intersector(corner_cell, search(corner_cell));

  std::cout << "internal list size: "<< weights_internal.size() - 1 << std::endl;
  std::cout << "boundary list size: "<< weights_boundary.size() -1 << std::endl;
  std::cout << "corner list size: "<< weights_corner.size() -1 << std::endl;
}