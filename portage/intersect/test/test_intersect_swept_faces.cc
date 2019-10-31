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
 * by a displacement vector v=(1,1) and we aim to check the volume
 * of each swept face for a given source cell.
 * Notice that mesh topology is unchanged, hence source and target
 * cells and faces indices are kept.
 *
 *    .............
 *   _:___:___:_  :     source mesh: plain
 *  | :.|.:.|.:.|.:     target mesh: dotted
 *  |_:_|_:_|_:_| :
 *  | :.|.:.|.:.|.:     cell indices ordering:
 *  |_:_|_:_|_:_| :     6 7 8
 *  | :.|.:.|.:.|.:     3 4 5
 *  |___|___|___|       0 1 2
 *  0   2   4   6
 */
class IntersectSweptTest : public testing::Test {

protected:
  // useful shortcuts
  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
    Wonton::Jali_Mesh_Wrapper,
    Wonton::Jali_State_Wrapper>;

  /**
   * @brief Setup each test-case.
   *
   * It initializes both source and target meshes and states,
   * then computes and assigns a density field on source mesh.
   */
  IntersectSweptTest()
    : source_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 6.0, 6.0, 3, 3)),
      target_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(1.0, 1.0, 7.0, 7.0, 3, 3)),
      source_state(Jali::State::create(source_mesh)),
      target_state(Jali::State::create(target_mesh)),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state_wrapper(*source_state),
      target_state_wrapper(*target_state)
  {
    num_tols.use_default();
    // save meshes sizes
    nb_source_cells = source_mesh_wrapper.num_owned_cells();
    nb_target_cells = target_mesh_wrapper.num_owned_cells();

  }

  /**
   * @brief Retrieve all cells incident to the faces of a given cell.
   *
   * @param cell: the current cell.
   * @return a list of cell neighbors incident to the faces of that cell.
   */
  std::vector<int> get_stencil(int current) {

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

protected:
  // useful constants and aliases
  static constexpr double upper_bound = std::numeric_limits<double>::max();
  static constexpr double lower_bound = std::numeric_limits<double>::min();
  static constexpr double epsilon = 1.E-10;
  static constexpr auto CELL = Wonton::Entity_kind::CELL;

  int nb_source_cells = 0;
  int nb_target_cells = 0;

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
};


TEST_F(IntersectSweptTest, SweptAreaCheck) {
  // todo
//  // retrieve the central cell index
//  std::for_each(source_mesh_wrapper.begin(CELL),
//                   source_mesh_wrapper.end(CELL),
//                   [&](int current) {
//    auto centroid = source_mesh->cell_centroid(current);
//    std::cout << "cell["<< current <<"]: "<< centroid << std::endl;
//  });

  using Intersector = Portage::IntersectSwept<2, Wonton::Entity_kind::CELL,
                                                 Wonton::Jali_Mesh_Wrapper,
                                                 Wonton::Jali_State_Wrapper,
                                                 Wonton::Jali_Mesh_Wrapper>;

  // pick source cell and retrieve its stencil
  int const source_id = 4;
  auto const stencil = get_stencil(source_id);

  ASSERT_EQ(stencil.size(), unsigned(5));
  ASSERT_EQ(stencil[1], 3);
  ASSERT_EQ(stencil[2], 7);
  ASSERT_EQ(stencil[3], 5);
  ASSERT_EQ(stencil[4], 1);

  Intersector swept_intersector(source_mesh_wrapper,
                                source_state_wrapper,
                                target_mesh_wrapper,
                                num_tols);

  auto const source_weights = swept_intersector(source_id, stencil);
  for (auto&& entry : source_weights) {
    std::cout << "weight.entity_id: " << entry.entityID << " [";
    for (auto&& weight : entry.weights)
      std::cout << weight << ", ";
    std::cout << "]" << std::endl;
  }
}


TEST_F(IntersectSweptTest, RemapCheck) {
  // todo
}