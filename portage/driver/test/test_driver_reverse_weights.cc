/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include <iostream>
#include <memory>
#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/driver/coredriver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

TEST(ReverseWeights, SingleMat) {

  using Remapper = Portage::CoreDriver<2, Wonton::CELL,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Wonton::Jali_State_Wrapper>;

  MPI_Comm comm = MPI_COMM_WORLD;
  auto source_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto target_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto forward_weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);
  auto reverse_weights = remapper.deduce_reverse_weights(forward_weights);

#ifdef DEBUG
  int const num_target_cells = target_mesh_wrapper.num_entities(Wonton::CELL,
                                                                Wonton::PARALLEL_OWNED);
  for (int t = 0; t < num_target_cells; ++t) {
    std::vector<Wonton::Weights_t> const& list = forward_weights[t];
    std::cout << "forward_weight[target: "<< t <<"]: (";
    for (auto const& weight : list) {
      std::cout << weight.entityID <<", ";
    }
    std::cout << std::endl;
  }
#endif

  int const num_source_cells = source_mesh_wrapper.num_entities(Wonton::CELL, Wonton::ALL);

#ifdef DEBUG
  for (int s = 0; s < num_source_cells; ++s) {
    std::vector<Wonton::Weights_t> const& list = reverse_weights[s];
    std::cout << "reverse_weight[source: "<< s <<"]: (";
    for (auto const& weight : list) {
      std::cout << weight.entityID <<", ";
    }
    std::cout << std::endl;
  }
#endif

  // check if the given source cell is contained in
  // the weight list of the target cell for forward remap.
  auto weights_matches = [&](int source, int target) -> bool {
    std::vector<Wonton::Weights_t> const& list = forward_weights[target];
    bool found = false;
    for (auto&& weight : list) {
      if (weight.entityID == source) {
        found = true;
        break;
      }
    }
    return found;
  };

  // check that the reverse weight sparse matrix is a
  // perfect transposition of the forward weight matrix.
  for (int s = 0; s < num_source_cells; ++s) {
    std::vector<Wonton::Weights_t> const& list = reverse_weights[s];
    for (auto const& weight : list) {
      int const& t = weight.entityID;
      ASSERT_TRUE(weights_matches(s, t));
    }
  }
}