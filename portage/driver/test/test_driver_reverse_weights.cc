/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "portage/support/portage.h"
#ifdef PORTAGE_HAS_TANGRAM
#include <iostream>
#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/driver/coredriver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"

#include "tangram/driver/driver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/intersect/split_r2d.h"

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

  int const num_source_cells = source_mesh_wrapper.num_entities(Wonton::CELL, Wonton::ALL);

#ifdef DEBUG
  int const num_target_cells = target_mesh_wrapper.num_entities(Wonton::CELL,
                                                                Wonton::PARALLEL_OWNED);
  for (int t = 0; t < num_target_cells; ++t) {
    std::vector<Wonton::Weights_t> const& list = forward_weights[t];
    std::cout << "forward_weight["<< t <<"]: (";
    int const num_weights = list.size();
    for (int i = 0; i < num_weights; ++i) {
      auto const& weight = list[i];
      std::cout << weight.entityID;
      if (i < num_weights - 1)
        std::cout <<", ";
    }
    std::cout << ")" << std::endl;
  }
  std::cout << "--------------" << std::endl;

  for (int s = 0; s < num_source_cells; ++s) {
    std::vector<Wonton::Weights_t> const& list = reverse_weights[s];
    std::cout << "reverse_weight["<< s <<"]: (";
    int const num_weights = list.size();
    for (int i = 0; i < num_weights; ++i) {
      auto const& weight = list[i];
      std::cout << weight.entityID;
      if (i < num_weights - 1)
        std::cout <<", ";
    }
    std::cout << ")" << std::endl;
  }
#endif

  // check if the given source cell is contained in
  // the weight list of the target cell for forward remap.
  auto weights_matches = [&](int source, int target) -> bool {
    std::vector<Wonton::Weights_t> const& list = forward_weights[target];
    for (auto&& weight : list) {
      if (weight.entityID == source) { return true; }
    }
    return false;
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

TEST(ReverseWeights, MultiMat) {

  using Remapper = Portage::CoreDriver<2, Wonton::CELL,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>;

  MPI_Comm comm = MPI_COMM_WORLD;
  auto source_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto target_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  int num_source_cells = source_mesh_wrapper.num_entities(Wonton::CELL, Wonton::ALL);

  // The material geometry in the overall domain will look like this
  // and we will put down a rectangular mesh that has multiple cells
  // in each direction on this domain so that we get some pure and
  // some mixed cells
  //
  // Note that only MOF type algorithms or VOF algorithms with the
  // material ordering 0,1,2 will get the T-junction geometry right
  //
  //    0,1           0.5,1         1,1
  //     *-------------:------------*
  //     |             :            |
  //     |             :        2   |
  //     |             :     mat2   |
  //     |             :            |
  //     |             :            |
  //     |     0       +............|1,0.5
  //     |   mat0      :            |
  //     |             :            |
  //     |             :     mat1   |
  //     |             :        1   |
  //     |             :            |
  //     *-------------:------------*
  //    0,0           0.5,0         1,0

  int const num_mats = 3;
  std::string matnames[] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  Wonton::Point<2> matlo[] = {{0.0, 0.0}, {0.5, 0.0}, {0.5, 0.5}};
  Wonton::Point<2> mathi[] = {{0.5, 1.0}, {1.0, 0.5}, {1.0, 1.0}};

  std::vector<int> matcells_src[num_mats];
  std::vector<double> matvf_src[num_mats];
  std::vector<Wonton::Point<2>> matcen_src[num_mats];
  std::vector<std::vector<int>> index_lookup(num_mats);

  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  for (int c = 0; c < num_source_cells; c++) {
    std::vector<Wonton::Point<2>> coords;
    source_mesh_wrapper.cell_get_coordinates(c, &coords);

    double cellvol = source_mesh_wrapper.cell_volume(c);

    Wonton::Point<2> box[2];
    BOX_INTERSECT::bounding_box<2>(coords, box, box+1);

    std::vector<double> moments;
    for (int m = 0; m < num_mats; m++) {
      bool intersected = BOX_INTERSECT::intersect_boxes<2>(matlo[m], mathi[m],box[0], box[1], &moments);
      if (intersected and moments[0] > 1.E-6) {  // non-trivial intersection
        matcells_src[m].emplace_back(c);
        matvf_src[m].emplace_back(moments[0] / cellvol);
        matcen_src[m].emplace_back(moments[1] / moments[0], moments[2] / moments[0]);
      }
    }
  }

  //-------------------------------------------------------------------
  // Now add the material and material cells to the source state
  //-------------------------------------------------------------------
  for (int m = 0; m < num_mats; m++) {
    source_state_wrapper.add_material(matnames[m], matcells_src[m]);
    source_state_wrapper.mat_add_celldata("mat_volfracs", m, matvf_src[m].data());
    source_state_wrapper.mat_add_celldata("mat_centroids", m, matcen_src[m].data());
    source_state_wrapper.mat_add_celldata<double>("density", m, 0.0);
    target_state_wrapper.add_material(matnames[m], {});
    target_state_wrapper.mat_add_celldata<double>("density", m, 0.0);
  }

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto forward_weights = remapper.intersect_materials<Portage::IntersectR2D>(candidates);
  auto reverse_weights = remapper.deduce_reverse_material_weights(forward_weights, index_lookup);

#ifdef DEBUG
  for (int m = 0; m < 3; ++m) {
    int num_material_cells = forward_weights[m].size();
    for (int t = 0; t < num_material_cells; ++t) {
      std::vector<Wonton::Weights_t> const& list = forward_weights[m][t];
      std::cout << "mat: "<< m << ", forward_weight["<< t <<"]: (";
      int const num_weights = list.size();
      for (int i = 0; i < num_weights; ++i) {
        auto const& weight = list[i];
        std::cout << weight.entityID;
        if (i < num_weights - 1)
          std::cout <<", ";
      }
      std::cout << ")" << std::endl;
    }
    std::cout << "--------------" << std::endl;

    num_material_cells = reverse_weights[m].size();
    for (int s = 0; s < num_material_cells; ++s) {
      std::vector<Wonton::Weights_t> const& list = reverse_weights[m][s];
      std::cout << "mat: " << m << ", reverse_weight[" << index_lookup[m][s] << "]: (";
      int const num_weights = list.size();
      for (int i = 0; i < num_weights; ++i) {
        auto const& weight = list[i];
        std::cout << weight.entityID;
        if (i < num_weights - 1)
          std::cout <<", ";
      }
      std::cout << ")" << std::endl;
    }
    std::cout << "=====================" << std::endl;
  }
#endif

  // check if the given source cell is contained in
  // the weight list of the target cell for forward remap.
  auto weights_matches = [&](int m, int s, int t) -> bool {
    std::vector<Wonton::Weights_t> const& list = forward_weights[m][t];
    for (auto&& weight : list) {
      if (weight.entityID == s) { return true; }
    }
    return false;
  };

  // check that the reverse weight sparse matrix is a
  // perfect transposition of the forward weight matrix.
  for (int m = 0; m < num_mats; ++m) {
    int const num_source_material_cells = reverse_weights[m].size();
    for (int i = 0; i < num_source_material_cells; ++i) {
      std::vector<Wonton::Weights_t> const& list = reverse_weights[m][i];
      for (auto const& weight : list) {
        int const& s = index_lookup[m][i];
        int const& t = weight.entityID;
        ASSERT_TRUE(weights_matches(m, s, t));
      }
    }
  }
}
#endif