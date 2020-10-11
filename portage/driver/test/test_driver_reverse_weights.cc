/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "portage/support/portage.h"
#if defined(PORTAGE_HAS_TANGRAM) && defined(WONTON_ENABLE_MPI)
// system
#include <iostream>
#include "gtest/gtest.h"

// wonton
#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/mesh/flat/flat_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/state/flat/flat_state_mm_wrapper.h"

// portage
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/driver/coredriver.h"
#include "portage/distributed/mpi_bounding_boxes.h"

// jali
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"

// tangram
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
    return std::any_of(list.begin(), list.end(),
                       [source](auto const& weight) { return weight.entityID == source; });
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
                                       Tangram::MOF, Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>;

  MPI_Comm comm = MPI_COMM_WORLD;
  auto source_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto target_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 7);
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  int const num_source_cells = source_mesh_wrapper.num_entities(Wonton::CELL, Wonton::ALL);

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

  std::string const materials[] = {"mat0", "mat1", "mat2"};
  Wonton::Point<2> const lower_bound[] = {{0.0, 0.0}, {0.5, 0.0}, {0.5, 0.5}};
  Wonton::Point<2> const upper_bound[] = {{0.5, 1.0}, {1.0, 0.5}, {1.0, 1.0}};

  int const num_mats = 3;
  std::vector<int> material_cells[num_mats];
  std::vector<double> volume_fractions[num_mats];
  std::vector<Wonton::Point<2>> centroids[num_mats];
  std::vector<std::vector<int>> index_lookup(num_mats);

  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  for (int c = 0; c < num_source_cells; c++) {
    std::vector<Wonton::Point<2>> coords;
    Wonton::Point<2> box[2];
    std::vector<double> moments;

    double cellvol = source_mesh_wrapper.cell_volume(c);
    source_mesh_wrapper.cell_get_coordinates(c, &coords);
    BOX_INTERSECT::bounding_box<2>(coords, box, box+1);

    for (int m = 0; m < num_mats; m++) {
      bool intersected = BOX_INTERSECT::intersect_boxes<2>(lower_bound[m], upper_bound[m], box[0], box[1], &moments);
      if (intersected and moments[0] > 1.E-6) {  // non-trivial intersection
        material_cells[m].emplace_back(c);
        volume_fractions[m].emplace_back(moments[0] / cellvol);
        centroids[m].emplace_back(moments[1] / moments[0], moments[2] / moments[0]);
      }
    }
  }

  //-------------------------------------------------------------------
  // Now add the material and material cells to the source state
  //-------------------------------------------------------------------
  for (int m = 0; m < num_mats; m++) {
    source_state_wrapper.add_material(materials[m], material_cells[m]);
    source_state_wrapper.mat_add_celldata("mat_volfracs", m, volume_fractions[m].data());
    source_state_wrapper.mat_add_celldata("mat_centroids", m, centroids[m].data());
    source_state_wrapper.mat_add_celldata<double>("density", m, 0.0);
    target_state_wrapper.add_material(materials[m], {});
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
    return std::any_of(list.begin(), list.end(),
                       [s](auto const& weight) { return weight.entityID == s; });
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

TEST(ReverseWeights, Redistributed) {

  using Remapper = Portage::CoreDriver<2, Wonton::CELL,
                                       Wonton::Flat_Mesh_Wrapper<>,
                                       Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>>,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Tangram::MOF, Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>;

  MPI_Comm comm = MPI_COMM_WORLD;
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.included_entities(Jali::Entity_kind::ALL_KIND);
  mesh_factory.partitioner(Jali::Partitioner_type::BLOCK);

  auto source_mesh  = mesh_factory(0.0, 0.0, 1.0, 1.0, 10, 10);
  auto target_mesh  = mesh_factory(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);
  Wonton::MPIExecutor_type executor(comm);

  int const num_source_cells = source_mesh_wrapper.num_entities(Wonton::CELL, Wonton::ALL);

  //    0,1            0.5,1.          1,1
  //     *------------------------------*
  //     |            0, mat0           |
  //     |                              |
  //     |    *....................*    |
  //     |    :       1, mat1      :    |
  //     |    :     *.........*    :    |
  //     |    :     :         :    :    |
  //     |    :     : 2, mat2 :    :    |
  //     |    :     :         :    :    | 1,0.5
  //     |    :     :         :    :    |
  //     |    :     *.........*    :    |
  //     |    :                    :    |
  //     |    *....................*    |
  //     |                              |
  //     |                              |
  //     *------------------------------*
  //    0,0            0.5,0           1,0

  std::string const materials[] = {"mat0", "mat1", "mat2"};

  std::vector<Wonton::Point<2>> const lower_bound[] = {
    {{0.0, 0.0}, {0.0, 0.2}, {0.8, 0.2}, {0.0, 0.8}},
    {{0.2, 0.2}, {0.2, 0.4}, {0.6, 0.4}, {0.2, 0.6}},
    {{0.4, 0.4}}
  };

  std::vector<Wonton::Point<2>> const upper_bound[] = {
    {{1.0, 0.2}, {0.2, 0.8}, {1.0, 0.8}, {1.0, 1.0}},
    {{0.8, 0.4}, {0.4, 0.6}, {0.8, 0.6}, {0.8, 0.8}},
    {{0.6, 0.6}},
  };

  int const num_boxes[] = {4, 4, 1};

  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  double const volume_tolerance = 1.E-6;

  int const num_mats = 3;
  std::vector<int> material_cells[num_mats];
  std::vector<double> volume_fractions[num_mats];
  std::vector<Wonton::Point<2>> centroids[num_mats];
  std::vector<std::vector<int>> index_lookup(num_mats);

  for (int c = 0; c < num_source_cells; c++) {
    std::vector<Wonton::Point<2>> coords;
    Wonton::Point<2> box[2];
    std::vector<double> moments;

    double const cellvol = source_mesh_wrapper.cell_volume(c);
    source_mesh_wrapper.cell_get_coordinates(c, &coords);
    BOX_INTERSECT::bounding_box<2>(coords, box, box+1);

    for (int m = 0; m < num_mats; m++) {
      double matvol = 0.;
      Wonton::Point<2> matcen(0.,0.);

      for (int nb = 0; nb < num_boxes[m]; nb++) {
        bool intersected = BOX_INTERSECT::intersect_boxes<2>(lower_bound[m][nb], upper_bound[m][nb], box[0], box[1], &moments);
        if (intersected and moments[0] > volume_tolerance) {  // non-trivial intersection
          matvol += moments[0];
          matcen[0] += moments[1];
          matcen[1] += moments[2];
        }
      } // num_boxes
      if (matvol > volume_tolerance) {
        material_cells[m].emplace_back(c);
        volume_fractions[m].emplace_back(matvol / cellvol);
        centroids[m].emplace_back(matcen[0] / matvol, matcen[1] / matvol);
      }
    }
  }

  //---------------------------------------------------------
  // add the material and material cells to the source state
  //---------------------------------------------------------
  for (int m = 0; m < num_mats; m++) {
    source_state_wrapper.add_material(materials[m], material_cells[m]);
    source_state_wrapper.mat_add_celldata("mat_volfracs", m, volume_fractions[m].data());
    source_state_wrapper.mat_add_celldata("mat_centroids", m, centroids[m].data());
    source_state_wrapper.mat_add_celldata<double>("density", m, 0.0);
    target_state_wrapper.add_material(materials[m], {});
    target_state_wrapper.mat_add_celldata<double>("density", m, 0.0);
  }

  //------------------------------------------------------------------------
  // redistribute source entities to ensure perfect overlap with target ones
  //------------------------------------------------------------------------
  Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
  source_mesh_flat.initialize(source_mesh_wrapper);
  Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);
  source_state_flat.initialize(source_state_wrapper, {"density"});

  Portage::MPI_Bounding_Boxes redistributor(&executor);
  redistributor.distribute(source_mesh_flat, source_state_flat,
                           target_mesh_wrapper, target_state_wrapper);

  //-------------------------------
  // compute interpolation weights
  //-------------------------------
  Remapper remapper(source_mesh_flat, source_state_flat,
                    target_mesh_wrapper, target_state_wrapper, &executor);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto forward_weights = remapper.intersect_materials<Portage::IntersectR2D>(candidates);
  auto reverse_weights = remapper.deduce_reverse_material_weights(forward_weights, index_lookup);

  // check if the given source cell is contained in
  // the weight list of the target cell for forward remap.
  auto weights_matches = [&](int m, int s, int t) -> bool {
    std::vector<Wonton::Weights_t> const& list = forward_weights[m][t];
    return std::any_of(list.begin(), list.end(),
                       [s](auto const& weight) { return weight.entityID == s; });
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