/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */


#include "gtest/gtest.h"

// wonton
#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/mesh/flat/flat_mesh_wrapper.h"
#include "wonton/state/flat/flat_state_mm_wrapper.h"
#include "wonton/distributed/mpi_ghost_manager.h"

// portage
#include "portage/support/portage.h"
#include "portage/driver/coredriver.h"
#include "portage/distributed/mpi_bounding_boxes.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/interpolate/interpolate_2nd_order.h"

// jali
#include "Mesh.hh"
#include "MeshFactory.hh"

#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/reconstruct/VOF.h"
#include "tangram/intersect/split_rNd.h"
#endif


TEST(GhostManager, UpdateValues) {

  using GhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                Wonton::Jali_State_Wrapper,
                                                Wonton::CELL>;

  using Remapper = Portage::CoreDriver<2, Wonton::CELL,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Wonton::Jali_State_Wrapper>;

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);
  Wonton::MPIExecutor_type executor(comm);

  auto jali_source_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto jali_target_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto jali_source_state = Jali::State::create(jali_source_mesh);
  auto jali_target_state = Jali::State::create(jali_target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh(*jali_source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh(*jali_target_mesh);
  Wonton::Jali_State_Wrapper source_state(*jali_source_state);
  Wonton::Jali_State_Wrapper target_state(*jali_target_state);

  int const num_source_cells = source_mesh.num_entities(Wonton::CELL, Wonton::ALL);

  double source_field[num_source_cells];
  for (int i = 0; i < num_source_cells; ++i) {
    Wonton::Point<2> centroid;
    source_mesh.cell_centroid(i, &centroid);
    source_field[i] = centroid[0] + 2 * centroid[1];
  }

  source_state.mesh_add_data<double>(Wonton::CELL, "temperature", source_field);
  target_state.mesh_add_data<double>(Wonton::CELL, "temperature", 0.0);

  // remap
  Remapper remapper(source_mesh, source_state, target_mesh, target_state, &executor);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights    = remapper.intersect_meshes<Portage::IntersectRnD>(candidates);
  auto gradients  = remapper.compute_source_gradient("temperature");

  remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "temperature", "temperature", weights, &gradients
  );

  // fill ghost values on target mesh
  GhostManager ghost_manager(target_mesh, target_state, comm);
  ghost_manager.update_ghost_values("temperature", true);

  // gather all sent/received values from all ranks to the master
  std::vector<std::vector<double>> send[num_ranks]; /* values of owned cells that are ghost for some rank */
  std::vector<std::vector<double>> take[num_ranks]; /* values of ghost cells that are owned by some rank */
  std::vector<int> num_sent[num_ranks];             /* number of owned cells sent by each rank */
  std::vector<int> num_take[num_ranks];             /* number of ghost cells receive from each rank */
  std::vector<MPI_Request> requests;                /* list of asynchronous MPI requests */

  send[rank] = ghost_manager.owned_values();
  take[rank] = ghost_manager.ghost_values();

  for (int i = 0; i < num_ranks; ++i) {
    num_sent[rank].emplace_back(send[rank][i].size());
    num_take[rank].emplace_back(take[rank][i].size());
  }

  if (rank > 0) {
    MPI_Request request;
    MPI_Isend(num_sent[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
    requests.emplace_back(request);
    MPI_Isend(num_take[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
    requests.emplace_back(request);
  } else {
    for (int i = 1; i < num_ranks; ++i) {
      MPI_Request request;
      num_sent[i].resize(num_ranks, 0);
      MPI_Irecv(num_sent[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
      requests.emplace_back(request);
      num_take[i].resize(num_ranks, 0);
      MPI_Irecv(num_take[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
      requests.emplace_back(request);
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
  requests.clear();

  if (rank > 0) {
    int const offset = num_ranks * num_ranks;
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        int tag = num_ranks * (rank - 1) + i;
        MPI_Request request;
        MPI_Isend(send[rank][i].data(), send[rank][i].size(), MPI_DOUBLE, 0, tag, comm, &request);
        requests.emplace_back(request);
        MPI_Isend(take[rank][i].data(), take[rank][i].size(), MPI_DOUBLE, 0, tag + offset, comm, &request);
        requests.emplace_back(request);
      }
    }
  } else {
    int const offset = num_ranks * num_ranks;
    for (int i = 1; i < num_ranks; ++i) {
      send[i].resize(num_ranks);
      take[i].resize(num_ranks);
      for (int j = 0; j < num_ranks; ++j) {
        if (i != j) {
          int tag = num_ranks * (i - 1) + j;
          MPI_Request request;
          send[i][j].resize(num_sent[i][j]);
          MPI_Irecv(send[i][j].data(), num_sent[i][j], MPI_DOUBLE, i, tag, comm, &request);
          requests.emplace_back(request);
          take[i][j].resize(num_take[i][j]);
          MPI_Irecv(take[i][j].data(), num_take[i][j], MPI_DOUBLE, i, tag + offset, comm, &request);
          requests.emplace_back(request);
        }
      }
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

  if (rank == 0) {
    for (int i = 0; i < num_ranks; ++i) {
      for (int j = 0; j < num_ranks; ++j) {
        if (i == j) { continue; }
        int const num_owned = send[i][j].size();
        int const num_ghost = take[j][i].size();
        ASSERT_EQ(num_owned, num_ghost);
        for (int k = 0; k < num_owned; ++k) {
          ASSERT_DOUBLE_EQ(send[i][j][k], take[j][i][k]);
        }
      }
    }
  }

  MPI_Barrier(comm);
}

#ifdef PORTAGE_HAS_TANGRAM
TEST(GhostManager, MultiMat) {

  using GhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                Wonton::Jali_State_Wrapper,
                                                Wonton::CELL>;

  using Remapper = Portage::CoreDriver<2, Wonton::CELL,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Tangram::MOF,
                                       Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>;

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);
  Wonton::MPIExecutor_type executor(comm);

  auto source_jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto target_jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto source_jali_state = Jali::State::create(source_jali_mesh);
  auto target_jali_state = Jali::State::create(target_jali_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh(*source_jali_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh(*target_jali_mesh);
  Wonton::Jali_State_Wrapper source_state(*source_jali_state);
  Wonton::Jali_State_Wrapper target_state(*target_jali_state);

  // - step 1: prepare source data
  int num_source_cells = source_mesh.num_entities(Wonton::CELL, Wonton::ALL);

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
  std::string const materials[num_mats] = {"mat0", "mat1", "mat2"};
  Wonton::Point<2> const lower_bound[num_mats] = {{0.0, 0.0}, {0.5, 0.0}, {0.5, 0.5}};
  Wonton::Point<2> const upper_bound[num_mats] = {{0.5, 1.0}, {1.0, 0.5}, {1.0, 1.0}};

  std::vector<int> cells[num_mats];
  std::vector<double> volume_fractions[num_mats];
  std::vector<Wonton::Point<2>> centroids[num_mats];
  std::vector<double> densities[num_mats];

  for (int c = 0; c < num_source_cells; c++) {
    std::vector<Wonton::Point<2>> coords;
    source_mesh.cell_get_coordinates(c, &coords);

    double cellvol = source_mesh.cell_volume(c);
    double const min_vol = 1.E-6;

    Wonton::Point<2> box[2];
    BOX_INTERSECT::bounding_box<2>(coords, box, box+1);

    std::vector<double> moments;
    for (int m = 0; m < num_mats; m++) {
      bool intersected = BOX_INTERSECT::intersect_boxes<2>(lower_bound[m], upper_bound[m], box[0], box[1], &moments);
      if (intersected and moments[0] > min_vol) {  // non-trivial intersection
        cells[m].emplace_back(c);
        volume_fractions[m].emplace_back(moments[0] / cellvol);
        Wonton::Point<2> p(moments[1] / moments[0], moments[2] / moments[0]);
        centroids[m].emplace_back(p);
        densities[m].emplace_back((m + 1) * (m + 1) * (p[0] + p[1]));
      }
    }
  }

  for (int m = 0; m < num_mats; m++) {
    source_state.add_material(materials[m], cells[m]);
    source_state.mat_add_celldata("mat_volfracs", m, volume_fractions[m].data());
    source_state.mat_add_celldata("mat_centroids", m, centroids[m].data());
    source_state.mat_add_celldata("density", m, densities[m].data());
    target_state.add_material(materials[m], {});
    target_state.mat_add_celldata<double>("density", m, 0.0);
  }

  // - step 2: remap material fields
  Remapper remapper(source_mesh, source_state, target_mesh, target_state, &executor);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights    = remapper.intersect_materials<Portage::IntersectRnD>(candidates);
  auto gradients  = remapper.compute_source_material_gradient("density");
  remapper.interpolate_mat_var<double, Portage::Interpolate_2ndOrder>(
    "density", "density", weights, &gradients
  );

  // fill ghost values on target mesh
  GhostManager ghost_manager(target_mesh, target_state, comm);
  ghost_manager.update_ghost_values("density", true);

  for (int m = 1; m <= num_mats; ++m) {

    std::vector<std::vector<double>> send[num_ranks]; /* values of owned cells that are ghost for some rank */
    std::vector<std::vector<double>> take[num_ranks]; /* values of ghost cells that are owned by some rank */
    std::vector<int> num_sent[num_ranks];             /* number of owned cells sent by each rank */
    std::vector<int> num_take[num_ranks];             /* number of ghost cells receive from each rank */
    std::vector<MPI_Request> requests;                /* list of asynchronous MPI requests */

    send[rank] = ghost_manager.owned_values(m);
    take[rank] = ghost_manager.ghost_values(m);

    for (int i = 0; i < num_ranks; ++i) {
      num_sent[rank].emplace_back(send[rank][i].size()); // (!) may be empty
      num_take[rank].emplace_back(take[rank][i].size());
    }

    if (rank > 0) {
      MPI_Request request;
      MPI_Isend(num_sent[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
      MPI_Isend(num_take[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
    } else {
      for (int i = 1; i < num_ranks; ++i) {
        MPI_Request request;
        num_sent[i].resize(num_ranks, 0);
        MPI_Irecv(num_sent[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
        num_take[i].resize(num_ranks, 0);
        MPI_Irecv(num_take[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    requests.clear();

    if (rank > 0) {
      int const offset = num_ranks * num_ranks;
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          MPI_Request request;
          int tag = num_ranks * (rank - 1) + i;
          // nb: can be empty for this material
          if (not send[rank][i].empty()) {
            MPI_Isend(send[rank][i].data(), send[rank][i].size(), MPI_DOUBLE, 0, tag, comm, &request);
            requests.emplace_back(request);
          }
          if (not take[rank][i].empty()) {
            MPI_Isend(take[rank][i].data(), take[rank][i].size(), MPI_DOUBLE, 0, tag + offset, comm, &request);
            requests.emplace_back(request);
          }
        }
      }
    } else {
      int const offset = num_ranks * num_ranks;
      for (int i = 1; i < num_ranks; ++i) {
        send[i].resize(num_ranks);
        take[i].resize(num_ranks);
        for (int j = 0; j < num_ranks; ++j) {
          if (i != j) {
            MPI_Request request;
            int tag = num_ranks * (i - 1) + j;
            // nb: can be empty for this material
            if (num_sent[i][j]) {
              send[i][j].resize(num_sent[i][j]);
              MPI_Irecv(send[i][j].data(), num_sent[i][j], MPI_DOUBLE, i, tag, comm, &request);
              requests.emplace_back(request);
            }
            if (num_take[i][j]) {
              take[i][j].resize(num_take[i][j]);
              MPI_Irecv(take[i][j].data(), num_take[i][j], MPI_DOUBLE, i, tag + offset, comm, &request);
              requests.emplace_back(request);
            }
          }
        }
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

    if (rank == 0) {
      for (int i = 0; i < num_ranks; ++i) {
        for (int j = 0; j < num_ranks; ++j) {
          if (i == j) { continue; }
          int const num_owned = send[i][j].size();
          int const num_ghost = take[j][i].size();
          ASSERT_EQ(num_owned, num_ghost);
          for (int k = 0; k < num_owned; ++k) {
            ASSERT_DOUBLE_EQ(send[i][j][k], take[j][i][k]);
          }
        }
      }
    }

    MPI_Barrier(comm);
  } // end for each material
}


TEST(GhostManager, MultiMatRedistributed) {
  
  using SourceFlatGhostManager =
      Wonton::MPI_GhostManager<Wonton::Flat_Mesh_Wrapper<>,
                               Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>>,
                               Wonton::CELL>;

  using TargetGhostManager =
      Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                               Wonton::Jali_State_Wrapper,
                               Wonton::CELL>;

  using Remapper =
      Portage::CoreDriver<2, Wonton::CELL,
                          Wonton::Flat_Mesh_Wrapper<>,
                          Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>>,
                          Wonton::Jali_Mesh_Wrapper,
                          Wonton::Jali_State_Wrapper,
                          Tangram::MOF,
                          Tangram::SplitRnD<2>,
                          Tangram::ClipRnD<2>>;

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);
  Wonton::MPIExecutor_type executor(comm);

  auto source_jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 4, 4);
  auto target_jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto source_jali_state = Jali::State::create(source_jali_mesh);
  auto target_jali_state = Jali::State::create(target_jali_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh(*source_jali_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh(*target_jali_mesh);
  Wonton::Jali_State_Wrapper source_state(*source_jali_state);
  Wonton::Jali_State_Wrapper target_state(*target_jali_state);

  // - step 1: prepare source data
  int num_source_cells = source_mesh.num_entities(Wonton::CELL, Wonton::ALL);

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
  std::string const materials[num_mats] = {"mat0", "mat1", "mat2"};
  Wonton::Point<2> const lower_bound[num_mats] = {{0.0, 0.0}, {0.5, 0.0}, {0.5, 0.5}};
  Wonton::Point<2> const upper_bound[num_mats] = {{0.5, 1.0}, {1.0, 0.5}, {1.0, 1.0}};

  std::vector<int> cells[num_mats];
  std::vector<double> volume_fractions[num_mats];
  std::vector<Wonton::Point<2>> centroids[num_mats];
  std::vector<double> densities[num_mats];

  // Initialize source mesh (but don't initialize density on ghost
  // cells - we will test doing that using a ghost exchange later)
  for (int c = 0; c < num_source_cells; c++) {
    std::vector<Wonton::Point<2>> coords;
    source_mesh.cell_get_coordinates(c, &coords);

    double cellvol = source_mesh.cell_volume(c);
    double const min_vol = 1.E-6;

    Wonton::Point<2> box[2];
    BOX_INTERSECT::bounding_box<2>(coords, box, box+1);

    std::vector<double> moments;
    for (int m = 0; m < num_mats; m++) {
      bool intersected = BOX_INTERSECT::intersect_boxes<2>(lower_bound[m], upper_bound[m], box[0], box[1], &moments);
      if (intersected and moments[0] > min_vol) {  // non-trivial intersection
        cells[m].emplace_back(c);

        volume_fractions[m].emplace_back(moments[0] / cellvol);
        Wonton::Point<2> p(moments[1] / moments[0], moments[2] / moments[0]);
        if (source_mesh.cell_get_type(c) == Wonton::PARALLEL_OWNED) {
          densities[m].emplace_back((m + 1) * (m + 1) * (p[0] + p[1]));
          centroids[m].emplace_back(p);
        } else {
          densities[m].emplace_back(0.0);
          centroids[m].emplace_back(Wonton::Point<2>(0.0, 0.0));
        }
      }
    }
  }

  for (int m = 0; m < num_mats; m++) {
    source_state.add_material(materials[m], cells[m]);
    source_state.mat_add_celldata("mat_volfracs", m, volume_fractions[m].data());
    source_state.mat_add_celldata("mat_centroids", m, centroids[m].data());
    source_state.mat_add_celldata("density", m, densities[m].data());
    target_state.add_material(materials[m], {});
    target_state.mat_add_celldata<double>("density", m, 0.0);
  }


  // - step 2: Create Flat mesh/state wrapper and redistribute
  
  std::vector<std::string> source_remap_var_names = {"density"};

  Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
  Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);

  source_mesh_flat.initialize(source_mesh);
  source_state_flat.initialize(source_state, source_remap_var_names);

  Portage::MPI_Bounding_Boxes distributor(&executor);
  distributor.distribute(source_mesh_flat, source_state_flat,
                         target_mesh, target_state);


  // - step 3: Do a ghost exchange on the flat, redistributed source mesh

  SourceFlatGhostManager src_flat_ghost_manager(source_mesh_flat,
                                                source_state_flat,
                                                MPI_COMM_WORLD);

  src_flat_ghost_manager.update_ghost_values("mat_volfracs");
  src_flat_ghost_manager.update_ghost_values("mat_centroids");
  src_flat_ghost_manager.update_ghost_values("density");
  
  // - step 3: remap material fields

  Remapper remapper(source_mesh_flat, source_state_flat,
                    target_mesh, target_state, &executor);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights    = remapper.intersect_materials<Portage::IntersectRnD>(candidates);
  auto gradients  = remapper.compute_source_material_gradient("density");
  remapper.interpolate_mat_var<double, Portage::Interpolate_2ndOrder>(
    "density", "density", weights, &gradients
  );

  // fill ghost values on target mesh
  TargetGhostManager ghost_manager(target_mesh, target_state, comm);
  ghost_manager.update_ghost_values("density", true);

  for (int m = 1; m <= num_mats; ++m) {

    std::vector<std::vector<double>> send[num_ranks]; /* values of owned cells that are ghost for some rank */
    std::vector<std::vector<double>> take[num_ranks]; /* values of ghost cells that are owned by some rank */
    std::vector<int> num_sent[num_ranks];             /* number of owned cells sent by each rank */
    std::vector<int> num_take[num_ranks];             /* number of ghost cells receive from each rank */
    std::vector<MPI_Request> requests;                /* list of asynchronous MPI requests */

    send[rank] = ghost_manager.owned_values(m);
    take[rank] = ghost_manager.ghost_values(m);

    for (int i = 0; i < num_ranks; ++i) {
      num_sent[rank].emplace_back(send[rank][i].size()); // (!) may be empty
      num_take[rank].emplace_back(take[rank][i].size());
    }

    if (rank > 0) {
      MPI_Request request;
      MPI_Isend(num_sent[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
      MPI_Isend(num_take[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
    } else {
      for (int i = 1; i < num_ranks; ++i) {
        MPI_Request request;
        num_sent[i].resize(num_ranks, 0);
        MPI_Irecv(num_sent[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
        num_take[i].resize(num_ranks, 0);
        MPI_Irecv(num_take[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    requests.clear();

    if (rank > 0) {
      int const offset = num_ranks * num_ranks;
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          MPI_Request request;
          int tag = num_ranks * (rank - 1) + i;
          // nb: can be empty for this material
          if (not send[rank][i].empty()) {
            MPI_Isend(send[rank][i].data(), send[rank][i].size(), MPI_DOUBLE, 0, tag, comm, &request);
            requests.emplace_back(request);
          }
          if (not take[rank][i].empty()) {
            MPI_Isend(take[rank][i].data(), take[rank][i].size(), MPI_DOUBLE, 0, tag + offset, comm, &request);
            requests.emplace_back(request);
          }
        }
      }
    } else {
      int const offset = num_ranks * num_ranks;
      for (int i = 1; i < num_ranks; ++i) {
        send[i].resize(num_ranks);
        take[i].resize(num_ranks);
        for (int j = 0; j < num_ranks; ++j) {
          if (i != j) {
            MPI_Request request;
            int tag = num_ranks * (i - 1) + j;
            // nb: can be empty for this material
            if (num_sent[i][j]) {
              send[i][j].resize(num_sent[i][j]);
              MPI_Irecv(send[i][j].data(), num_sent[i][j], MPI_DOUBLE, i, tag, comm, &request);
              requests.emplace_back(request);
            }
            if (num_take[i][j]) {
              take[i][j].resize(num_take[i][j]);
              MPI_Irecv(take[i][j].data(), num_take[i][j], MPI_DOUBLE, i, tag + offset, comm, &request);
              requests.emplace_back(request);
            }
          }
        }
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

    if (rank == 0) {
      for (int i = 0; i < num_ranks; ++i) {
        for (int j = 0; j < num_ranks; ++j) {
          if (i == j) { continue; }
          int const num_owned = send[i][j].size();
          int const num_ghost = take[j][i].size();
          ASSERT_EQ(num_owned, num_ghost);
          for (int k = 0; k < num_owned; ++k) {
            ASSERT_DOUBLE_EQ(send[i][j][k], take[j][i][k]);
          }
        }
      }
    }

    MPI_Barrier(comm);
  } // end for each material
}
#endif
