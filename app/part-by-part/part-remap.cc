/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "remap.h"

/**
 * @brief Run the application.
 *
 * @param argc: arguments count
 * @param argv: arguments values
 * @return status code
 */
int main(int argc, char* argv[]) {

  auto tic = timer::now();

  Params params;
  auto& my_rank  = params.rank;
  auto& nb_ranks = params.nb_ranks;
  auto& comm     = params.comm;

  // init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nb_ranks);
  MPI_Comm_rank(comm, &my_rank);

  if (my_rank == 0) {
    std::printf(" ---------------------------------------------------------- \n");
    std::printf("  Demomstration application for part-by-part remap.         \n");
    std::printf("  It is intended to be used to:                             \n");
    std::printf("  - preserve pure cells in multi-material context           \n");
    std::printf("  - avoid diffusion effects on discontinuous fields         \n");
    std::printf("  - speedup remap when only a small delimited area changed  \n");
    std::printf(" ---------------------------------------------------------- \n");
  }

  // check and parse parameters
  if (params.parse(argc, argv) == EXIT_FAILURE)
    return EXIT_FAILURE;

  /* ------------------------------------------------------------------------ */
  if (my_rank == 0)
    std::printf("Initializing distributed meshes ... ");

  // source and target meshes and states
  std::shared_ptr<Jali::Mesh>  source_mesh;
  std::shared_ptr<Jali::Mesh>  target_mesh;
  std::shared_ptr<Jali::State> source_state;
  std::shared_ptr<Jali::State> target_state;

  // mesh regions, geometric models and factories
  std::set<int> source_regions_indices;
  std::set<int> target_regions_indices;
  JaliGeometry::GeometricModelPtr source_model = nullptr;
  JaliGeometry::GeometricModelPtr target_model = nullptr;
  Jali::MeshFactory source_factory(comm);
  Jali::MeshFactory target_factory(comm);

  // concatenate all source regions
  for (auto&& field : params.parts) {
    for (auto&& part : field.second) {
      for (auto&& i : part.source.blocks) { source_regions_indices.insert(i); }
      for (auto&& i : part.target.blocks) { target_regions_indices.insert(i); }
    }
  }

  // check if we need to create a geometric model so that we could load
  // part definitions (cell sets) from the file.
  if (not source_regions_indices.empty()) {
    JaliGeometry::RegionVector source_regions;
    for (auto& id : source_regions_indices) {
      auto region = new JaliGeometry::LabeledSetRegion("part_" + std::to_string(id),
                                                       id, "Entity_kind::CELL",
                                                       params.source, "Exodus II",
                                                       std::to_string(id));
      source_regions.emplace_back(region);
    }
    source_regions_indices.clear();
    source_model = new JaliGeometry::GeometricModel(params.dimension, source_regions);
  }

  if (not target_regions_indices.empty()) {
    JaliGeometry::RegionVector target_regions;
    for (auto& id : target_regions_indices) {
      auto region = new JaliGeometry::LabeledSetRegion("part_" + std::to_string(id),
                                                       id, "Entity_kind::CELL",
                                                       params.target, "Exodus II",
                                                       std::to_string(id));
      target_regions.emplace_back(region);
    }
    target_regions_indices.clear();
    target_model = new JaliGeometry::GeometricModel(params.dimension, target_regions);
  }

  // load both distributed meshes
  source_factory.included_entities(Jali::Entity_kind::ALL_KIND);
  source_factory.partitioner(Jali::Partitioner_type::METIS);
  if (source_model)
    source_factory.geometric_model(source_model);

  target_factory.included_entities(Jali::Entity_kind::ALL_KIND);
  target_factory.partitioner(Jali::Partitioner_type::METIS);
  if (target_model)
    target_factory.geometric_model(target_model);

  source_mesh  = source_factory(params.source);
  target_mesh  = target_factory(params.target);
  source_mesh->init_sets_from_geometric_model(); // needed for querying sets
  target_mesh->init_sets_from_geometric_model();
  source_state = Jali::State::create(source_mesh);
  target_state = Jali::State::create(target_mesh);
  source_state->init_from_mesh(); // import any state data from mesh


  // interfaces with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  Wonton::MPIExecutor_type mpi_executor(comm);
  Wonton::Executor_type *executor = (nb_ranks > 1 ? &mpi_executor : nullptr);

  // ensure that source and target mesh have the same dimension,
  // and that it corresponds to the one specified by the user.
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  assert(unsigned(params.dimension) == source_mesh->space_dimension());

  // retrieve mesh resolutions
  // nb: only consider owned cells for source mesh to avoid errors.
  long const nb_source_cells = source_mesh_wrapper.num_owned_cells();
  long const nb_target_cells = target_mesh_wrapper.num_owned_cells();

  MPI_Barrier(comm);

  if (my_rank == 0)
    std::printf(" done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic, true));

  /* ------------------------------------------------------------------------ */
  if (my_rank == 0)
    std::printf("\nRunning part-by-part remap ... \n");

  long total_source_cells = 0;
  long total_target_cells = 0;
  MPI_Reduce(&nb_source_cells, &total_source_cells, 1, MPI_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(&nb_target_cells, &total_target_cells, 1, MPI_LONG, MPI_SUM, 0, comm);

  // for formatting
  int format = params.get_number_digit(std::max(total_source_cells, total_target_cells));

  // print some infos for the user
  if (my_rank == 0) {
    std::printf(" - source mesh has %*ld cells.\n", format, total_source_cells);
    std::printf(" - target mesh has %*ld cells.\n", format, total_target_cells);
    std::printf(" - specified numerical fields: \n");
    for (auto&& field : params.fields)
      std::printf("   \u2022 %s: \e[32m%s\e[0m\n", field.first.data(), field.second.data());

    std::printf("\n");
  }

  // assign scalar fields - or use existing ones
  for (auto&& field : params.fields) {
    double* field_data = nullptr;
    // add scalar field to both meshes
    target_state_wrapper.mesh_add_data<double>(Wonton::CELL, field.first, 0.);

    if (field.second != "internal") {
      source_state_wrapper.mesh_add_data<double>(Wonton::CELL, field.first, 0.);
      source_state_wrapper.mesh_get_data(Wonton::CELL, field.first, &field_data);

      // evaluate field expression and assign it
      user_field_t source_field;
      if(source_field.initialize(params.dimension, field.second)) {
        for (long c = 0; c < nb_source_cells; c++) {
          field_data[c] = source_field(source_mesh->cell_centroid(c));
          #if DEBUG_PART_BY_PART
            std::printf("field_data[%d]: %.3f\n", c, field_data[c]);
          #endif
        }
      } else
        return params.abort("cannot parse numerical field "+ field.second, false);
    } // no need to evaluate otherwise
  }

  MPI_Barrier(comm);

  std::vector<std::vector<int>> source_cells;
  std::vector<std::vector<int>> target_cells;

  // remap each field
  for (auto&& entry : params.fields) {
    std::string field = entry.first;
    int const nb_parts = params.parts[field].size();
    assert(nb_parts > 1);

    // first determine source and target parts for this field.
    filter_t filter;
    source_cells.resize(nb_parts);
    target_cells.resize(nb_parts);

    if (my_rank == 0)
      std::printf(" - Remapping \e[32m%s\e[0m field on %d parts:\n", field.data(), nb_parts);

    long max_source_parts = 0;
    long max_target_parts = 0;
    long total_source_part[nb_parts];
    long total_target_part[nb_parts];
    std::fill(total_source_part, total_source_part + nb_parts, 0);
    std::fill(total_target_part, total_target_part + nb_parts, 0);

    // Filter source and target cells for each part
    // with respect to user-defined predicates.
    for (int i = 0; i < nb_parts; ++i) {
      auto const& part = params.parts[field][i];
      source_cells[i].reserve(nb_source_cells);
      target_cells[i].reserve(nb_target_cells);

      if (not part.source.blocks.empty()) {
        for (auto&& region : part.source.blocks) {
          std::vector<int> cells;
          source_mesh->get_set_entities("part_"+std::to_string(region),
                                        Jali::Entity_kind::CELL,
                                        Jali::Entity_type::PARALLEL_OWNED, &cells);
          for (auto&& c : cells) { source_cells[i].emplace_back(c); }
        }
      } else if(not part.source.expr.empty()) {
        // populate source part entities
        if (filter.initialize(params.dimension, part.source.expr)) {
          for (long s = 0; s < nb_source_cells; ++s) {
            if (filter(source_mesh->cell_centroid(s)))
              source_cells[i].push_back(s);
          }
        } else
          return params.abort("cannot filter source part cells for field " + field, false);
      } else
        return params.abort("undefined source part");

      if (not part.target.blocks.empty()) {
        for (auto&& region : part.target.blocks) {
          std::vector<int> cells;
          target_mesh->get_set_entities("part_"+std::to_string(region),
                                        Jali::Entity_kind::CELL,
                                        Jali::Entity_type::PARALLEL_OWNED, &cells);
          for (auto&& c : cells) { target_cells[i].emplace_back(c); }
        }
      } else if(not part.target.expr.empty()) {
        // populate target part entities
        if (filter.initialize(params.dimension, part.target.expr)) {
          for (long t = 0; t < nb_target_cells; ++t) {
            if (filter(target_mesh->cell_centroid(t)))
              target_cells[i].push_back(t);
          }
        } else
          return params.abort("cannot filter target part cells for field "+field, false);
      } else
        return params.abort("undefined target part");

      source_cells[i].shrink_to_fit();
      target_cells[i].shrink_to_fit();

      long local_source_part = source_cells[i].size();
      long local_target_part = target_cells[i].size();
      MPI_Reduce(&local_source_part, total_source_part + i, 1, MPI_LONG, MPI_SUM, 0, comm);
      MPI_Reduce(&local_target_part, total_target_part + i, 1, MPI_LONG, MPI_SUM, 0, comm);

      max_source_parts = std::max(total_source_part[i], max_source_parts);
      max_target_parts = std::max(total_target_part[i], max_target_parts);
    }

    MPI_Barrier(comm);

    if (my_rank == 0) {
      int const source_digits = params.get_number_digit(max_source_parts);
      int const target_digits = params.get_number_digit(max_target_parts);

      for (int i = 0; i < nb_parts; ++i) {
        std::printf(
          "   \u2022 ["
          " source: %*ld cells \e[32m(%5.2f %%)\e[0m,"
          " target: %*ld cells \e[32m(%5.2f %%)\e[0m ]\n",
          source_digits, total_source_part[i],
          100. * double(total_source_part[i]) / double(total_source_cells),
          target_digits, total_target_part[i],
          100. * double(total_target_part[i]) / double(total_target_cells)
        );
      }
      std::printf("\n");
    }

    // then process part-by-part remapping.
    // need to explicitly instantiate the driver due to template arguments
    // forwarding failure when instantiating search and intersection kernels.
    // part-by-part remap each field
    switch (params.dimension) {
      case 2: remap<2>(field, nb_parts, source_mesh, target_mesh,
                       source_mesh_wrapper, target_mesh_wrapper,
                       source_state_wrapper, target_state_wrapper,
                       executor, source_cells, target_cells, params); break;

      case 3: remap<3>(field, nb_parts, source_mesh, target_mesh,
                       source_mesh_wrapper, target_mesh_wrapper,
                       source_state_wrapper, target_state_wrapper,
                       executor, source_cells, target_cells, params); break;

      default: return params.abort("invalid dimension", false);
    }

    source_cells.clear();
    target_cells.clear();
    MPI_Barrier(comm);

    if (my_rank == 0) {
      std::fflush(stdout);
      std::printf(" %s \n", std::string(58,'-').data());
    }
  }

  if (my_rank == 0)
    std::printf("Remap done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic, true));

  /* ------------------------------------------------------------------------ */
  if (params.kind == Wonton::Entity_kind::CELL) {
    // compute error for each field
    for (auto&& field : params.fields) {

      /*
       * skip error computation for internal fields,
       * there is no way to assess the exact value
       * at a given target cell centroid.
       */
      if (field.second == "internal") { continue; }

      double error_l1 = 0;
      double error_l2 = 0;
      double relative_error = 0;
      user_field_t exact_value;

      if (my_rank == 0)
        std::printf("\nComputing error for \e[32m%s\e[0m field ... ", field.first.data());

      if (exact_value.initialize(params.dimension, field.second)) {
        // keep track of min and max field values on both meshes
        double min_source_val = std::numeric_limits<double>::max();
        double max_source_val = std::numeric_limits<double>::min();
        double min_target_val = min_source_val;
        double max_target_val = max_source_val;
        double source_extents[] = { min_source_val, max_source_val };
        double target_extents[] = { min_target_val, max_target_val };
        double source_mass = 0;
        double target_mass = 0;
        double total_mass[] = { 0, 0 };
        double global_error[] = { 0, 0, 0 };

        // retrieve field data on both meshes
        double* source_field_data = nullptr;
        double* target_field_data = nullptr;
        source_state_wrapper.mesh_get_data(Wonton::CELL, field.first, &source_field_data);
        target_state_wrapper.mesh_get_data(Wonton::CELL, field.first, &target_field_data);

        // compute total mass on the source mesh to check conservation
        for (int s = 0; s < nb_source_cells; ++s) {
          min_source_val = std::min(source_field_data[s], min_source_val);
          max_source_val = std::max(source_field_data[s], max_source_val);
          source_mass += source_field_data[s] * source_mesh_wrapper.cell_volume(s);
        }

        // compute cell error
        for (int t = 0; t < nb_target_cells; ++t) {
          min_target_val = std::min(target_field_data[t], min_target_val);
          max_target_val = std::max(target_field_data[t], max_target_val);
          // compute difference between exact and remapped value
          auto const& centroid = target_mesh->cell_centroid(t);
          auto const error = exact_value(centroid) - target_field_data[t];
          auto const cell_volume = target_mesh_wrapper.cell_volume(t);
          // update L^p norm error and target mass
          error_l1 += std::abs(error) * cell_volume;
          error_l2 += error * error * cell_volume;
          relative_error += std::abs(target_field_data[t]) * cell_volume;
          target_mass += target_field_data[t] * cell_volume;
        }

        // accumulate all local values on rank 0
        MPI_Reduce(&error_l1, global_error, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&error_l2, global_error+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&relative_error, global_error+2, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&min_source_val, source_extents, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&max_source_val, source_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&min_target_val, target_extents, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&max_target_val, target_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&source_mass, total_mass, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&target_mass, total_mass+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

        // update local values then
        error_l1 = global_error[0];
        error_l2 = sqrt(global_error[1]);
        relative_error = error_l1 / global_error[2];
        source_mass = total_mass[0];
        target_mass = total_mass[1];

        MPI_Barrier(comm);

        if (my_rank == 0) {
          std::printf( " done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic, true));
          std::printf(" \u2022 L1-norm error     = %lf\n", error_l1);
          std::printf(" \u2022 L2-norm error     = %lf\n", error_l2);
          std::printf(" \u2022 relative L1 error = %.15f\n", relative_error);
          std::printf(" \u2022 source values     = [%.15f, %.15f]\n",
            source_extents[0], source_extents[1]);
          std::printf(" \u2022 target values     = [%.15f, %.15f]\n",
            target_extents[0], target_extents[1]);
          std::printf(" \u2022 source total mass = %.15f\n", source_mass);
          std::printf(" \u2022 target total mass = %.15f\n", target_mass);
          std::printf(" \u2022 mass discrepancy  = %.15f\n",
            std::abs(source_mass - target_mass));
        }
      } else
        return params.abort("cannot parse numerical field "+ field.second, false);
    }
  } else
    return params.abort("part-by-part node remap is not supported", false);

  MPI_Barrier(comm);

  /* ------------------------------------------------------------------------ */
  if (params.dump) {
    if (my_rank == 0)
      std::printf("\nDump data to exodus files ... ");

    // dump meshes with attached data
    source_state->export_to_mesh();
    target_state->export_to_mesh();
    source_mesh->write_to_exodus_file("source.exo");
    target_mesh->write_to_exodus_file("target.exo");

    MPI_Barrier(comm);

    if (my_rank == 0)
      std::printf(" done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));
  }

  if (my_rank == 0)
    std::printf("\nDumping remapped fields ... ");

  std::vector<int> index_helper;
  std::vector<int> local_index;
  std::vector<int> global_index;
  std::vector<double> local_field;
  std::vector<double> global_field;
  double* field_data = nullptr;

  local_index.resize(nb_target_cells);
  local_field.resize(nb_target_cells);

  // dump each field in a separate file
  for (auto&& field : params.fields) {
    // retrieve cell global indices and field values for current rank
    for (auto c = 0; c < nb_target_cells; ++c)
      local_index[c] = target_mesh->GID(c, Jali::Entity_kind::CELL);

    target_state_wrapper.mesh_get_data(Wonton::CELL, field.first, &field_data);
    std::copy(field_data, field_data + nb_target_cells, local_field.begin());

    // append local index and values lists to master rank global lists
    Portage::collate(comm, my_rank, nb_ranks, local_index, global_index);
    Portage::collate(comm, my_rank, nb_ranks, local_field, global_field);

    if (my_rank == 0) {
      // sort field values by global ID
      Portage::argsort(global_index, index_helper);
      Portage::reorder(global_index, index_helper);
      Portage::reorder(global_field, index_helper);

      // dump sorted data eventually
      std::ofstream file(params.results+"_"+ field.first +".dat");
      file << std::scientific;
      file.precision(17);

      for (long c = 0; c < nb_target_cells; ++c)
        file << global_index[c] << "\t" << global_field[c] << std::endl;

      index_helper.clear();
      global_index.clear();
      global_field.clear();
    }

    MPI_Barrier(comm);
  }

  if (my_rank == 0)
    std::printf(" done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));

  delete source_model;  // no effect if nullptr
  delete target_model;

  MPI_Finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
