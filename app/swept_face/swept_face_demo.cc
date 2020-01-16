/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <limits>
#include <mpi.h>

#include "Mesh.hh"         // see https://github.com/lanl/jali
#include "MeshFactory.hh"
#include "JaliState.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/driver/coredriver.h"
#include "portage/search/search_kdtree.h"
#include "portage/search/search_swept_face.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/support/portage.h"
#include "portage/support/timer.h"
#include "user_field.h" // parsing and evaluating user defined expressions

#define ENABLE_TIMINGS 1

#ifdef ENABLE_PROFILE
  #include "ittnotify.h"
#endif
#if ENABLE_TIMINGS
  #include "portage/support/timer.h"
#endif

/**
 * @file swept_face_demo.cc
 *
 * @brief A simple application that remaps using a swept face algorithm
 *        for single material field between two meshes which can be internally
 *        generated cartesian grids or externally read unstructured meshes.
 */

#define DEBUG_PRINT 0

static int rank = 0;
static int numpe = 1;
static MPI_Comm comm = MPI_COMM_WORLD;

/**
 * @brief Identify boundary nodes of the given mesh.
 *
 * @param mesh the current mesh pointer.
 * @return a flag array to check if a given mesh node is a boundary one.
 */
std::vector<bool> identify_exterior_boundary_nodes(std::shared_ptr<Jali::Mesh> mesh);

/**
 * @brief Move the given target mesh points.
 *
 * It moves the target points to obtain a target mesh with same
 * connectivity but different point positions: loop over all the
 * boundary vertices assuming that we are only dealing with internally
 * generated meshes.
 * x_new = x_old + x_velocity(tcur) * deltaT
 * y_new = y_old + y_velocity(tcur) * deltaT
 *
 * @param mesh    current mesh pointer
 * @param iter    current timestep
 * @param deltaT  displacement step
 * @param periodT displacement period
 * @param scale   scaling factor
 */
template<int dim>
void move_target_mesh_nodes(std::shared_ptr<Jali::Mesh> mesh,
                            int iter, double& deltaT, double& periodT, int& scale);

/**
 * @brief Compute a single vortex velocity value at the given point.
 *
 * @param[in]  coords  coordinates of the current point.
 * @param[in]  tcur    scaled displacement step
 * @param[in]  periodT displacement period.
 * @param[out] veloc   computed velocity for current point.
 */
void compute_single_vortex_velocity(double* coords, double& tcur,
                                    double& periodT, double* veloc);

/**
 * @brief Remap the analytically imposed field and output related errors.
 *
 * @tparam dim             dimension of the problem.
 * @param source_mesh      source mesh pointer.
 * @param target_mesh      target mesh pointer.
 * @param field_expression analytical expression of the field to remap.
 * @param interp_order     order of accuracy for interpolation.
 * @param limiter          gradient limiter to use for internal cells.
 * @param bnd_limiter      gradient limiter to use for boundary cells.
 * @param mesh_output      dump meshes or not?
 * @param field_filename   field file name for imported meshes.
 * @param iteration        current iteration.
 * @param L1_error         L1-norm error.
 * @param L2_error         L2-norm error.
 * @param comm             the MPI communicator to use.
 * @param profiler         profiler object pointer.
 */
template<int dim>
void remap(std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh,
           std::string field_expression,
           int interp_order,
           Portage::Limiter_type limiter,
           Portage::Boundary_Limiter_type bnd_limiter,
           bool intersect_based, bool mesh_output,
           std::string field_filename, int iteration,
           double& L1_error, double& L2_error,
           std::shared_ptr<Profiler> profiler);

/**
 * @brief Print an error message followed by command-line usage, and exits.
 *
 * @param message the error message to be displayed
 * @return status code
 */
int abort(std::string message);

/**
 * @brief Print some infos to the user.
 *
 * @param source_mesh      source mesh pointer
 * @param target_mesh      target mesh pointer
 * @param field_expression expression of the field to remap
 * @param interp_order     order of interpolation
 * @param intersect_based  use intersection-based method
 * @param limiter          gradient limiter for internal cells
 * @param bnd_limiter      gradient limiter for boundary cells
 */
void print_infos(std::shared_ptr<Jali::Mesh> source_mesh,
                 std::shared_ptr<Jali::Mesh> target_mesh,
                 std::string field_expression,
                 int interp_order, bool intersect_based,
                 Portage::Limiter_type limiter,
                 Portage::Boundary_Limiter_type bnd_limiter);

/**
 * @brief Print command-line usage.
 *
 */
void print_usage() {

  std::cerr << "Usage: swept_face_demo [options]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "\t--help                  show this help message and exit"    << std::endl;
  std::cerr << "\t--dim           INT     dimension of the problem"           << std::endl;
  std::cerr << "\t--ncells        INT     cells per axis for generated grids" << std::endl;
  std::cerr << "\t--source_file   STRING  source mesh file to import"         << std::endl;
  std::cerr << "\t--target_file   STRING  target mesh file to import"         << std::endl;
  std::cerr << "\t--remap_order   INT     order of interpolation"             << std::endl;
  std::cerr << "\t--intersect     CHAR    use intersection-based remap [y|n]" << std::endl;
  std::cerr << "\t--field         STRING  numerical field to remap"           << std::endl;
  std::cerr << "\t--field_file    STRING  numerical field file to import"     << std::endl;
  std::cerr << "\t--ntimesteps    INT     number of timesteps"                << std::endl;
  std::cerr << "\t--scale_by      FLOAT   displacement scaling factor"        << std::endl;
  std::cerr << "\t--limiter       STRING  gradient limiter to use"            << std::endl;
  std::cerr << "\t--bnd_limiter   STRING  gradient limiter for boundary"      << std::endl;
  std::cerr << "\t--output_meshes CHAR    dump meshes [y|n]"                  << std::endl;
#if ENABLE_TIMINGS
  std::cerr << "\t--scaling       STRING  scaling study [strong|weak]"        << std::endl;
  std::cerr << "\t--only_threads  CHAR    thread scaling profiling [y|n]"     << std::endl;
#endif
}


/**
 * @brief Main method.
 *
 * @param argc number of arguments
 * @param argv arguments values
 * @return status code
 */
int main(int argc, char** argv) {
  // Pause profiling until main loop
#if ENABLE_PROFILE
  __itt_pause();
#endif

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &numpe);
  MPI_Comm_rank(comm, &rank);

  if (argc == 1) {
    print_usage();
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  int dim = 2;
  int ncells = 0;
  int interp_order = 1;
  int ntimesteps = 4;
  int scale = 20;
  bool mesh_output = false;
  bool intersect_based = false;

  std::string source_file;
  std::string target_file;
  std::string field_expression;
  std::string field_path;  // No default;

  auto limiter     = Portage::Limiter_type::NOLIMITER;
  auto bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;

#if ENABLE_TIMINGS
  bool only_threads = false;
  std::string scaling_type = "strong";
#endif

  // Parse the input
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of("=");
    std::string keyword = arg.substr(keyword_beg, keyword_end-keyword_beg);
    std::string valueword = arg.substr(keyword_end+1, len-(keyword_end+1));

    if (keyword == "dim") {
      dim = stoi(valueword);
      assert(dim == 2 || dim == 3);
    } else if (keyword == "ncells") {
      ncells = stoi(valueword);
    } else if (keyword == "source_file") {
      source_file = valueword;
    } else if (keyword == "target_file") {
      target_file = valueword;
    } else if (keyword == "scale_by") {
      scale = stoi(valueword);
    } else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
    } else if (keyword == "intersect") {
      intersect_based = (valueword == "y");
    } else if (keyword == "field") {
      field_expression = valueword;
    } else if (keyword == "field_file") {
      field_path = valueword;
    } else if (keyword == "limiter") {
      if (valueword == "barth_jespersen" or valueword == "BARTH_JESPERSEN")
        limiter = Portage::Limiter_type::BARTH_JESPERSEN;
    } else if (keyword == "bnd_limiter") {
      if (valueword == "zero_gradient" or valueword == "ZERO_GRADIENT")
        bnd_limiter = Portage::Boundary_Limiter_type::BND_ZERO_GRADIENT;
      else if (valueword == "barth_jespersen" or valueword == "BARTH_JESPERSEN")
        bnd_limiter = Portage::Boundary_Limiter_type::BND_BARTH_JESPERSEN;
    } else if (keyword == "ntimesteps") {
      ntimesteps = stoi(valueword);
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
#if ENABLE_TIMINGS
    } else if (keyword == "only_threads"){
      only_threads = (numpe == 1 and valueword == "y");
    } else if (keyword == "scaling") {
      assert(valueword == "strong" or valueword == "weak");
      scaling_type = valueword;
#endif
    } else if (keyword == "help") {
      return abort("");
    } else
      return abort("unrecognized option "+ keyword);
  }

  // Some input error checking
  bool invalid_source = (ncells > 0 and not source_file.empty()) or (not ncells and source_file.empty());
  bool invalid_target = (ncells > 0 and not target_file.empty()) or (not ncells and target_file.empty());
  bool no_field_expression = ncells > 0 and field_expression.empty();

  if (invalid_source)
    return abort("please choose either generated or imported source mesh");

  if (invalid_target)
    return abort("please choose either generated or imported target mesh");

  if (no_field_expression)
    return abort("no field imposed on internally generated source mesh");

#if ENABLE_TIMINGS
  auto profiler = std::make_shared<Profiler>();
  // save params for after
  profiler->params.ranks   = numpe;
  profiler->params.nsource = std::pow(ncells, dim);
  profiler->params.ntarget = std::pow(ncells, dim);
  profiler->params.order   = interp_order;
  profiler->params.nmats   = 1;
  profiler->params.output += scaling_type + "_scaling_" +
                             std::string(only_threads ? "omp.dat": "mpi.dat");
  #if defined(_OPENMP)
    profiler->params.threads = omp_get_max_threads();
  #endif
  // start timers here
  auto start = timer::now();
  auto tic = start;
#else
  std::shared_ptr<Profiler> profiler = nullptr;
#endif

  if (rank == 0) {
    std::cout << " ------------------------------------------------------------ " << std::endl;
    std::cout << "  A simple application for single material field remap using  " << std::endl;
    std::cout << "  a swept-face algorithm. Meshes can be internally generated  " << std::endl;
    std::cout << "  cartesian grids or externally imported unstructured meshes. " << std::endl;
    std::cout << " ------------------------------------------------------------ " << std::endl;
    std::cout << (ncells > 0 ? "Generate" : "Import") << " meshes ... " << std::flush;
  }

  std::shared_ptr<Jali::Mesh> source_mesh;
  std::shared_ptr<Jali::Mesh> target_mesh;

  Jali::MeshFactory factory(comm);
  factory.included_entities(Jali::Entity_kind::ALL_KIND);

  if (ncells == 0) {
    // import meshes while using a graph partitioner to partition them.
    factory.partitioner(Jali::Partitioner_type::METIS);
    source_mesh = factory(source_file);
    target_mesh = factory(target_file);
  } else {
    // generate distributed source and target meshes
    // using a block partitioner to partition them.
    factory.partitioner(Jali::Partitioner_type::BLOCK);
    // no boundary mismatch
    double x_min = 0.0, x_max = 1.0;
    double y_min = 0.0, y_max = 1.0;
    // Create the source and target meshes
    if (dim == 2) {
      source_mesh = factory(x_min, y_min, x_max, y_max, ncells, ncells);
      target_mesh = factory(x_min, y_min, x_max, y_max, ncells, ncells);
    }
    else if (dim == 3) {
      double z_min = 0.0, z_max = 1.0;
      source_mesh = factory(x_min, y_min, z_min, x_max, y_max, z_max, ncells, ncells, ncells);
      target_mesh = factory(x_min, y_min, z_min, x_max, y_max, z_max, ncells, ncells, ncells);
    }
  }

  // make sure that we have the right dimension and that source and
  // target mesh dimensions match (important when one or both of the
  // meshes are read in)
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  assert(dim == source_mesh->space_dimension());

#if ENABLE_TIMINGS
  profiler->time.mesh_init = timer::elapsed(tic);
  if (rank == 0)
    std::printf("done. \e[32m(%.3f s)\e[0m\n", profiler->time.mesh_init);
#else
  if (rank == 0)
    std::cout << "done" << std::endl;
#endif

  // Output some information for the user
  print_infos(source_mesh, target_mesh,
              field_expression, interp_order,
              intersect_based, limiter, bnd_limiter);

  double periodT = 2.0;
  double deltaT = periodT/ntimesteps;

  std::vector<double> l1_err(ntimesteps, 0.0);
  std::vector<double> l2_err(ntimesteps, 0.0);

  for (int i = 1; i < ntimesteps; i++) {
    if (rank == 0)
      std::cout << "------------- timestep "<< i << " -------------" << std::endl;

    // move nodes of the target mesh then run the remap and output related errors
    try {
      switch (dim) {
        case 2:
          move_target_mesh_nodes<2>(target_mesh, i, deltaT, periodT, scale);
          remap<2>(source_mesh, target_mesh, field_expression, interp_order,
                   limiter, bnd_limiter, intersect_based, mesh_output, field_path,
                   i, l1_err[i], l2_err[i], profiler); break;
        case 3:
          move_target_mesh_nodes<3>(target_mesh, i, deltaT, periodT, scale);
          remap<3>(source_mesh, target_mesh, field_expression, interp_order,
                   limiter, bnd_limiter, intersect_based, mesh_output, field_path,
                   i, l1_err[i], l2_err[i], profiler); break;
        default:
          return abort("invalid dimension");
      }
    } catch (std::exception const& exception) {
      return abort(exception.what());
    }

    if (rank == 0)
      std::cout << std::endl;
  }

#if ENABLE_TIMINGS
  profiler->time.total = timer::elapsed(start);

  // dump timing data
  if (rank == 0) {
    profiler->dump();
  }
#endif

  MPI_Finalize();
  return EXIT_SUCCESS;
}

/**
 * @brief Remap the analytically imposed field and output related errors.
 *
 * @tparam dim             dimension of the problem.
 * @param source_mesh      source mesh pointer.
 * @param target_mesh      target mesh pointer.
 * @param field_expression analytical expression of the field to remap.
 * @param interp_order     order of accuracy for interpolation.
 * @param limiter          gradient limiter to use for internal cells.
 * @param bnd_limiter      gradient limiter to use for boundary cells.
 * @param mesh_output      dump meshes or not?
 * @param field_filename   field file name for imported meshes.
 * @param iteration        current iteration.
 * @param L1_error         L1-norm error.
 * @param L2_error         L2-norm error.
 * @param comm             the MPI communicator to use.
 * @param profiler         profiler object pointer.
 */
template<int dim>
void remap(std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh,
           std::string field_expression,
           int interp_order,
           Portage::Limiter_type limiter,
           Portage::Boundary_Limiter_type bnd_limiter,
           bool intersect_based, bool mesh_output,
           std::string field_filename, int iteration,
           double& L1_error, double& L2_error,
           std::shared_ptr<Profiler> profiler) {

  if (rank == 0)
    std::cout << "Remap field ... " << std::flush;

  // the remapper to use
  using Remapper = Portage::CoreDriver<dim,
                                       Wonton::Entity_kind::CELL,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Wonton::Jali_State_Wrapper>;
  // Executor
  Wonton::MPIExecutor_type mpi_executor(comm);
  Wonton::Executor_type* executor = (numpe > 1 ? &mpi_executor : nullptr);

  // Native jali state managers for source and target
  std::shared_ptr<Jali::State> source_state(Jali::State::create(source_mesh));
  std::shared_ptr<Jali::State> target_state(Jali::State::create(target_mesh));

  // wrappers for interfacing with underlying mesh data structures,
  // and for source and target fields.
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  const int nsrccells = source_mesh_wrapper.num_owned_cells() +
                        source_mesh_wrapper.num_ghost_cells();
  const int ntarcells = target_mesh_wrapper.num_owned_cells() +
                        target_mesh_wrapper.num_ghost_cells();

#if ENABLE_TIMINGS
  auto tic = timer::now();
#endif

  // Compute field data on source, and add to state manager.
  user_field_t exact_value;
  if (not exact_value.initialize(dim, field_expression))
    MPI_Abort(comm, EXIT_FAILURE);

  std::vector<double> source_data(nsrccells);
  for (int c = 0; c < nsrccells; ++c) {
    source_data[c] = exact_value(source_mesh->cell_centroid(c));
  }

  source_state->add("density", source_mesh,
                    Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL,
                    source_data.data());

  target_state->add<double, Jali::Mesh, Jali::UniStateVector>("density", target_mesh,
                                                              Jali::Entity_kind::CELL,
                                                              Jali::Entity_type::ALL,
                                                              0.0);
  MPI_Barrier(comm);

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper, executor);

  auto candidates = (intersect_based ? remapper.template search<Portage::SearchKDTree>()
                                     : remapper.template search<Portage::SearchSweptFace>());
  auto weights = (intersect_based
    ? (dim == 2 ? remapper.template intersect_meshes<Portage::IntersectR2D>(candidates)
                : remapper.template intersect_meshes<Portage::IntersectR3D>(candidates))
    : (dim == 2 ? remapper.template intersect_meshes<Portage::IntersectSweptFace2D>(candidates)
                : remapper.template intersect_meshes<Portage::IntersectSweptFace3D>(candidates)));

  switch (interp_order) {
    case 1:
      remapper.template interpolate_mesh_var<double, Portage::Interpolate_1stOrder>("density",
                                                                                    "density",
                                                                                    weights);
      break;
    case 2: {
      auto gradients = remapper.compute_source_gradient("density", limiter, bnd_limiter);
      remapper.template interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>("density",
                                                                                    "density",
                                                                                    weights, &gradients);
      break;
    }
    default: break;
  }

#if ENABLE_TIMINGS
  profiler->time.remap += timer::elapsed(tic);

  if (rank == 0)
    std::printf("done. \e[32m(%.3f s)\e[0m\n", profiler->time.remap);
#else
  if (rank == 0)
    std::cout << "done" << std::endl;
#endif

  // dump meshes if requested
  if (mesh_output) {
    if (rank == 0)
      std::cout << "Dump data ... " << std::flush;

    std::string suffix = std::to_string(rank) + std::to_string(iteration)+ ".exo";
    source_state->export_to_mesh();
    target_state->export_to_mesh();
    source_mesh->write_to_exodus_file("source_" + suffix);
    target_mesh->write_to_exodus_file("target_" + suffix);

    if (rank == 0)
      std::cout << "done." << std::endl;
  }

  if (rank == 0)
    std::cout << "Extract stats ... " << std::flush;

  // Compute error
  L1_error = 0.0;
  L2_error = 0.0;

  //double error;
  double const* target_data_field;
  double const* source_data_field;
  double min_source_val = std::numeric_limits<double>::max();
  double max_source_val = std::numeric_limits<double>::min();
  double min_target_val = min_source_val;
  double max_target_val = max_source_val;
  double source_extents[] = { min_source_val, max_source_val };
  double target_extents[] = { min_target_val, max_target_val };
  double target_mass = 0.;
  double source_mass = 0.;
  double total_mass[] = { 0., 0. };
  double global_error[] = { 0., 0. };
  double total_volume = 0. ;

  source_state_wrapper.mesh_get_data<double>(Portage::CELL, "density", &source_data_field);
  target_state_wrapper.mesh_get_data<double>(Portage::CELL, "density", &target_data_field);

  // Compute total mass on the source mesh to check conservation
  for (int c = 0; c < nsrccells; ++c) {
    min_source_val = std::min(min_source_val, source_data_field[c]);
    max_source_val = std::max(max_source_val, source_data_field[c]);
    double cellvol2 = source_mesh_wrapper.cell_volume(c);
    source_mass += source_data_field[c] * cellvol2;
  }

  // cell error computation
  for (int c = 0; c < ntarcells; ++c) {
    // skip ghost cells to avoid duplicated values
    if (target_mesh_wrapper.cell_get_type(c) == Portage::Entity_type::PARALLEL_OWNED) {
      // update field values bound
      min_target_val = std::min(min_target_val, target_data_field[c]);
      max_target_val = std::max(max_target_val, target_data_field[c]);
      // compute difference between exact and remapped value
      auto const centroid = target_mesh->cell_centroid(c);
      auto const cellvol = target_mesh_wrapper.cell_volume(c);
      auto const error = exact_value(centroid) - target_data_field[c];
      // update L^p norm error and target mass
      L1_error += std::abs(error) * cellvol;
      L2_error += error * error * cellvol;
      total_volume += cellvol;
      target_mass += target_data_field[c] * cellvol;
    }
  }

  // accumulate all local value on master rank
  //MPI_Barrier(comm);
  MPI_Reduce(&L1_error, global_error+0, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&L2_error, global_error+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&min_source_val, source_extents, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&max_source_val, source_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&min_target_val, target_extents, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&max_target_val, target_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&source_mass, total_mass, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&target_mass, total_mass+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  // update local value then
  L1_error = global_error[0];
  L2_error = sqrt(global_error[1]);
  min_source_val = source_extents[0];
  max_source_val = source_extents[1];
  min_target_val = target_extents[0];
  max_target_val = target_extents[1];
  source_mass = total_mass[0];
  target_mass = total_mass[1];

  if (rank == 0) {
    std::printf("done\n");
    std::printf(" \u2022 L1-norm error     = %lf\n", L1_error);
    std::printf(" \u2022 L2-norm error     = %lf\n", L2_error);
    std::printf(" \u2022 source values     = [%.15f, %.15f]\n", min_source_val, max_source_val);
    std::printf(" \u2022 target values     = [%.15f, %.15f]\n", min_target_val, max_target_val);
    std::printf(" \u2022 source total mass = %.15f\n", source_mass);
    std::printf(" \u2022 target total mass = %.15f\n", target_mass);
    std::printf(" \u2022 mass discrepancy  = %.15f\n", std::abs(source_mass - target_mass));
  }

#if DEBUG_PRINT
  Wonton::Point<dim> centroid;

  std::cout << std::endl;
  std::cout << "\tsource owned cell global id's on rank " << rank << ":" << std::endl;

  for (int ic = 0; ic < source_mesh_wrapper.num_owned_cells(); ic++) {
    source_mesh_wrapper.cell_centroid(ic, &centroid);
    int const id = source_mesh_wrapper.get_global_id(ic, Wonton::Entity_kind::CELL);
    std::cout << "source cell: " << id << ", centroid: " << centroid << std::endl;
  }

  std::cout << std::endl;
  std::cout << "\ttarget owned cell global id's on rank " << rank << ":" << std::endl;

  for (int ic = 0; ic < target_mesh_wrapper.num_owned_cells(); ic++) {
    target_mesh_wrapper.cell_centroid(ic, &centroid);
    int const id = target_mesh_wrapper.get_global_id(ic, Wonton::Entity_kind::CELL);
    std::cout << "target cell: " << id << ", centroid: " << centroid << std::endl;
  }
#endif

} //remap

/**
 * @brief Move the given target mesh points.
 *
 * It moves the target points to obtain a target mesh with same
 * connectivity but different point positions: loop over all the
 * boundary vertices assuming that we are only dealing with internally
 * generated meshes.
 * x_new = x_old + x_velocity(tcur) * deltaT
 * y_new = y_old + y_velocity(tcur) * deltaT
 *
 * @param mesh    current mesh pointer
 * @param iter    current timestep
 * @param deltaT  displacement step
 * @param periodT displacement period
 * @param scale   scaling factor
 */
template<int dim>
void move_target_mesh_nodes(std::shared_ptr<Jali::Mesh> mesh,
                            int iter, double& deltaT, double& periodT, int& scale) {

  if (rank == 0)
    std::cout << "Move target mesh ... " << std::flush;

  const int nb_nodes = mesh->num_entities(Jali::Entity_kind::NODE,
                                          Jali::Entity_type::ALL);
  double tcur = iter * deltaT;

  auto boundary_node = identify_exterior_boundary_nodes(mesh);

  for (int i = 0; i < nb_nodes; i++) {
    std::array<double, dim> coord{};
    std::array<double, dim> veloc{};
    mesh->node_get_coordinates(i, &coord);

    if (not boundary_node[i]) {
      compute_single_vortex_velocity(coord.data(), tcur, periodT, veloc.data());
      coord[0] = coord[0] + veloc[0] * deltaT / scale;
      coord[1] = coord[1] + veloc[1] * deltaT / scale;
      mesh->node_set_coordinates(i, coord.data());
    }
  }
  if (rank == 0)
    std::cout << "done" << std::endl;
}

/**
 * @brief Compute a single vortex velocity value at the given point.
 *
 * @param[in]  coords  coordinates of the current point.
 * @param[in]  tcur    scaled displacement step
 * @param[in]  periodT displacement period.
 * @param[out] veloc   computed velocity for current point.
 */
void compute_single_vortex_velocity(double* coords, double& tcur,
                                    double& periodT, double* veloc) {
  double const& x = coords[0];
  double const& y = coords[1];
  veloc[0] = -2*sin(M_PI*tcur/periodT)*sin(M_PI*x)*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*y);
  veloc[1] =  2*sin(M_PI*tcur/periodT)*sin(M_PI*x)*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*y);
}

/**
 * @brief Identify boundary nodes of the given mesh.
 *
 * @param mesh the current mesh pointer.
 * @return a flag array to check if a given mesh node is a boundary one.
 */
std::vector<bool> identify_exterior_boundary_nodes(std::shared_ptr<Jali::Mesh> mesh) {

  int nfaces = mesh->num_entities(Jali::Entity_kind::FACE, Jali::Entity_type::ALL);
  int nnodes = mesh->num_entities(Jali::Entity_kind::NODE, Jali::Entity_type::ALL);

  std::vector<bool> boundary_node(nnodes, false);

  for (int f = 0; f < nfaces; f++) {
    std::vector<int> fcells;
    mesh->face_get_cells(f, Jali::Entity_type::ALL, &fcells);

    if (fcells.size() == 1) {
      // if the face is on an exterior boundary, all its nodes are too
      std::vector<int> fnodes;
      mesh->face_get_nodes(f, &fnodes);
      for (auto const &n : fnodes)
        boundary_node[n] = true;
    }
  }
  return boundary_node;
}

/**
 * @brief Print an error message followed by command-line usage, and exits.
 *
 * @param message the error message to be displayed
 * @return status code
 */
int abort(std::string message) {
  if (rank == 0) {
    std::fprintf(stderr, "\e[31mError: %s\e[0m\n", message.data());
    print_usage();
  }
  MPI_Finalize();
  return EXIT_FAILURE;
}

/**
 * @brief Print some infos to the user.
 *
 * @param source_mesh      source mesh pointer
 * @param target_mesh      target mesh pointer
 * @param field_expression expression of the field to remap
 * @param interp_order     order of interpolation
 * @param intersect_based  use intersection-based method
 * @param limiter          gradient limiter for internal cells
 * @param bnd_limiter      gradient limiter for boundary cells
 */
void print_infos(std::shared_ptr<Jali::Mesh> source_mesh,
                 std::shared_ptr<Jali::Mesh> target_mesh,
                 std::string field_expression,
                 int interp_order, bool intersect_based,
                 Portage::Limiter_type limiter,
                 Portage::Boundary_Limiter_type bnd_limiter) {

  // get actual number of cells on both meshes
  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  int total_count[] = {0, 0};
  int nb_source_cells = source_mesh_wrapper.num_owned_cells();
  int nb_target_cells = target_mesh_wrapper.num_owned_cells();

  MPI_Reduce(&nb_source_cells, total_count+0, 1, MPI_INT, MPI_SUM, 0, comm);
  MPI_Reduce(&nb_target_cells, total_count+1, 1, MPI_INT, MPI_SUM, 0, comm);

  if (rank == 0) {
    nb_source_cells = total_count[0];
    nb_target_cells = total_count[1];
    std::cout << " \u2022 source mesh has " << nb_source_cells << " cells" << std::endl;
    std::cout << " \u2022 target mesh has " << nb_target_cells << " cells" << std::endl;
    std::cout << " \u2022 material field: \""<< field_expression << "\""<< std::endl;
    std::cout << " \u2022 interpolation order: " << interp_order << std::endl;
    if (interp_order == 2) {
      std::cout << " \u2022 internal gradient limiter: " << to_string(limiter) << std::endl;
      std::cout << " \u2022 boundary gradient limiter: " << to_string(bnd_limiter) << std::endl;
    }
    std::cout << " \u2022 method: "<< (intersect_based ? "\"intersect-based\"" : "\"swept-face\"");
    std::cout << std::endl << std::endl;
  }
}