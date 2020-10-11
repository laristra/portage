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

// wonton includes
#include "Mesh.hh"         // see https://github.com/lanl/jali
#include "MeshFactory.hh"
#include "JaliState.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// portage includes
#include "portage/driver/coredriver.h"
#include "portage/search/search_kdtree.h"
#include "portage/search/search_swept_face.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/support/portage.h"
#include "portage/support/mpi_collate.h"
#include "portage/support/timer.h"
#include "user_field.h" // parsing and evaluating user defined expressions

#ifdef PORTAGE_HAS_TANGRAM
  #include "tangram/intersect/split_r2d.h"
  #include "tangram/intersect/split_r3d.h"
  #include "tangram/reconstruct/MOF.h"
  #include "tangram/reconstruct/VOF.h"
#endif

#include <cctype>
#include <string>


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

static int rank = 0;
static int nb_ranks = 1;
static MPI_Comm comm = MPI_COMM_WORLD;

/**
 * @brief Move the given point.
 *
 * 2D: the point is updated using a periodic vortex velocity:
 *   t: current timestep
 *   dT = T/n with n: number of timesteps
 *   u = -2*sin(PI*t/T) * sin(PI*x)^2 * sin(PI*y)*cos(PI*y);
 *   v =  2*sin(PI*t/T) * sin(PI*x)*cos(PI*x) * sin(PI*y)^2;
 *   x += x_old + u*t * dT
 *   y += y_old + v*t * dT
 *
 * 3D: the point is updated as follow:
 *   t: current timestep
 *   x += 0.1 * sin(2*PI*x/(n-t))
 *   y += 0.1 * sin(2*PI*y/(n-t))
 *   z += 0.1 * sin(2*PI*z/(n-t))
 *
 * @param coords     coordinates of the given point
 * @param iter       current iteration
 * @param ntimesteps number of timesteps
 * @param scale      displacement scaling factor
 */
template<int dim>
void rotate_vortex(double* coords, int iter, int ntimesteps, int scale=1);

template<int dim>
void shift_point(double* coords, double delta);

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
 * @param simple  use simple displacement scheme
 */
template<int dim>
void move_points(std::shared_ptr<Jali::Mesh> mesh,
                 double p_min, double p_max,
                 int iter, int ntimesteps, int scale, bool simple);

/**
 * @brief Update the source mesh to the previous target mesh.
 *
 * @param source_mesh: the source mesh to update.
 * @param target_mesh: the target mesh to copy from.
 */
template<int dim>
void update_mesh(std::shared_ptr<Jali::Mesh> source_mesh,
                 std::shared_ptr<Jali::Mesh> target_mesh);

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
 * @param result_file   field file name for imported meshes.
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
           bool intersect_based, bool mesh_output, bool dump_field,
           std::string result_file, int iteration,
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
                 Portage::Boundary_Limiter_type bnd_limiter,
                 int ntimesteps, bool keep_source, bool simple_move);

/**
 * @brief Print command-line usage.
 *
 */
void print_usage() {

  std::cout << "Usage: swept_face_demo [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "\t--help                  show this help message and exit"    << std::endl;
  std::cout << "\t--dim           INT     dimension of the problem"           << std::endl;
  std::cout << "\t--ncells        INT     cells per axis for generated grids" << std::endl;
  std::cout << "\t--source_file   STRING  source mesh file to import"         << std::endl;
  std::cout << "\t--target_file   STRING  target mesh file to import"         << std::endl;
  std::cout << "\t--remap_order   INT     order of interpolation"             << std::endl;
  std::cout << "\t--intersect     BOOL    use intersection-based remap"       << std::endl;
  std::cout << "\t--field         STRING  numerical field to remap"           << std::endl;
  std::cout << "\t--ntimesteps    INT     number of timesteps"                << std::endl;
  std::cout << "\t--keep_source   BOOL    keep source mesh for all timesteps" << std::endl;
  std::cout << "\t--simple        BOOL    use simple point displacement"      << std::endl;
  std::cout << "\t--scale_by      FLOAT   displacement scaling factor"        << std::endl;
  std::cout << "\t--limiter       STRING  gradient limiter to use"            << std::endl;
  std::cout << "\t--bnd_limiter   STRING  gradient limiter for boundary"      << std::endl;
  std::cout << "\t--output_meshes BOOL    dump meshes"                        << std::endl;
  std::cout << "\t--result_file   STRING  dump remapped field to file"        << std::endl;
#if ENABLE_TIMINGS
  std::cout << "\t--scaling       STRING  scaling study [strong|weak]"        << std::endl;
  std::cout << "\t--only_threads  BOOL    thread scaling profiling"           << std::endl;
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

  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nb_ranks);
  MPI_Comm_rank(comm, &rank);

  if (argc == 1)
    return abort("");

  int dim = 2;
  int ncells = 0;
  int interp_order = 1;
  int ntimesteps = 4;
  int scale = 10;
  bool mesh_output = false;
  bool intersect_based = false;
  bool update_source = true;
  bool simple_shift = false;
  // point coordinates bounds
  double const p_min = 0.0;
  double const p_max = 1.0;

  std::string source_file;
  std::string target_file;
  std::string field_expression;
  std::string result_file;

  auto limiter     = Portage::Limiter_type::NOLIMITER;
  auto bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;

#if ENABLE_TIMINGS
  bool only_threads = false;
  std::string scaling_type = "strong";
#endif

  // convert to lower case and check if true
  auto is_true = [](std::string input) -> bool {
    std::transform(input.begin(), input.end(), input.begin(),
                   [](auto c) { return std::tolower(c); });
    return input == "true" or input == "y" or input == "yes";
  };

  // parse the input
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of('=');
    std::string keyword = arg.substr(keyword_beg, keyword_end-keyword_beg);
    std::string valueword = arg.substr(keyword_end+1, len-(keyword_end+1));

    if (keyword == "dim") {
      dim = stoi(valueword);
      assert(dim == 2 || dim == 3);
    } else if (keyword == "ncells") {
      ncells = stoi(valueword);
      assert(ncells > 0);
    } else if (keyword == "source_file") {
      source_file = valueword;
    } else if (keyword == "target_file") {
      target_file = valueword;
    } else if (keyword == "scale_by") {
      scale = std::stoi(valueword);
      assert(scale > 0);
    } else if (keyword == "remap_order") {
      interp_order = std::stoi(valueword);
      assert(interp_order > 0 and interp_order < 3);
    } else if (keyword == "intersect") {
      intersect_based = is_true(valueword);
    } else if (keyword == "field") {
      field_expression = valueword;
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
      assert(ntimesteps > 0);
    } else if (keyword == "keep_source") {
      update_source = not is_true(valueword);
    } else if (keyword == "simple") {
      simple_shift = is_true(valueword);
    } else if (keyword == "output_meshes") {
      mesh_output = is_true(valueword);
    } else if (keyword == "result_file") {
      result_file = valueword;
#if ENABLE_TIMINGS
    } else if (keyword == "only_threads"){
      only_threads = (nb_ranks == 1 and is_true(valueword));
    } else if (keyword == "scaling") {
      scaling_type = valueword;
      assert(valueword == "strong" or valueword == "weak");
#endif
    } else if (keyword == "help") {
      return abort("");
    } else
      return abort("unrecognized option "+ keyword);
  }

  // some input error checking
  bool invalid_source = (ncells and not source_file.empty()) or (not ncells and source_file.empty());
  bool invalid_target = (ncells and not target_file.empty()) or (not ncells and target_file.empty());
  bool no_field_expression = ncells and field_expression.empty();
  bool const generated_grids = ncells > 0;

  if (invalid_source)
    return abort("please choose either generated or imported source mesh");

  if (invalid_target)
    return abort("please choose either generated or imported target mesh");

  if (no_field_expression)
    return abort("no field imposed on internally generated source mesh");

  if (not generated_grids and ntimesteps > 1)
    return abort("multiple timesteps only valid for generated grids");

#if ENABLE_TIMINGS
  auto profiler = std::make_shared<Profiler>();
  // save params for after
  profiler->params.ranks   = nb_ranks;
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
    std::cout << (generated_grids ? "Generate" : "Import") << " meshes ... " << std::flush;
  }

  std::shared_ptr<Jali::Mesh> source_mesh;
  std::shared_ptr<Jali::Mesh> target_mesh;

  Jali::MeshFactory factory(comm);
  factory.included_entities(Jali::Entity_kind::ALL_KIND);

  if (not generated_grids) {
    // import meshes while using a graph partitioner to partition them.
    factory.partitioner(Jali::Partitioner_type::METIS);
    source_mesh = factory(source_file);
    target_mesh = factory(target_file);
  } else {
    // generate distributed source and target meshes
    // using a block partitioner to partition them.
    // (!) no boundary mismatch
    factory.partitioner(Jali::Partitioner_type::BLOCK);
    if (dim == 2) {
      source_mesh = factory(p_min, p_min, p_max, p_max, ncells, ncells);
      target_mesh = factory(p_min, p_min, p_max, p_max, ncells, ncells);
    } else if (dim == 3) {
      source_mesh = factory(p_min, p_min, p_min, p_max, p_max, p_max, ncells, ncells, ncells);
      target_mesh = factory(p_min, p_min, p_min, p_max, p_max, p_max, ncells, ncells, ncells);
    }
  }

  // make sure that we have the right dimension and that source and
  // target mesh dimensions matches.
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  assert(unsigned(dim) == source_mesh->space_dimension());

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
              intersect_based, limiter, bnd_limiter,
              ntimesteps, not update_source, simple_shift);

  if (not generated_grids) {
    switch (dim) {
      case 2: remap<2>(source_mesh, target_mesh, field_expression, interp_order,
                       limiter, bnd_limiter, intersect_based, mesh_output, true,
                       result_file, 1, profiler); break;
      case 3: remap<3>(source_mesh, target_mesh, field_expression, interp_order,
                       limiter, bnd_limiter, intersect_based, mesh_output, true,
                       result_file, 1, profiler); break;
      default: break;
    }
  } else {
    for (int i = 1; i <= ntimesteps; i++) {
      if (rank == 0) {
        std::cout << std::endl;
        std::cout << " ------------- ";
        std::cout << "timestep "<< i;
        std::cout << " ------------- ";
        std::cout << std::endl;
      }

      switch (dim) {
        case 2:
          // update source mesh if necessary
          if (update_source and i > 1)
            update_mesh<2>(source_mesh, target_mesh);
          // move target points for internally generated grid
          move_points<2>(target_mesh, p_min, p_max, i, ntimesteps, scale, simple_shift);
          // perform actual remap and output related errors
          remap<2>(source_mesh, target_mesh, field_expression, interp_order,
                   limiter, bnd_limiter, intersect_based, mesh_output, i == ntimesteps,
                   result_file, i, profiler);
          break;
        case 3:
          // update source mesh if necessary
          if (update_source and i > 1)
            update_mesh<3>(source_mesh, target_mesh);
          // move target points for internally generated grid
          move_points<3>(target_mesh, p_min, p_max, i, ntimesteps, scale, simple_shift);
          // perform actual remap and output related errors
          remap<3>(source_mesh, target_mesh, field_expression, interp_order,
                   limiter, bnd_limiter, intersect_based, mesh_output,i == ntimesteps,
                   result_file, i, profiler);
          break;
        default: break;
      }
    }
  }



#if ENABLE_TIMINGS
  profiler->time.total = timer::elapsed(start);

  // dump timing data
  if (rank == 0)
    profiler->dump();
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
 * @param result_file   field file name for imported meshes.
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
           bool intersect_based, bool mesh_output, bool dump_field,
           std::string result_file, int iteration,
           std::shared_ptr<Profiler> profiler) {

  if (rank == 0)
    std::cout << "Remap field ... " << std::flush;

  // parallel executor
  Wonton::MPIExecutor_type mpi_executor(comm);
  Wonton::Executor_type* executor = (nb_ranks > 1 ? &mpi_executor : nullptr);

  // create native jali state managers for source and target
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  // create wrappers to access underlying mesh data structures and fields
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  int const nb_source_cells = source_mesh_wrapper.num_owned_cells() +
                              source_mesh_wrapper.num_ghost_cells();
  int const nb_target_cells = target_mesh_wrapper.num_owned_cells();

  // initialize the source field and add to state
  user_field_t exact_value;
  if (not exact_value.initialize(dim, field_expression))
    throw std::runtime_error("expression parsing failure");

  double source_field[nb_source_cells];

  for (int c = 0; c < nb_source_cells; ++c) {
    Wonton::Point<dim> centroid;
    source_mesh_wrapper.cell_centroid(c, &centroid);
    source_field[c] = exact_value(centroid);
  }

  source_state_wrapper.mesh_add_data(Portage::CELL, "density", source_field);
  target_state_wrapper.mesh_add_data(Portage::CELL, "density", 0.0);
  MPI_Barrier(comm);

  // the remapper to use
  if (dim == 2) { //2D
#ifndef PORTAGE_HAS_TANGRAM
    Portage::CoreDriver<2,
                        Wonton::Entity_kind::CELL,
                        Wonton::Jali_Mesh_Wrapper,
                        Wonton::Jali_State_Wrapper>
      remapper(source_mesh_wrapper, source_state_wrapper,
               target_mesh_wrapper, target_state_wrapper, executor);
#else
    Portage::CoreDriver<2,
                        Wonton::Entity_kind::CELL,
                        Wonton::Jali_Mesh_Wrapper,
                        Wonton::Jali_State_Wrapper,
                        Wonton::Jali_Mesh_Wrapper,
                        Wonton::Jali_State_Wrapper,
                        Tangram::MOF,
                        Tangram::SplitRnD<2>,
                        Tangram::ClipRnD<2>>
      remapper(source_mesh_wrapper, source_state_wrapper,
               target_mesh_wrapper, target_state_wrapper, executor);
#endif

#if ENABLE_TIMINGS
    auto tic = timer::now();
#endif

    // search for neighboring cells candidates
    auto candidates = (intersect_based ? remapper.template search<Portage::SearchKDTree>()
                                       : remapper.template search<Portage::SearchSweptFace>());
#if ENABLE_TIMINGS
    profiler->time.search += timer::elapsed(tic, true);
#endif

    // compute interpolation weights using intersection or swept regions moments.
    auto weights = (intersect_based
      ? remapper.template intersect_meshes<Portage::IntersectR2D>(candidates)
      : remapper.template intersect_meshes<Portage::IntersectSweptFace2D>(candidates));

#if ENABLE_TIMINGS
    profiler->time.intersect += timer::elapsed(tic, true);
#endif

    // interpolate the field eventually
    if (interp_order == 1) {
      remapper.template interpolate_mesh_var<double, Portage::Interpolate_1stOrder>("density",
                                                                                    "density",
                                                                                    weights);
    } else if (interp_order == 2) {
      auto gradients = remapper.compute_source_gradient("density", limiter, bnd_limiter);

#if ENABLE_TIMINGS
      profiler->time.gradient += timer::elapsed(tic);
#endif

     remapper.template interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>("density",
                                                                                    "density",
                                                                                    weights, &gradients);
										    
    }

#if ENABLE_TIMINGS
     profiler->time.interpolate += timer::elapsed(tic, true);
#endif

  } else { //3D
#ifndef PORTAGE_HAS_TANGRAM
    Portage::CoreDriver<3,
                        Wonton::Entity_kind::CELL,
                        Wonton::Jali_Mesh_Wrapper,
                        Wonton::Jali_State_Wrapper>
      remapper(source_mesh_wrapper, source_state_wrapper,
               target_mesh_wrapper, target_state_wrapper, executor);
#else
    Portage::CoreDriver<3,
                        Wonton::Entity_kind::CELL,
                        Wonton::Jali_Mesh_Wrapper,
                        Wonton::Jali_State_Wrapper,
                        Wonton::Jali_Mesh_Wrapper,
                        Wonton::Jali_State_Wrapper,
                        Tangram::MOF,
                        Tangram::SplitRnD<3>,
                        Tangram::ClipRnD<3>>
      remapper(source_mesh_wrapper, source_state_wrapper,
               target_mesh_wrapper, target_state_wrapper, executor);
#endif

#if ENABLE_TIMINGS
    auto tic = timer::now();
#endif

    // search for neighboring cells candidates
    auto candidates = (intersect_based ? remapper.template search<Portage::SearchKDTree>()
                                       : remapper.template search<Portage::SearchSweptFace>());

#if ENABLE_TIMINGS
    profiler->time.search += timer::elapsed(tic, true);
#endif

    // compute interpolation weights using intersection or swept regions moments.
    auto weights = (intersect_based
      ? remapper.template intersect_meshes<Portage::IntersectR3D>(candidates)
      : remapper.template intersect_meshes<Portage::IntersectSweptFace3D>(candidates));

#if ENABLE_TIMINGS
    profiler->time.intersect += timer::elapsed(tic, true);
#endif

    // interpolate the field eventually
    if (interp_order == 1) {
      remapper.template interpolate_mesh_var<double, Portage::Interpolate_1stOrder>("density",
                                                                                    "density",
                                                                                    weights);
    } else if (interp_order == 2) {
      auto gradients = remapper.compute_source_gradient("density", limiter, bnd_limiter);

#if ENABLE_TIMINGS
      profiler->time.gradient += timer::elapsed(tic);
#endif

      remapper.template interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>("density",
                                                                                    "density",
                                                                                    weights, &gradients);
    }

#if ENABLE_TIMINGS
     profiler->time.interpolate += timer::elapsed(tic, true);
#endif

  }

#if ENABLE_TIMINGS
  profiler->time.remap = profiler->time.search + 
                         profiler->time.intersect + 
                         profiler->time.interpolate ;

  if (rank == 0)
    std::printf("done. \e[32m(%.3f s)\e[0m\n", profiler->time.remap);
#else
  if (rank == 0)
    std::cout << "done" << std::endl;
#endif

  if (rank == 0)
    std::cout << "Extract stats ... " << std::flush;

  double const* target_field;
  double min_source_val   = std::numeric_limits<double>::max();
  double max_source_val   = std::numeric_limits<double>::min();
  double min_target_val   = min_source_val;
  double max_target_val   = max_source_val;
  double source_extents[] = { min_source_val, max_source_val };
  double target_extents[] = { min_target_val, max_target_val };
  double target_mass      = 0.0;
  double source_mass      = 0.0;
  double total_mass[]     = {0.0, 0.0};
  double global_error[]   = {0.0, 0.0};
  double total_volume     = 0.0;
  double L1_error         = 0.0;
  double L2_error         = 0.0;

  target_state_wrapper.mesh_get_data<double>(Portage::CELL, "density", &target_field);

  // compute total mass on the source mesh to check conservation
  for (int c = 0; c < nb_source_cells; ++c) {
    if (source_mesh_wrapper.cell_get_type(c) == Portage::Entity_type::PARALLEL_OWNED) {
      min_source_val = std::min(min_source_val, source_field[c]);
      max_source_val = std::max(max_source_val, source_field[c]);
      double const volume = source_mesh_wrapper.cell_volume(c);
      source_mass += source_field[c] * volume;
    }
  }

  for (int c = 0; c < nb_target_cells; ++c) {
    // skip ghost cells to avoid duplicated values
    if (target_mesh_wrapper.cell_get_type(c) == Portage::Entity_type::PARALLEL_OWNED) {
      // update field values bound
      min_target_val = std::min(min_target_val, target_field[c]);
      max_target_val = std::max(max_target_val, target_field[c]);
      // compute difference between exact and remapped value
      Wonton::Point<dim> centroid;
      target_mesh_wrapper.cell_centroid(c, &centroid);
      double const volume = target_mesh_wrapper.cell_volume(c);
      double const error = exact_value(centroid) - target_field[c];
      // update L^p norm error and target mass
      L1_error += std::abs(error) * volume;
      L2_error += error * error * volume;
      total_volume += volume;
      target_mass += target_field[c] * volume;
    }
  }

  // accumulate all local value on master rank
  MPI_Barrier(comm);
  MPI_Reduce(&L1_error, global_error+0, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&L2_error, global_error+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&min_source_val, source_extents+0, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&max_source_val, source_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&min_target_val, target_extents+0, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&max_target_val, target_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&source_mass, total_mass+0, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&target_mass, total_mass+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  if (rank == 0) {
    // update local value then
    L1_error = global_error[0];
    L2_error = sqrt(global_error[1]);
    min_source_val = source_extents[0];
    max_source_val = source_extents[1];
    min_target_val = target_extents[0];
    max_target_val = target_extents[1];
    source_mass = total_mass[0];
    target_mass = total_mass[1];

    std::printf("done\n");
    std::printf(" \u2022 L1-norm error     = %.15f\n", L1_error);
    std::printf(" \u2022 L2-norm error     = %.15f\n", L2_error);
    std::printf(" \u2022 source values     = [%.15f, %.15f]\n", min_source_val, max_source_val);
    std::printf(" \u2022 target values     = [%.15f, %.15f]\n", min_target_val, max_target_val);
    std::printf(" \u2022 source total mass = %.15f\n", source_mass);
    std::printf(" \u2022 target total mass = %.15f\n", target_mass);
    std::printf(" \u2022 mass discrepancy  = %.15f\n", std::abs(source_mass - target_mass));
  }

  // dump meshes if requested
  if (mesh_output) {
    if (rank == 0)
      std::cout << "Dump meshes ... " << std::flush;

    std::string suffix = "time_" + std::to_string(iteration)+ ".exo";
    source_state->export_to_mesh();
    target_state->export_to_mesh();
    source_mesh->write_to_exodus_file("source_" + suffix);
    target_mesh->write_to_exodus_file("target_" + suffix);

    if (rank == 0)
      std::cout << "done." << std::endl;
  }

  // dump remapped field
  if (dump_field and not result_file.empty()) {

    if (rank == 0)
      std::cout << "Dump field ... " << std::flush;

    std::vector<int> index(nb_target_cells);
    std::vector<double> value(nb_target_cells);

    for (int i=0; i < nb_target_cells; i++) {
      index[i] = target_mesh->GID(i, Jali::Entity_kind::CELL);
      value[i] = target_field[i];
    }

    // sort the field values by global ID
    std::vector<int> swap;
    Portage::argsort(index, swap);   // find sorting indices based on global IDS
    Portage::reorder(index, swap);   // sort the global ids
    Portage::reorder(value, swap);  // sort the values

    std::string result_dist = result_file +"_dist.txt";
    result_file += ".txt";

    if (nb_ranks > 1) {
      int width = static_cast<int>(std::ceil(std::log10(nb_ranks)));
      char ext[10];
      std::snprintf(ext, sizeof(ext), "%0*d", width, rank);
      result_file = result_file + "." + std::string(ext);
      result_dist = result_dist + "." + std::string(ext);
    }

    // write out the values
    std::ofstream file(result_file);
    file << std::scientific;
    file.precision(17);

    // for compare app
    for (int i = 0; i < nb_target_cells; i++)
      file << index[i] <<"\t"<< value[i] << std::endl;

    file.close();
    file.clear();
    file.open(result_dist);
    file << std::scientific;
    file.precision(17);

    for (int i = 0; i < nb_target_cells; i++)
      file << index[i] <<"\t"<< 0 <<"\t"<< value[i] << std::endl;

    file.close();

    MPI_Barrier(comm);
    if (rank == 0)
      std::cout << "done." << std::endl;
  }

  MPI_Barrier(comm);
} //remap

/**
 * @brief Move the given target mesh points.
 *
 * It moves points to obtain a target mesh with same connectivity
 * but with different point positions: loop over all the boundary
 * nodes assuming that we are only dealing with internally
 * generated grids.
 *
 * @tparam dim       dimension of the problem
 * @param mesh       current mesh pointer
 * @param p_min      point coordinates lower bound
 * @param p_max      point coordinates upper bound
 * @param iter       current iteration
 * @param ntimesteps number of timesteps
 * @param scale      scaling factor
 * @param simple     simple displacement scheme
 */
template<int dim>
void move_points(std::shared_ptr<Jali::Mesh> mesh,
                 double p_min, double p_max,
                 int iter, int ntimesteps, int scale, bool simple) {
  if (rank == 0)
    std::cout << "Move target mesh ... " << std::flush;

  int const nb_nodes = mesh->num_nodes<Jali::Entity_type::ALL>();
  int const nb_cells = mesh->num_nodes<Jali::Entity_type::ALL>();
  double const epsilon = 1.E-16;

  // --------------------------------------
  // check if the given point should be skipped or not.
  auto skip = [&](auto const& p) -> bool {
    // check if two real numbers are equal
    auto equals = [&](double const& u, double const& v) -> bool {
      return std::abs(u - v) < epsilon;
    };

    double const& x = p[0];
    double const& y = p[1];
    if (dim == 3) {
      double const& z = p[2];
      // skip only corners
      return (equals(x, p_min) and equals(y, p_min) and equals(z, p_min))
          or (equals(x, p_min) and equals(y, p_min) and equals(z, p_max))
          or (equals(x, p_min) and equals(y, p_max) and equals(z, p_min))
          or (equals(x, p_min) and equals(y, p_max) and equals(z, p_max))
          or (equals(x, p_max) and equals(y, p_min) and equals(z, p_min))
          or (equals(x, p_max) and equals(y, p_min) and equals(z, p_max))
          or (equals(x, p_max) and equals(y, p_max) and equals(z, p_min))
          or (equals(x, p_max) and equals(y, p_max) and equals(z, p_max));
    } else {
      // skip boundary points
      return equals(x, p_min) or equals(x, p_max)
          or equals(y, p_min) or equals(y, p_max);
    }
  };
  // --------------------------------------
  if (not simple) {
    for (int i = 0; i < nb_nodes; i++) {
      std::array<double, dim> point{};
      mesh->node_get_coordinates(i, &point);
      if (not skip(point)) {
        rotate_vortex<dim>(point.data(), iter, ntimesteps, scale);
        mesh->node_set_coordinates(i, point.data());
      }
    }
  } else {
    double const delta_shift = (p_max - p_min) / (scale * nb_cells);

    for (int i = 0; i < nb_nodes; i++) {
      std::array<double, dim> point{};
      mesh->node_get_coordinates(i, &point);
      if (not skip(point)) {
        point[0] += delta_shift;
        mesh->node_set_coordinates(i, point.data());
      }
    }
  }

  MPI_Barrier(comm);
  if (rank == 0)
    std::cout << "done" << std::endl;
}

/**
 * @brief Update the source mesh to the previous target mesh.
 *
 * @param source_mesh: the source mesh to update.
 * @param target_mesh: the target mesh to copy from.
 */
template<int dim>
void update_mesh(std::shared_ptr<Jali::Mesh> source_mesh,
                 std::shared_ptr<Jali::Mesh> target_mesh) {
  if (rank == 0)
    std::cout << "Update source mesh ... " << std::flush;

  int const nb_nodes = source_mesh->num_nodes<Jali::Entity_type::ALL>();

  for (int i = 0; i < nb_nodes; i++) {
    std::array<double, dim> point {}, p_old {}, p_new {};
    target_mesh->node_get_coordinates(i, &point);
    source_mesh->node_get_coordinates(i, &p_old);
    source_mesh->node_set_coordinates(i, point.data());
    source_mesh->node_get_coordinates(i, &p_new);
  }

  MPI_Barrier(comm);
  if (rank == 0)
    std::cout << "done" << std::endl;
}

/**
 * @brief Move the two-dimensional point.
 * Its coordinates is updated using a periodic vortex velocity:
 *   t: current timestep
 *   dT = T/n with n: number of timesteps
 *   u = -2*sin(PI*t/T) * sin(PI*x)^2 * sin(PI*y) * cos(PI*y);
 *   v =  2*sin(PI*t/T) * sin(PI*x) * cos(PI*x) * sin(PI*y)^2;
 *   x += u * t * dT
 *   y += v * t * dT
 *
 * @param coords     coordinates of the given point
 * @param iter       current iteration
 * @param ntimesteps number of timesteps
 * @param scale      displacement scaling factor
 */
template<>
void rotate_vortex<2>(double* coords, int iter, int ntimesteps, int scale) {

  double const periodT = 2.0;
  double const deltaT = (periodT/ntimesteps)/scale;
  double const tcur = iter * deltaT;
  double const step = deltaT;

  double& x = coords[0];
  double& y = coords[1];
  double vx = -2*sin(M_PI*tcur/periodT)*sin(M_PI*x)*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*y);
  double vy =  2*sin(M_PI*tcur/periodT)*sin(M_PI*x)*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*y);

  x += vx * step;
  y += vy * step;
}

/**
 * @brief Move the three-dimensional point.
 * Its coordinate is updated as follow:
 *   t: current timestep
 *   dT = T/n with n: number of timesteps
 *   u = -2*sin(PI*t/T) * sin(PI*x)^2 * sin(PI*y) * cos(PI*y);
 *   v =  2*sin(PI*t/T) * sin(PI*x) * cos(PI*x) * sin(PI*y)^2;
 *   w = -2*sin(PI*t/T) * sin(PI*x*y) * cos(PI*x*y) * sin(PI*z)^2;
 *   x += u * t * dT
 *   y += v * t * dT
 *   z += w * t * dT
 *
 * @param coords     coordinates of the given point
 * @param iter       current iteration (unused)
 * @param ntimesteps number of timesteps (unused)
 * @param scale      displacement scaling factor (unused)
 */
template<>
void rotate_vortex<3>(double* coords, int iter, int ntimesteps, int scale) {

  double const periodT = 2.0;
  double const deltaT = (periodT/ntimesteps)/scale;
  double const tcur = iter * deltaT;
  double const step = deltaT;

  double& x = coords[0];
  double& y = coords[1];
  double& z = coords[2];
  double vx = -2*sin(M_PI*tcur/periodT) * pow(sin(M_PI*x),2) * sin(M_PI*y)*cos(M_PI*y);
  double vy =  2*sin(M_PI*tcur/periodT) * sin(M_PI*x)*cos(M_PI*x) * pow(sin(M_PI*y),2);
  double vz = -2*sin(M_PI*tcur/periodT) * sin(M_PI*x*y)*cos(M_PI*x*y) * pow(sin(M_PI*z),2);

  x += vx * step;
  y += vy * step;
  z += vz * step;
//  static double const factor = 0.075;
//  coords[0] += factor * sin(2 * M_PI * coords[0]);
//  coords[1] += factor * sin(2 * M_PI * coords[1]);
//  coords[2] += factor * sin(2 * M_PI * coords[2]);
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
                 Portage::Boundary_Limiter_type bnd_limiter,
                 int ntimesteps, bool keep_source, bool simple_move) {

  // get actual number of cells on both meshes
  int total_count[] = {0, 0};
  int nb_source_cells = Wonton::Jali_Mesh_Wrapper(*source_mesh).num_owned_cells();
  int nb_target_cells = Wonton::Jali_Mesh_Wrapper(*target_mesh).num_owned_cells();
  MPI_Reduce(&nb_source_cells, total_count+0, 1, MPI_INT, MPI_SUM, 0, comm);
  MPI_Reduce(&nb_target_cells, total_count+1, 1, MPI_INT, MPI_SUM, 0, comm);

  if (rank == 0) {
    assert(nb_source_cells == nb_target_cells);

    std::string limiter_name;
    std::string bnd_limiter_name;

    switch(limiter) {
      case Portage::NOLIMITER: limiter_name = "none"; break;
      case Portage::BARTH_JESPERSEN: limiter_name = "barth-jespersen"; break;
      default: limiter_name = "undefined"; break;
    }

    switch (bnd_limiter) {
      case Portage::BND_NOLIMITER: bnd_limiter_name = "none"; break;
      case Portage::BND_ZERO_GRADIENT: bnd_limiter_name = "zero-gradient"; break;
      case Portage::BND_BARTH_JESPERSEN: bnd_limiter_name = "barth-jespersen"; break;
      default: bnd_limiter_name = "undefined"; break;
    }

    std::cout << " \u2022 source mesh has " << total_count[0] << " cells" << std::endl;
    std::cout << " \u2022 target mesh has " << total_count[1] << " cells" << std::endl;
    std::cout << " \u2022 material field: \""<< field_expression << "\""<< std::endl;
    std::cout << " \u2022 interpolation order: " << interp_order << std::endl;
    if (interp_order == 2) {
      std::cout << " \u2022 internal gradient limiter: " << limiter_name << std::endl;
      std::cout << " \u2022 boundary gradient limiter: " << bnd_limiter_name << std::endl;
    }
    std::cout << " \u2022 method: "<< (intersect_based ? "\"intersect-based\"" : "\"swept-face\"") << std::endl;
    std::cout << " \u2022 timesteps: "<< ntimesteps << std::endl;
    std::cout << " \u2022 keep source: "<< std::boolalpha << keep_source << std::endl;
    std::cout << " \u2022 displacement: "<< (simple_move ? "simple shift" : "vortex") << std::endl;
    std::cout << std::endl;
  }
}
