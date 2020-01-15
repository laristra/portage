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
#include "portage/search/search_swept_face.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/support/portage.h"
#include "portage/support/timer.h"
#include "user_field.h" // parsing and evaluating user defined expressions

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
 * for single material field between two meshes - the meshes can be
 * internally generated rectangular meshes or externally read unstructured
 * meshes.
 */

//////////////////////////////////////////////////////////////////////
// Forward declarations

void print_usage() {

  std::cout << "Usage: swept_face_demo [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "\t--help                  show this help message and exit" << std::endl;
  std::cout << "\t--dim           INT     dimension of the problem"        << std::endl;
  std::cout << "\t--ncells        INT     number of cells per axis"        << std::endl;
  std::cout << "\t--remap_order   INT     order of interpolation"          << std::endl;
  std::cout << "\t--field         STRING  numerical field to remap"        << std::endl;
  std::cout << "\t--ntimesteps    INT     number of timesteps"             << std::endl;
  std::cout << "\t--scale_by      FLOAT   displacement scaling factor"     << std::endl;
  std::cout << "\t--limiter       STRING  gradient limiter to use"         << std::endl;
  std::cout << "\t--bnd_limiter   STRING  gradient limiter for boundary"   << std::endl;
  std::cout << "\t--output_meshes CHAR    dump meshes [y|n]"               << std::endl;
#if ENABLE_TIMINGS
  std::cout << "\t--scaling       STRING  scaling study [strong|weak]"     << std::endl;
  std::cout << "\t--only_threads  CHAR    thread scaling profiling [y|n]"  << std::endl;
#endif
}

// Forward declaration of function to run remap on two meshes and
// return the L1 and L2 error norm in the remapped field w.r.t. to an
// analytically imposed field. If no field was imposed, the errors are
// returned as 0

template<int dim>
void remap(std::shared_ptr<Jali::Mesh> source_mesh,                     
           std::shared_ptr<Jali::Mesh> target_mesh,                     
           std::string field_expression,                                
           int interp_order,                                            
           Portage::Limiter_type limiter,                               
           Portage::Boundary_Limiter_type bnd_limiter,
           bool mesh_output,                                            
           std::string field_filename,                                  
           int rank, int numpe, int iteration,                          
           double& L1_error, double& L2_error,
           MPI_Comm comm,
           std::shared_ptr<Profiler> profiler);

// Forward declaration of function to find nodes on exterior boundary
void find_nodes_on_exterior_boundary(std::shared_ptr<Jali::Mesh> mesh,
                                     std::vector<bool>& node_on_bnd);

void move_target_mesh_nodes(std::shared_ptr<Jali::Mesh> mesh,
                            int i, double& deltaT, double& periodT, int& scale);

void single_vortex_velocity_function(double* coords, double& tcur,
                                     double& periodT, double* disp);

//////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  // Initialize MPI
  int numpe, rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &numpe);
  MPI_Comm_rank(comm, &rank);

  if (argc == 1) {
    print_usage();
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  int dim = 2;
  int ncells = 10;
  int interp_order = 1;
  int ntimesteps = 4;
  int scale = 20;
  bool mesh_output = false; // (!) 'write_to_gmv' segfaults in parallel

  std::string field_expression;
  std::string field_filename;  // No default;

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
    } else if (keyword == "scale_by") {
      scale = stoi(valueword);
    } else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
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
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
    }
#if ENABLE_TIMINGS
      else if (keyword == "only_threads"){
      only_threads = (numpe == 1 and valueword == "y");
    } else if (keyword == "scaling") {
      assert(valueword == "strong" or valueword == "weak");
      scaling_type = valueword;
    }
#endif
    else if (keyword == "help") {
      print_usage();
      MPI_Finalize();
      return EXIT_SUCCESS;
    } else {
      std::cerr << "Unrecognized option " << keyword << " !\n";
      MPI_Finalize();
      return EXIT_FAILURE;
    }
  }

#if ENABLE_TIMINGS
  auto profiler = std::make_shared<Profiler>();
  // save params for after
  profiler->params.ranks   = numpe;
  profiler->params.nsource = std::pow(nsourcecells, dim);
  profiler->params.ntarget = std::pow(ntargetcells, dim);
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

  if (rank == 0)
    std::cout << "Starting swept_face_demo ... " << std::endl << std::endl;

  // bounds of generated mesh in each dir
  // no mismatch between source and target meshes.
  double x_min = 0.0, x_max = 1.0;
  double y_min = 0.0, y_max = 1.0;

  std::shared_ptr<Jali::Mesh> source_mesh, target_mesh;

  // generate distributed source and target meshes
  // using a BLOCK partitioner to partition them.
  Jali::MeshFactory factory(comm);
  factory.included_entities(Jali::Entity_kind::ALL_KIND);
  factory.partitioner(Jali::Partitioner_type::BLOCK);

  // Create the source and target meshes
  if (dim == 2) {
    source_mesh = factory(x_min, y_min, x_max, y_max, ncells, ncells);
    target_mesh = factory(x_min, y_min, x_max, y_max, ncells, ncells);
  }
  else if (dim == 3) {
    double z_min = 0.0, z_max = 1.0;
    source_mesh = factory(x_min, y_min, z_min, x_max, y_max, z_max,
                          ncells, ncells, ncells);
    target_mesh = factory(x_min, y_min, z_min, x_max, y_max, z_max,
                          ncells, ncells, ncells);
  }


  // Make sure we have the right dimension and that source and
  // target mesh dimensions match (important when one or both of the
  // meshes are read in)
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  assert(dim == source_mesh->space_dimension());

#if ENABLE_TIMINGS
  profiler->time.mesh_init = timer::elapsed(tic);

  if (rank == 0) {
        std::cout << "Mesh Initialization Time: " << profiler->time.mesh_init << std::endl;
    }
#endif

  double periodT = 2.0;
  double deltaT = periodT/ntimesteps;

  std::vector<double> l1_err(ntimesteps, 0.0);
  std::vector<double> l2_err(ntimesteps, 0.0);

  for (int i = 1; i < ntimesteps; i++) {

    std::cout << "------------- timestep "<< i << " -------------" << std::endl;
    // move nodes of the target mesh
    move_target_mesh_nodes(target_mesh, i, deltaT, periodT, scale);

    // Now run the remap on the meshes and get back the L1-L2 errors
    switch (dim) {
      case 2: remap<2>(source_mesh, target_mesh, field_expression, interp_order,
                       limiter, bnd_limiter, mesh_output, field_filename,
                       rank, numpe, i, l1_err[i], l2_err[i], comm, profiler); break;
      case 3: remap<3>(source_mesh, target_mesh, field_expression, interp_order,
                       limiter, bnd_limiter, mesh_output, field_filename,
                       rank, numpe, i, l1_err[i], l2_err[i], comm, profiler); break;
      default:
        std::cerr << "Invalid dimension" << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

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
//////////////////////////////////////////////////////////////////////////////
// Run a remap between two meshes and return the L1 and L2 error norms
// with respect to the specified field.

template<int dim>
void remap(std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh,
           std::string field_expression,
           int interp_order,
           Portage::Limiter_type limiter,
           Portage::Boundary_Limiter_type bnd_limiter,
           bool mesh_output,
           std::string field_filename,
           int rank, int numpe, int iteration,
           double& L1_error, double& L2_error,
           MPI_Comm comm,
           std::shared_ptr<Profiler> profiler) {

  std::cout << "Remap field using swept-face algorithm ... "<< std::endl;

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

  // Output some information for the user
  if (rank == 0) {
    std::cout << " \u2022 source mesh has " << nsrccells << " cells" << std::endl;
    std::cout << " \u2022 target mesh has " << ntarcells << " cells" << std::endl;
    std::cout << " \u2022 single material field: "<< field_expression << std::endl;
    std::cout << " \u2022 interpolation order: " << interp_order << std::endl;
    if (interp_order == 2) {
      std::cout << " \u2022 internal slope limiter is: " << limiter << std::endl;
      std::cout << " \u2022 boundary slope limiter is: " << bnd_limiter << std::endl;
    }
  }

#if ENABLE_TIMINGS
  auto tic = timer::now();
#endif

  // Compute field data on source, and add to state manager.
  user_field_t source_field;
  if (not source_field.initialize(dim, field_expression))
    MPI_Abort(comm, EXIT_FAILURE);

  std::vector<double> source_data(nsrccells);
  std::vector<double> target_data;

  for (int c = 0; c < nsrccells; ++c) {
    source_data[c] = source_field(source_mesh->cell_centroid(c));
  }

  source_state->add("density", source_mesh,
                    Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL,
                    source_data.data());

  target_state->add<double, Jali::Mesh, Jali::UniStateVector>("density", target_mesh,
                                                              Jali::Entity_kind::CELL,
                                                              Jali::Entity_type::ALL,
                                                              0.0);

#if ENABLE_TIMINGS
  profiler->time.remap = timer::elapsed(tic, true);
#endif

  std::vector<std::string> fieldnames { "density" };

  if (rank == 0) {
    std::cout << " \u2022 registered fields: ";
    for (auto&& name: source_state_wrapper.names())
      std::cout << name << " ";
    std::cout << std::endl;
  }

  MPI_Barrier(comm);

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper, executor);

  auto candidates = remapper.template search<Portage::SearchSweptFace>();
  auto weights = remapper.template intersect_meshes<Portage::IntersectSweptFace2D>(candidates);

  if (interp_order == 1) {
    remapper.template interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      "density", "density", weights
    );
  } else if (interp_order == 2) {
    auto gradients = remapper.compute_source_gradient("density", limiter, bnd_limiter);
    remapper.template interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
      "density", "density", weights, &gradients
    );
  }

#if !ENABLE_TIMINGS
  // Output some information for the user
#else
  profiler->time.remap += timer::elapsed(tic);

  if (rank == 0) {
    std::cout << "Remap Time: " << profiler->time.remap << std::endl;
  }
#endif

  // Write out the meshes if requested
  if (mesh_output) {
    if (rank == 0)
      std::cout << "Dumping data to Exodus files ... ";

    std::string suffix = std::to_string(rank) + std::to_string(iteration)+ ".exo";
    source_state->export_to_mesh();
    target_state->export_to_mesh();
    source_mesh->write_to_exodus_file("source_" + suffix);
    target_mesh->write_to_exodus_file("target_" + suffix);

    if (rank == 0)
      std::cout << "done." << std::endl;
  }

  // Compute error
  L1_error = 0.0;
  L2_error = 0.0;

  double error;
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

  if (numpe == 1 && ntarcells < 10)
    std::cout << "celldata vector on target mesh after remapping is:" << std::endl;

  target_data.resize(ntarcells);

  // Compute total mass on the source mesh to check conservation
  for (int c = 0; c < nsrccells; ++c) {
    min_source_val = std::min(min_source_val, source_data_field[c]);
    max_source_val = std::max(max_source_val, source_data_field[c]);
    double cellvol2 = source_mesh_wrapper.cell_volume(c);
    source_mass += source_data_field[c] * cellvol2;
  }

  // cell error computation
  Portage::Point<dim> ccen;
  for (int c = 0; c < ntarcells; ++c) {
    if (target_mesh_wrapper.cell_get_type(c) == Portage::Entity_type::PARALLEL_OWNED) {

      target_mesh_wrapper.cell_centroid(c, &ccen);
      double cellvol = target_mesh_wrapper.cell_volume(c);
      target_data[c] = target_data_field[c];
      min_target_val = std::min(min_target_val, target_data_field[c]);
      max_target_val = std::max(max_target_val, target_data_field[c]);
      error = source_data_field[c] - target_data_field[c];
      L1_error += std::abs(error) * cellvol;
      L2_error += error * error * cellvol;
      total_volume += cellvol;
      target_mass += target_data_field[c] * cellvol;

      if (ntarcells < 10) {
        std::printf("rank %d\n", rank);
        std::printf("cell=%4d centroid = (%8.5lf,%8.5lf)", c, ccen[0], ccen[1]);
        std::printf("  value = %10.6lf  L2 error = %lf\n", target_data_field[c], error);
      }
    }
  }

  // accumulate all local value on master rank
  MPI_Barrier(comm);
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

  std::printf(" \u2022 L1-norm error     = %lf\n", L1_error);
  std::printf(" \u2022 L2-norm error     = %lf\n", L2_error);
  std::printf(" \u2022 source values     = [%.15f, %.15f]\n", min_source_val, max_source_val);
  std::printf(" \u2022 target values     = [%.15f, %.15f]\n", min_target_val, max_target_val);
  std::printf(" \u2022 source total mass = %.15f\n", source_mass);
  std::printf(" \u2022 target total mass = %.15f\n", target_mass);
  std::printf(" \u2022 mass discrepancy  = %.15f\n", std::abs(source_mass - target_mass));

#ifdef DEBUG_PRINT
  // debug diagnostics
  std::cout << "\n----source owned cell global id's on rank " << rank << ":\n";
  for (int ic = 0; ic < sourceMeshWrapper.num_owned_cells(); ic++) {
    Wonton::Point<dim> centroid;
    sourceMeshWrapper.cell_centroid(ic,&centroid);
    std::cout << "source cell: " << sourceMeshWrapper.get_global_id(ic, Wonton::Entity_kind::CELL) <<
      "  centroid:(" <<centroid[0]<< ", " << centroid[1] << ")\n";
  }

  std::cout << "\n----target owned cell global id's on rank " << rank << ":\n";
  for (int ic = 0; ic < targetMeshWrapper.num_owned_cells(); ic++) {
    Wonton::Point<dim> centroid;
    targetMeshWrapper.cell_centroid(ic,&centroid);
    std::cout << "target cell: "  << targetMeshWrapper.get_global_id(ic, Wonton::Entity_kind::CELL) <<
      "  centroid:(" <<centroid[0]<< ", " << centroid[1] << ")\n";
  }
#endif

} //remap

void move_target_mesh_nodes(std::shared_ptr<Jali::Mesh> mesh,
                            int iter, double& deltaT, double& periodT, int& scale)
{
  std::cout << "Moving target mesh points ... ";
  double tcur = iter * deltaT;

  // Move the target nodes to obtain a target mesh with same
  // connectivity but different point positions.  Loop over all the
  // boundary vertices Assume that we are only doing internally
  // generated meshes.
  // xnew = xold + Xvelocity(tcur)*deltaT
  // ynew = yold + Yvelocity(tcur)*deltaT

  const int ntarnodes = mesh->num_entities(Jali::Entity_kind::NODE,
                                           Jali::Entity_type::ALL);

  std::vector<bool> node_on_bnd;
  find_nodes_on_exterior_boundary(mesh, node_on_bnd);

  if (ntarnodes <= 20)
    std::cout << std::endl;

  for (int i = 0; i < ntarnodes; i++)
  {
    std::array<double, 2> coords {0.0, 0.0};
    mesh->node_get_coordinates(i, &coords);

    if (ntarnodes <= 20)
      std::cout<<"Target Node = "<<i<<" Original Coords = {"<<coords[0]
               <<", "<<coords[1]<<"}"<<std::endl;

    if (!node_on_bnd[i]) {

      std::array<double, 2> disp {0.0, 0.0};
      single_vortex_velocity_function(&coords[0], tcur, periodT, &disp[0]);
      coords[0] = coords[0] + disp[0]*deltaT/scale;
      coords[1] = coords[1] + disp[1]*deltaT/scale;

      mesh->node_set_coordinates(i, &coords[0]);

      if (ntarnodes <= 20) {
        mesh->node_get_coordinates(i, &coords);
        std::cout<<"Target Node = "<<i<<" Modified Coords = {"<< coords[0]
                 <<", "<< coords[1]<<"}"<<std::endl;
      }
    }
  }

  if (ntarnodes <= 20)
    std::cout << std::endl;

  std::cout << "done" << std::endl;

} //move_target_mesh_nodes

void single_vortex_velocity_function(double* coords, double& tcur,
                                     double& periodT, double* disp)
{
  double x = coords[0];
  double y = coords[1];
  disp[0] = -2*sin(M_PI*tcur/periodT)*sin(M_PI*x)*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*y);
  disp[1] = 2*sin(M_PI*tcur/periodT)*sin(M_PI*x)*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*y);


} //velocity_function

void find_nodes_on_exterior_boundary(std::shared_ptr<Jali::Mesh> mesh,
                                     std::vector<bool>& node_on_bnd) {

  int nfaces = mesh->num_entities(Jali::Entity_kind::FACE,
                                  Jali::Entity_type::ALL);
  int nnodes = mesh->num_entities(Jali::Entity_kind::NODE,
                                  Jali::Entity_type::ALL);

  node_on_bnd.assign(nnodes, false);

  for (int f = 0; f < nfaces; f++) {
    std::vector<int> fcells;
    mesh->face_get_cells(f, Jali::Entity_type::ALL, &fcells);

    if (fcells.size() == 1) {
      // if the face is on an exterior boundary, all its nodes are too
      std::vector<int> fnodes;
      mesh->face_get_nodes(f, &fnodes);
      for (auto const &n : fnodes)
        node_on_bnd[n] = true;
    }
  }
} //find_nodes_on_exterior_boundary
