/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <sys/time.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <utility>
#include <cmath>
#include <set>
#include <numeric>
#include <iomanip>
#include <limits>

#include <mpi.h>

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

// Jali mesh infrastructure library
// See https://github.com/lanl/jali
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

// Wonton
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// Portage
#include "portage/driver/coredriver.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/search/search_swept_face.h"
#include "portage/search/search_kdtree.h"
#include "portage/support/portage.h"
#include "portage/support/timer.h"
#include "portage/support/mpi_collate.h"

#if ENABLE_TIMINGS
#include "portage/support/timer.h"
#endif

//#define DEBUG_PRINT

// For parsing and evaluating user defined expressions in apps
#include "user_field.h"

using Wonton::Jali_Mesh_Wrapper;

# define PI 3.14159265358979323846  /* pi */


/*!
  @file sweptfaceapp_jali.cc

  @brief A simple application that remaps using a swept face algorithm
  for single material field between two meshes - the meshes can be
  internally generated rectangular meshes or externally read unstructured
  meshes.
*/

//////////////////////////////////////////////////////////////////////
// Forward declarations

int print_usage();

// Forward declaration of function to run remap on two meshes and
// return the L1 and L2 error norm in the remapped field w.r.t. to an
// analytically imposed field. If no field was imposed, the errors are
// returned as 0

template<int dim> void remap(std::shared_ptr<Jali::Mesh> sourceMesh,
                             std::shared_ptr<Jali::Mesh> targetMesh,
                             std::string field_expression,
                             Jali::Entity_kind entityKind,
                             int interp_order,
                             Portage::Limiter_type limiter,
                             Portage::Boundary_Limiter_type bnd_limiter,
                             bool mesh_output,
                             std::string field_filename,
                             int rank, int numpe, int iteration,
                             double& L1_error, double& L2_error,
                             bool sweptface,
                             std::shared_ptr<Profiler> profiler);

// Forward declaration of function to find nodes on exterior boundary
void find_nodes_on_exterior_boundary(std::shared_ptr<Jali::Mesh> mesh,
                                     std::vector<bool>& node_on_bnd);

void move_target_mesh_nodes(std::shared_ptr<Jali::Mesh> mesh,
                            double& tcur, double& deltaT, double& periodT, int& scale);

void single_vortex_velocity_function(double* coords, double& tcur,
                                     double& periodT, double* disp);

//////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc == 1) {
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int dim = 2; //dim 3 is not currently supported
  int ncells = 10;

  std::string field_expression;
  std::string field_filename="";  // No default;

  int interp_order = 1;
  bool sweptface = true;
  int ntimesteps = 4;
  int scale = 20;
  // since write_to_gmv segfaults in parallel, default to false and force the
  // user to output in serial
  bool mesh_output = false;

  //Node-based remap not supported for swept-face algorithm
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  Portage::Limiter_type limiter = Portage::Limiter_type::NOLIMITER;
  Portage::Boundary_Limiter_type bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;

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
    } else if (keyword == "sweptface") {
      sweptface = (valueword == "y");
    } else if (keyword == "field_expression") {
      field_expression = valueword;
    } else if (keyword == "limiter") {
      if (valueword == "barth_jespersen" || valueword == "BARTH_JESPERSEN")
        limiter = Portage::Limiter_type::BARTH_JESPERSEN;
    } else if (keyword == "bnd_limiter") {
      if (valueword == "zero_gradient" || valueword == "ZERO_GRADIENT")
        bnd_limiter = Portage::Boundary_Limiter_type::BND_ZERO_GRADIENT;
      else if (valueword == "barth_jespersen" || valueword == "BARTH_JESPERSEN")
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
      MPI_Abort(MPI_COMM_WORLD, -1);
    } else if (keyword == "field_filename") {
      field_filename=valueword;
    } else {
      std::cerr << "Unrecognized option " << keyword << " !\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
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

  // The mesh factory and mesh setup
  int nsourcecells = ncells, ntargetcells = ncells;
  double srclo = 0.0, srchi = 1.0;  // bounds of generated mesh in each dir
  double trglo = srclo, trghi = srchi;  // bounds of generated mesh in each dir

  std::shared_ptr<Jali::Mesh> source_mesh, target_mesh;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});

  // partitioner for internally generated meshes is a BLOCK partitioner
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  // Create the source and target meshes
  if (dim == 2) {
    source_mesh = mf(srclo, srclo, srchi, srchi, nsourcecells, nsourcecells);
    target_mesh = mf(trglo, trglo, trghi, trghi, ntargetcells, ntargetcells);
  }
  else if (dim == 3) {
    source_mesh = mf(srclo, srclo, srclo, srchi, srchi, srchi,
                     nsourcecells, nsourcecells, nsourcecells);
    target_mesh = mf(trglo, trglo, trglo, trghi, trghi, trghi,
                     ntargetcells, ntargetcells, ntargetcells);
  }


  // Make sure we have the right dimension and that source and
  // target mesh dimensions match (important when one or both of the
  // meshes are read in)
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  dim = source_mesh->space_dimension();

#if ENABLE_TIMINGS
  profiler->time.mesh_init = timer::elapsed(tic);

  if (rank == 0) {
        std::cout << "Mesh Initialization Time: " << profiler->time.mesh_init << std::endl;
    }
#endif

  double periodT = 2.0;
  double deltaT = periodT/ntimesteps;

  std::vector<double> l1_err(ntimesteps, 0.0), l2_err(ntimesteps, 0.0);

  for (int i = 1; i < ntimesteps; i++) {
    // Move nodes of the target mesh
    std::cout<<"Start moving mesh for timestep "<<i<<" scaled by "<<scale<<std::endl;

    double tcur = i*deltaT;
    move_target_mesh_nodes(target_mesh, tcur, deltaT, periodT, scale);

    std::cout<<"Completed moving mesh for timestep "<<i<<std::endl;
    // Now run the remap on the meshes and get back the L1-L2 errors
    switch (dim) {
      case 2:
        remap<2> (source_mesh, target_mesh, field_expression, entityKind,
                  interp_order, limiter, bnd_limiter, mesh_output,
                  field_filename, rank, numpe, i, l1_err[i], l2_err[i],
                  sweptface, profiler);
        break;
      case 3:
        //remap<3> (source_mesh, target_mesh, field_expression, entityKind,
        //        interp_order, limiter, bnd_limiter, mesh_output,
        //        field_filename, rank, numpe, l1_err[i], l2_err[i],
        //        sweptface, profiler);

        break;
      default:
        std::cerr << "Dimension not 2 or 3" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    std::cout << "L1 norm of error for timestep " << i << " is " <<
              l1_err[i] << std::endl;
    std::cout << "L2 norm of error for timestep " << i << " is " <<
              l2_err[i] << std::endl;

  }

  MPI_Finalize();

#if ENABLE_TIMINGS
  profiler->time.total = timer::elapsed(start);

  // dump timing data
  if (rank == 0) {
    profiler->dump();
  }
#endif
}
//////////////////////////////////////////////////////////////////////////////
// Run a remap between two meshes and return the L1 and L2 error norms
// with respect to the specified field.

template<int dim> void remap(std::shared_ptr<Jali::Mesh> sourceMesh,
                             std::shared_ptr<Jali::Mesh> targetMesh,
                             std::string field_expression,
                             Jali::Entity_kind entityKind,
                             int interp_order,
                             Portage::Limiter_type limiter,
                             Portage::Boundary_Limiter_type bnd_limiter,
                             bool mesh_output,
                             std::string field_filename,
                             int rank, int numpe, int iteration,
                             double& L1_error, double& L2_error,
                             bool sweptface,
                             std::shared_ptr<Profiler> profiler) {
  if (rank == 0)
    std::cout << "starting sweptfacedemo_jali_singlemat...\n";

  if (entityKind != Jali::Entity_kind::CELL) {
    std::cerr << "Sweptfacedemo not supported for node based fields!\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // Wrappers for interfacing with the underlying mesh data structures.
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  const int nsrccells = sourceMeshWrapper.num_owned_cells() +
                        sourceMeshWrapper.num_ghost_cells();
  const int ntarcells = targetMeshWrapper.num_owned_cells() +
                        targetMeshWrapper.num_ghost_cells();

  // Native jali state managers for source and target
  std::shared_ptr<Jali::State> sourceState(Jali::State::create(sourceMesh));
  std::shared_ptr<Jali::State> targetState(Jali::State::create(targetMesh));

  // Portage wrappers for source and target fields
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  // Output some information for the user
  if (rank == 0) {
    std::cout << "Source mesh has " << nsrccells << " cells\n";
    std::cout << "Target mesh has " << ntarcells << " cells\n";
    std::cout << " Specified single material field is ";
    std::cout << "    " << field_expression << ", ";
    std::cout << "\n";
    std::cout << "   Interpolation order is " << interp_order << "\n";
    if (interp_order == 2) {
      std::cout << "   Limiter type is " << limiter << "\n";
      std::cout << "   Boundary limiter type is " << bnd_limiter << "\n";
    }
  }

#if ENABLE_TIMINGS
  auto tic = timer::now();
#endif

  // Compute field data on source and add to state manager.
  user_field_t source_field;
  if (!source_field.initialize(dim, field_expression))
    MPI_Abort(MPI_COMM_WORLD, -1);

  std::vector<double> sourceData, targetData;
  sourceData.resize(nsrccells);
  for (unsigned int c = 0; c < nsrccells; ++c)
    sourceData[c] = source_field(sourceMesh->cell_centroid(c));

  sourceState->add("density", sourceMesh, Jali::Entity_kind::CELL,
                   Jali::Entity_type::ALL, &(sourceData[0]));

  targetState->add<double, Jali::Mesh, Jali::UniStateVector>("density",
                                                             targetMesh,
                                                             Jali::Entity_kind::CELL,
                                                             Jali::Entity_type::ALL, 0.0);

#if ENABLE_TIMINGS
  profiler->time.remap = timer::elapsed(tic, true);
#endif

  // Executor
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  Wonton::Executor_type *executor = (numpe > 1) ? &mpiexecutor : nullptr;

  if (rank == 0) {
    std::cout << "***Registered fieldnames:\n";
    for (auto & field_name: sourceStateWrapper.names())
      std::cout << " registered fieldname: " << field_name << std::endl;
  }

  std::vector<std::string> fieldnames;
  fieldnames.push_back("density");

  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);

  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
    cdriver(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper, executor);

  if (sweptface) {
    auto candidates = cdriver.search<Portage::SearchSweptFace>();
    auto weights = cdriver.intersect_meshes<Portage::IntersectSweptFace2D>(candidates);
    bool has_mismatch = cdriver.check_mesh_mismatch(weights);

    double dblmin = -std::numeric_limits<double>::max();
    double dblmax =  std::numeric_limits<double>::max();

    if (interp_order == 1) {
      cdriver.interpolate_mesh_var<double,
        Portage::Interpolate_1stOrder>("density",
                                       "density",
                                       weights,
                                       0.0, dblmax);
    } else if (interp_order == 2) {
      cdriver.interpolate_mesh_var<double,
        Portage::Interpolate_2ndOrder>("density",
                                       "density",
                                       weights,
                                       0.0, dblmax);
    }
  }
  else {
    auto candidates = cdriver.search<Portage::SearchKDTree>();
    auto weights = cdriver.intersect_meshes<Portage::IntersectR2D>(candidates);
    bool has_mismatch = cdriver.check_mesh_mismatch(weights);

    double dblmin = -std::numeric_limits<double>::max();
    double dblmax =  std::numeric_limits<double>::max();

    if (interp_order == 1) {
      cdriver.interpolate_mesh_var<double,
        Portage::Interpolate_1stOrder>("density",
                                       "density",
                                       weights,
                                       0.0, dblmax);
    } else if (interp_order == 2) {
      cdriver.interpolate_mesh_var<double,
        Portage::Interpolate_2ndOrder>("density",
                                       "density",
                                       weights,
                                       0.0, dblmax);
    }
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
      std::cout << "Dumping data to Exodus files..." << std::endl;

    std::string filename = "source_" + std::to_string(rank) + std::to_string(iteration)+ ".exo";
    sourceState->export_to_mesh();
    sourceMesh->write_to_exodus_file(filename);

    filename = "target_" + std::to_string(rank) + std::to_string(iteration)+ ".exo";
    targetState->export_to_mesh();
    targetMesh->write_to_exodus_file(filename);

    if (rank == 0)
      std::cout << "...done." << std::endl;
  }


  /*if (mesh_output) {
    std::string filename = "source_" + std::to_string(rank) + ".gmv";
    Portage::write_to_gmv<dim>(sourceMeshWrapper, sourceStateWrapper,
             fieldnames,
             filename);

    filename = "target_" + std::to_string(rank) + ".gmv";
    Portage::write_to_gmv<dim>(targetMeshWrapper, targetStateWrapper,
             fieldnames,
             filename);
  }*/

  // Compute error
  L1_error = 0.0; L2_error = 0.0;

  double error, toterr = 0.0;
  double const * cellvecout;
  double const * cellvecin;
  double minin =  1.0e50, minout =  1.0e50;
  double maxin = -1.0e50, maxout = -1.0e50;
  double err_l1 = 0.;
  double err_norm = 0.;
  double target_mass = 0.;
  double source_mass = 0.;
  double totvolume = 0. ;

  sourceStateWrapper.mesh_get_data<double>(Portage::CELL, "density",
                                           &cellvecin);
  targetStateWrapper.mesh_get_data<double>(Portage::CELL, "density",
                                           &cellvecout);

  if (numpe == 1 && ntarcells < 10)
    std::cout << "celldata vector on target mesh after remapping is:"
              << std::endl;

  targetData.resize(ntarcells);

  // Compute total mass on the source mesh to check conservation
  for (int c = 0; c < nsrccells; ++c) {
    minin = fmin( minin, cellvecin[c] );
    maxin = fmax( maxin, cellvecin[c] );
    double cellvol2 = sourceMeshWrapper.cell_volume(c);
    source_mass += cellvecin[c] * cellvol2;
  }

  // Cell error computation
  Portage::Point<dim> ccen;
  for (int c = 0; c < ntarcells; ++c) {
    if (targetMeshWrapper.cell_get_type(c) == Portage::Entity_type::PARALLEL_OWNED) {

      targetMeshWrapper.cell_centroid(c, &ccen);
      double cellvol = targetMeshWrapper.cell_volume(c);
      targetData[c] = cellvecout[c];
      minout = fmin( minout, cellvecout[c] );
      maxout = fmax( maxout, cellvecout[c] );
      error = cellvecin[c] - cellvecout[c];
      L1_error += fabs(error)*cellvol;
      L2_error += error*error*cellvol;
      totvolume += cellvol;
      target_mass += cellvecout[c]*cellvol;

      if (ntarcells < 10) {
        std::printf("Rank %d\n", rank);
        std::printf("Cell=% 4d Centroid = (% 8.5lf,% 8.5lf)", c,
                    ccen[0], ccen[1]);
        std::printf("  Value = % 10.6lf  L2 Err = % lf\n",
                    cellvecout[c], error);
      }
    }
  }

  if (numpe > 1) {
    std::cout << std::flush << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    double globalerr;
    MPI_Reduce(&L1_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    L1_error = globalerr;

    MPI_Reduce(&L2_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    L2_error = globalerr;
  }
  L2_error = sqrt(L2_error);

  std::printf("\n\nL1 NORM OF ERROR (excluding boundary) = %lf\n", L1_error);
  std::printf("L2 NORM OF ERROR (excluding boundary) = %lf\n\n", L2_error);
  std::printf("===================================================\n");
  std::printf("ON RANK %d\n", rank);
  std::printf("Relative L1 error = %.5e \n", err_norm);
  std::printf("Source min/max    = %.15e %.15e \n", minin, maxin);
  std::printf("Target min/max    = %.15e %.15e \n", minout, maxout);
  std::printf("Source total mass = %.15e \n", source_mass);
  std::printf("Target total mass = %.15e \n", target_mass);
  std::printf("===================================================\n");

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
                            double& tcur, double& deltaT, double& periodT, int& scale)
{
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

  for (int i = 0; i <ntarnodes; i++)
  {
    std::array<double, 2> coords;
    mesh->node_get_coordinates(i, &coords);

    if (ntarnodes <= 20)
      std::cout<<"Target Node = "<<i<<" Original Coords = {"<<coords[0]
               <<", "<<coords[1]<<"}"<<std::endl;

    if (!node_on_bnd[i]) {

      std::array<double, 2> disp;
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
} //move_target_mesh_nodes

void single_vortex_velocity_function(double* coords, double& tcur,
                                     double& periodT, double* disp)
{
  double x = coords[0];
  double y = coords[1];
  disp[0] = -2*sin(PI*tcur/periodT)*sin(PI*x)*sin(PI*x)*sin(PI*y)*cos(PI*y);
  disp[1] = 2*sin(PI*tcur/periodT)*sin(PI*x)*cos(PI*x)*sin(PI*y)*sin(PI*y);


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

int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: sweptfacedemo_jali_singlemat " <<
            "--dim=2|3 --ncells=N --remap_order=1|2\n"<<
            "--sweptface=y|n --field=expression\n" <<
            "--limiter=barth_jespersen --bnd_limiter=zero_gradient \n"
            "--output_meshes=y|n --field_filename=string \n" <<
            "--only_threads=y|n --scaling=strong|weak \n\n";

  std::cout << "--dim (default = 2): spatial dimension of mesh\n\n";

  std::cout << "--ncells (DEFAULT 10): Internally generated rectangular "
            << "SOURCE and TARGET mesh with num cells in each coord dir\n\n";

  std::cout << "--remap_order (default = 1): " <<
            "order of accuracy of interpolation\n\n";

  std::cout << "--sweptface (default = y): " <<
            "use sweptface based remap, if no use intersection based remap\n\n";

  std::cout << "--field-expression (NO DEFAULT): A quoted math expressions \n" <<
            " expressed in terms of \n" << "x, y and z following the syntax of " <<
            " the expression parser package ExprTk \n" <<
            "(http://www.partow.net/programming/exprtk/)\n";
  std::cout << "The syntax is generally predictable, e.g. " <<
            " \"24*x + y*y\" or \"x + log(abs(y)+1)\"\n";

  std::cout << "--limiter (default = NOLIMITER): " <<
            "slope limiter for a piecewise linear reconstrution\n\n";

  std::cout << "--bnd_limiter (default = NOLIMITER): " <<
            "slope limiter on the boundary for a piecewise linear reconstruction\n\n";

  std::cout << "--output_meshes (default = y)\n";
  std::cout << "  If 'y', the source and target meshes are output with the " <<
            "remapped field attached as input.exo and output.exo, \n" <<
            "and interface reconstruction results are output as source_mm.gmv and target_mm.gmv \n\n";

  std::cout << "--field_filename\n\n";
  std::cout << "  If defined, the field output filename. Rank is appended if parallel\n\n";

#if ENABLE_TIMINGS
  std::cout << "--only_threads (default = n)\n";
  std::cout << " enable if you want to profile only threads scaling\n\n";

  std::cout << "--scaling (default = strong)\n";
  std::cout << " specify the scaling study type [strong|weak]\n\n";
#endif

  return 0;
} //print_usage