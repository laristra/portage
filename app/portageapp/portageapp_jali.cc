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

// portage includes
#include "portage/support/portage.h"
#include "portage/support/mpi_collate.h"
#include "portage/driver/mmdriver.h"

// wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/support/Point.h" 
#include "wonton/support/wonton.h"

// For parsing and evaluating user defined expressions in apps

#include "user_field.h"

using Wonton::Jali_Mesh_Wrapper;
using Portage::argsort;
using Portage::reorder;

/*!  @file portageapp_jali.cc
  @brief A simple application that remaps fields between two meshes -
  the meshes can be internally generated rectangular meshes or
  externally read unstructured meshes

  This program is used to showcase our capabilities with various types
  of remap operations (e.g. interpolation order) on various types of
  meshes (2d or 3d; node-centered or zone-centered) of some simple
  linear or quadratic data.  For the cases of remapping linear data
  with a second-order interpolator, the L2 norm output at the end
  should be identically zero.
*/

//////////////////////////////////////////////////////////////////////


int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp " <<
      "--dim=2|3 --nsourcecells=N --ntargetcells=M --conformal=y|n \n" <<
      "--entity_kind=cell|node --field=\"your_math_expression\" --remap_order=1|2 \n" <<
      "--limiter=barth_jespersen --mesh_min=0. --mesh_max=1. \n" <<
      "--output_meshes=y|n --results_file=filename --convergence_study=NREF \n\n";

  std::cout << "--dim (default = 2): spatial dimension of mesh\n\n";

  std::cout << "--nsourcecells (NO DEFAULT): Internally generated rectangular "
            << "SOURCE mesh with num cells in each coord dir\n\n";
  std::cout << "        -----  OR  ----- \n";
  std::cout << "--source_file=srcfilename: file name of source mesh (Exodus II format only)\n\n";

  std::cout << "--ntargetcells (NO DEFAULT): Internally generated rectangular "
            << "TARGET mesh with num cells in each coord dir\n\n";
  std::cout << "        -----  OR  ----- \n";
  std::cout << "--target_file=trgfilename: file name of target mesh (Exodus II format only)\n\n";

  std::cout << "--conformal (default = y): 'y' means mesh boundaries match\n" <<
      "If 'n', target mesh is shifted in all directions by half-domain width\n";
  std::cout << " ONLY APPLICABLE TO INTERNALLY GENERATED TARGET MESH\n\n";

  std::cout << "--field (NO DEFAULT): A math expression enclosed in \"\" " <<
      " expressed in terms of \n" << "x, y and z following the syntax of " <<
      " the expression parser package ExprTk \n" <<
      "(http://www.partow.net/programming/exprtk/)\n";
  std::cout << "The syntax is generally predictable, e.g. " <<
      " \"24*x + y*y\" or \"x + log(abs(y)+1)\"\n";
  std::cout << "Advanced features include if-then-else and loop-statements\n\n";

  std::cout << "--entity_kind (default = cell): entities on which " <<
      "remapping field is to be done\n\n";

  std::cout << "--remap order (default = 1): " <<
      "order of accuracy of interpolation\n\n";

  std::cout << "--limiter (default = 0): " <<
      "slope limiter for a piecewise linear reconstrution\n\n";

  std::cout << "--mesh_min (default = 0.): " <<
      "coordinates (same in x, y, and z) of the lower corner of a mesh\n\n";

  std::cout << "--mesh_max (default = 1.): " <<
      "coordinates (same in x, y, and z) of the upper corner of a mesh\n\n";

  std::cout << "--convergence_study (default = 1): provide the number of times "
            << "you want to double source and target mesh sizes \n";
  std::cout << "  ONLY APPLICABLE IF BOTH MESHES ARE INTERNALLY GENERATED\n\n";

  std::cout << "--output_meshes (default = y)\n";
  std::cout << "  If 'y', the source and target meshes are output with the " <<
      "remapped field attached as input.exo and output.exo. \n\n";

  std::cout << "--results_file=results_filename (default = output.txt)\n";
  std::cout << "  If a filename is specified, the target field values are " <<
      "output to the file given by 'results_filename' in ascii format\n\n";
  return 0;
}
//////////////////////////////////////////////////////////////////////


// Forward declaration of function to run remap on two meshes and
// return the L1 and L2 error norm in the remapped field w.r.t. to an
// analytically imposed field. If no field was imposed, the errors are
// returned as 0

template<int dim> void run(std::shared_ptr<Jali::Mesh> sourceMesh,
                           std::shared_ptr<Jali::Mesh> targetMesh,
                           Portage::LimiterType limiter,
                           int interp_order,
                           std::string field_expression,
                           std::string field_filename,
                           bool mesh_output, int rank, int numpe,
                           Jali::Entity_kind entityKind,
                           double *L1_error, double *L2_error,
                           bool remap_back,
                           std::vector<double> &sourceData,
                           std::vector<double> &targetData);

int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  if (argc == 1) {
    print_usage();
    return -1;
  }

  struct timeval begin, end, diff;

  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int nsourcecells = 0, ntargetcells = 0;  // No default
  int dim = 2;
  bool conformal = true;
  std::string field_expression;
  std::string srcfile, trgfile;  // No default

  int interp_order = 1;
  bool mesh_output = true;
  int n_converge = 1;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  Portage::LimiterType limiter = Portage::LimiterType::NOLIMITER;
  double srclo = 0.0, srchi = 1.0;  // bounds of generated mesh in each dir
  bool remap_back = false;           // enable for cyclic  remap, i.e. remap
                                    // back to the original mesh

  std::string field_filename;  // No default

  // Parse the input

  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of("=");
    std::string keyword = arg.substr(keyword_beg, keyword_end-keyword_beg);
    std::string valueword = arg.substr(keyword_end+1, len-(keyword_end+1));

    if (keyword == "entity_kind") {
      if (valueword == "cell" || valueword == "CELL")
        entityKind = Jali::Entity_kind::CELL;
      else if (valueword == "node" || valueword == "NODE")
        entityKind = Jali::Entity_kind::NODE;
      else {
        std::cerr <<
            "Only node and cell based remapping supported at this time" <<
            std::endl;
        exit(-1);
      }
    } else if (keyword == "dim") {
      dim = stoi(valueword);
      assert(dim == 2 || dim == 3);
    } else if (keyword == "nsourcecells")
      nsourcecells = stoi(valueword);
    else if (keyword == "ntargetcells")
      if( stoi(valueword) == -1 )
        ntargetcells = nsourcecells + 1;
      else
        ntargetcells = stoi(valueword);
    else if (keyword == "source_file")
      srcfile = valueword;
    else if (keyword == "target_file")
      trgfile = valueword;
    else if (keyword == "conformal")
      conformal = (valueword == "y");
    else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
    } else if (keyword == "field") {
      field_expression = valueword;
    } else if (keyword == "limiter") {
      if (valueword == "barth_jespersen" || valueword == "BARTH_JESPERSEN")
        limiter = Portage::LimiterType::BARTH_JESPERSEN;
    } else if (keyword == "mesh_min") {
      srclo = stof(valueword);
    } else if (keyword == "mesh_max") {
      srchi = stof(valueword);
    } else if (keyword == "cyclic") {
      remap_back = (valueword == "y");
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
    } else if (keyword == "results_file") {
      field_filename = valueword;
    } else if (keyword == "convergence_study") {
      n_converge = stoi(valueword);
      if (n_converge <= 0) {
        std::cerr << "Number of meshes for convergence study should be greater than 0" << std::endl;
        throw std::exception();
      }
    } else if (keyword == "help") {
      print_usage();
    } else
      std::cerr << "Unrecognized option " << keyword << std::endl;
  }


  // Some input error checking

  if (nsourcecells > 0 && srcfile.length() > 0) {
    std::cout << "Cannot request internally generated source mesh "
              << "(--nsourcecells) and external file read (--source_file)\n\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if (!nsourcecells && srcfile.length() == 0) {
    std::cout << "Must specify one of the two options --nsourcecells "
              << "or --source_file\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (ntargetcells > 0 && trgfile.length() > 0) {
    std::cout << "Cannot request internally generated target mesh "
              << "(--ntargetcells) and external file read (--target_file)\n\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if (!ntargetcells && trgfile.length() == 0) {
    std::cout << "Must specify one of the two options --ntargetcells "
              << "or --target_file\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if ((srcfile.length() || trgfile.length()) && n_converge > 1) {
    std::cout <<
        "Convergence study possible only for internally generated meshes\n";
    std::cout << "Will do single remap and exit\n";
    n_converge = 1;
  }
  if (nsourcecells > 0 && field_expression.length() == 0) {
    std::cout << "No field imposed on internally generated source mesh\n";
    std::cout << "Nothing to remap. Exiting...";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }


  gettimeofday(&begin, 0);


  // The mesh factory and mesh setup
  std::shared_ptr<Jali::Mesh> source_mesh, target_mesh, cyclic_mesh;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});

  double trglo = srclo, trghi = srchi;  // bounds of generated mesh in each dir
  if (!conformal) {
    double dx = (trghi-trglo)/static_cast<double>(ntargetcells);
    trghi += 1.5*dx;
  }

  // If a mesh is being read from file, do it outside convergence loop
  if (srcfile.length() > 0) {
    mf.partitioner(Jali::Partitioner_type::METIS);
    source_mesh = mf(srcfile);
    cyclic_mesh = mf(srcfile);
  }
  if (trgfile.length() > 0) {
    mf.partitioner(Jali::Partitioner_type::METIS);
    target_mesh = mf(trgfile);
  }

  // partitioner for internally generated meshes is a BLOCK partitioner
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  std::vector<double> l1_err(n_converge, 0.0), l2_err(n_converge, 0.0);
  for (int i = 0; i < n_converge; i++) {

    // If a file is being internally generated, do it inside convergence loop
    if (nsourcecells)
      if (dim == 2) {
        source_mesh = mf(srclo, srclo, srchi, srchi, nsourcecells,
                         nsourcecells);
        cyclic_mesh = mf(srclo, srclo, srchi, srchi, nsourcecells,
                         nsourcecells);
      }
      else if (dim == 3) {
        source_mesh = mf(srclo, srclo, srclo, srchi, srchi, srchi,
                         nsourcecells, nsourcecells, nsourcecells);
        cyclic_mesh = mf(srclo, srclo, srclo, srchi, srchi, srchi,
                         nsourcecells, nsourcecells, nsourcecells);
      }
    if (ntargetcells)
      if (dim == 2)
        target_mesh = mf(trglo, trglo, trghi, trghi,
                         ntargetcells, ntargetcells);
      else if (dim == 3)
        target_mesh = mf(trglo, trglo, trglo, trghi, trghi, trghi,
                         ntargetcells, ntargetcells, ntargetcells);

    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    if (rank == 0) {
      if (n_converge == 1)
        std::cout << "Mesh Initialization Time: " << seconds << std::endl;
      else
        std::cout << "Mesh Initialization Time (Iteration i): " <<
            seconds << std::endl;
    }
    gettimeofday(&begin, 0);

    // Make sure we have the right dimension and that source and
    // target mesh dimensions match (important when one or both of the
    // meshes are read in)
    assert(source_mesh->space_dimension() == target_mesh->space_dimension());
    dim = source_mesh->space_dimension();

    // Store fields on initial and final meshes
    std::vector<double> sourceData, targetData, initialData;

    // Now run the remap on the meshes and get back the L2 error
    switch (dim) {
      case 2:
        run<2>(source_mesh, target_mesh, limiter, interp_order,
               field_expression, field_filename, mesh_output,
               rank, numpe, entityKind, &(l1_err[i]), &(l2_err[i]),
               false,sourceData,targetData);
        if(remap_back) {
          initialData = sourceData;
          run<2>(target_mesh, source_mesh, limiter, interp_order,
                 field_expression, field_filename, mesh_output,
                 rank, numpe, entityKind, &(l1_err[i]), &(l2_err[i]),
                 remap_back,targetData,sourceData);
        }
        break;
      case 3:
        run<3>(source_mesh, target_mesh, limiter, interp_order,
               field_expression, field_filename, mesh_output,
               rank, numpe, entityKind, &(l1_err[i]), &(l2_err[i]),
               false,sourceData,targetData);
        if(remap_back) {
          initialData = sourceData;
          run<3>(target_mesh, source_mesh, limiter, interp_order,
                 field_expression, field_filename, mesh_output,
                 rank, numpe, entityKind, &(l1_err[i]), &(l2_err[i]),
                 remap_back,targetData,sourceData);
        }

        break;
      default:
        std::cerr << "Dimension not 2 or 3" << std::endl;
        return 2;
    }

    if (remap_back) {
      Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
      double cyclic_error = 0.0, cyclic_error_norm = 0.0,
             max_initial = -1.0e50, max_final = -1.0e50,
             min_vol = 1.0e50, max_vol = -1.0e50;
      for (int j=0; j < initialData.size(); j++) {
        cyclic_error += fabs( initialData[j] - sourceData[j] ) *
                        sourceMeshWrapper.cell_volume(j);
        cyclic_error_norm += fabs( initialData[j] ) *
                        sourceMeshWrapper.cell_volume(j);
        max_initial = fmax( max_initial, initialData[j] );
        max_final = fmax( max_final, sourceData[j] );
        min_vol = fmin( min_vol, fabs( sourceMeshWrapper.cell_volume(j) - 1./double(nsourcecells*nsourcecells) ) );
        max_vol = fmax( max_vol, fabs( sourceMeshWrapper.cell_volume(j) - 1./double(nsourcecells*nsourcecells) ) );
      }
      std::cout << std::endl << std::endl << "CYCLIC REMAP RELATIVE ERROR = "
                << cyclic_error/cyclic_error_norm << " , NORM =  " <<
                cyclic_error_norm << " , MAX BOUND VIOLATION (RELATIVE) = " <<
                (max_final - max_initial)/max_initial << std::endl << std::endl;
      min_vol = min_vol * double(nsourcecells*nsourcecells);
      max_vol = max_vol * double(nsourcecells*nsourcecells);
      std::cout << "min_vol, max_vol = " << min_vol << " " << max_vol << std::endl;
    }
    else {
      std::cout << "L1 norm of error for iteration " << i << " is " <<
        l1_err[i] << std::endl;

      std::cout << "L2 norm of error for iteration " << i << " is " <<
        l2_err[i] << std::endl;
    }

    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    if (rank == 0) {
      if (n_converge == 1)
        std::cout << "Remap Time: " << seconds << std::endl;
      else
        std::cout << "Remap Time (Iteration i): " << seconds << std::endl;
    }
    gettimeofday(&begin, 0);


    // if convergence study, double the mesh resolution
    nsourcecells *= 2;
    ntargetcells *= 2;
  }

  for (int i = 1; i < n_converge; i++) {
    std::cout << "Error ratio L1(" << i-1 << ")/L1(" << i << ") is " <<
        l1_err[i-1]/l1_err[i] << std::endl;
    std::cout << "Error ratio L2(" << i-1 << ")/L2(" << i << ") is " <<
        l2_err[i-1]/l2_err[i] << std::endl;
  }

  MPI_Finalize();
}


// Run a remap between two meshes and return the L1 and L2 error norms
// with respect to the specified field. If a field was not specified and
//

template<int dim> void run(std::shared_ptr<Jali::Mesh> sourceMesh,
                           std::shared_ptr<Jali::Mesh> targetMesh,
                           Portage::LimiterType limiter,
                           int interp_order, std::string field_expression,
                           std::string field_filename, bool mesh_output,
                           int rank, int numpe, Jali::Entity_kind entityKind,
                           double *L1_error, double *L2_error, bool remap_back,
                           std::vector<double> &sourceData,
                           std::vector<double> &targetData) {

  // Wrappers for interfacing with the underlying mesh data structures.
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  const int nsrccells = sourceMeshWrapper.num_owned_cells() +
      sourceMeshWrapper.num_ghost_cells();
  const int ntarcells = targetMeshWrapper.num_owned_cells();

  const int nsrcnodes = sourceMeshWrapper.num_owned_nodes() +
      sourceMeshWrapper.num_ghost_nodes();
  const int ntarnodes = targetMeshWrapper.num_owned_nodes();

  if (rank == 0) {
    std::cout << "starting portageapp_jali...\n";
    std::cout << "All fields on" << entityKind << "including "
              << "the command line field (if specified) will be remapped\n";
    if (field_expression.length())
      std::cout << " Specified field is " << field_expression << "\n";
    std::cout << "   Interpolation order is " << interp_order << "\n";
    if (interp_order == 2)
      std::cout << "   Limiter type is " << limiter << "\n";
  }

  std::cout << "\nSource mesh on rank "<< rank <<" has " << nsrccells << " cells\n";
  std::cout << "Source mesh on rank "<< rank <<" has " << sourceMeshWrapper.num_owned_cells() << " owned cells\n";
  std::cout << "Target mesh on rank "<< rank <<" has " << ntarcells << " cells\n";
/*  
  for (int i = 0; i<nsrccells; i++){
    std::vector<Wonton::Point<2>> coords;
  	sourceMeshWrapper.cell_get_coordinates(i,&coords);
  	std::cout << "Src rank " << rank << " cell " << i ;
  	std::cout << " has gid " << sourceMeshWrapper.get_global_id(i, Wonton::Entity_kind::CELL);
  	std::cout << " : ";
  	std::cout << "(" << coords[0][0] << ", " << coords[0][1] << ") ";
  	std::cout << "(" << coords[1][0] << ", " << coords[1][1] << ") ";
  	std::cout << "(" << coords[2][0] << ", " << coords[2][1] << ") ";
  	std::cout << "(" << coords[3][0] << ", " << coords[3][1] << ")"<<std::endl;
  }
  
  for (int i = 0; i<ntarcells; i++){
    std::vector<Wonton::Point<2>> coords;
  	targetMeshWrapper.cell_get_coordinates(i,&coords);
  	std::cout << "tgt rank " << rank << " cell " << i ;
  	std::cout << " : ";
  	std::cout << "(" << coords[0][0] << ", " << coords[0][1] << ") ";
  	std::cout << "(" << coords[1][0] << ", " << coords[1][1] << ") ";
  	std::cout << "(" << coords[2][0] << ", " << coords[2][1] << ") ";
  	std::cout << "(" << coords[3][0] << ", " << coords[3][1] << ")"<<std::endl;
  }
*/  
  // Native jali state managers for source and target
  std::shared_ptr<Jali::State> sourceState(Jali::State::create(sourceMesh));
  std::shared_ptr<Jali::State> targetState(Jali::State::create(targetMesh));

  user_field_t source_field;
  if (!source_field.initialize(dim, field_expression))
    MPI_Abort(MPI_COMM_WORLD, -1);

  std::vector<std::string> remap_fields;

  // Cell-centered remaps
  if (entityKind == Jali::Entity_kind::CELL) {
    sourceData.resize(nsrccells);

    if (!remap_back) {
      for (unsigned int c = 0; c < nsrccells; ++c)
        sourceData[c] = source_field(sourceMesh->cell_centroid(c));
    }

    sourceState->add("celldata", sourceMesh, Jali::Entity_kind::CELL,
                     Jali::Entity_type::ALL, &(sourceData[0]));

    targetState->add<double, Jali::Mesh, Jali::UniStateVector>("celldata",
                                                            targetMesh,
                                                Jali::Entity_kind::CELL,
                                                Jali::Entity_type::ALL, 0.0);

    // Register the variable name and interpolation order with the driver
    remap_fields.push_back("celldata");

  } else {  // node-centered
    sourceData.resize(nsrcnodes);

    /*!
      @todo make node_get_coordinates be consistent in data type with
      cell_centroid?
    */

    Portage::Point<dim> node;
    for (int i = 0; i < nsrcnodes; ++i) {
      sourceMeshWrapper.node_get_coordinates(i, &node);
      sourceData[i] = source_field(node);
    }

    sourceState->add("nodedata", sourceMesh, Jali::Entity_kind::NODE,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    targetState->add<double, Jali::Mesh, Jali::UniStateVector>("nodedata",
                                                            targetMesh,
                                                Jali::Entity_kind::NODE,
                                                Jali::Entity_type::ALL, 0.0);

    // Register the variable name and remap order with the driver
    remap_fields.push_back("nodedata");

  }

  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);

  // Portage wrappers for source and target fields
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  if (dim == 2) {
    if (interp_order == 1) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR2D,
        Portage::Interpolate_1stOrder,
        2,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    } else if (interp_order == 2) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR2D,
        Portage::Interpolate_2ndOrder,
        2,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields, limiter);
      d.run(numpe > 1);
    }
  } else {  // 3D
    if (interp_order == 1) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_1stOrder,
        3,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    } else {  // 2nd order & 3D
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_2ndOrder,
        3,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields, limiter);
      d.run(numpe > 1);
    }
  }

  // Dump some timing information
  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);

  // Output results for small test cases
  double error, toterr = 0.0;
  double const * cellvecout;
  double const * cellvecin;
  double const * nodevecout;
  double const * nodevecin;
  double minin =  1.0e50, minout =  1.0e50;
  double maxin = -1.0e50, maxout = -1.0e50;
  double err_l1 = 0.;
  double err_norm = 0.;
  double target_mass = 0.;
  double source_mass = 0.;
  double totvolume = 0.;
  if (entityKind == Jali::Entity_kind::CELL)  {  // CELL error computation
    targetStateWrapper.mesh_get_data<double>(Portage::CELL, "celldata",
                                        &cellvecout);
    sourceStateWrapper.mesh_get_data<double>(Portage::CELL, "celldata",
                                        &cellvecin);
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
      targetMeshWrapper.cell_centroid(c, &ccen);
      error = source_field(ccen) - cellvecout[c];
      targetData[c] = cellvecout[c];
      minout = fmin( minout, cellvecout[c] );
      maxout = fmax( maxout, cellvecout[c] );
      double cellvol = targetMeshWrapper.cell_volume(c);
      err_l1 += fabs(error)*cellvol;
      err_norm  += fabs( cellvecout[c] ) * cellvol;
      target_mass += cellvecout[c] * cellvol;

      if (!targetMeshWrapper.on_exterior_boundary(Portage::Entity_kind::CELL, c)) {
        double cellvol = targetMeshWrapper.cell_volume(c);
        totvolume += cellvol;
        *L1_error += fabs(error)*cellvol;
        *L2_error += error*error*cellvol;
      }
      if (ntarcells < 10) {
    	std::printf("Rank %d\n", rank);
        std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                    ccen[0], ccen[1]);
        std::printf("  Value = % 10.6lf  L2 Err = % lf\n",
                    cellvecout[c], error);
      }
    }
  } else {  // NODE error computation
    targetStateWrapper.mesh_get_data<double>(Portage::NODE, "nodedata",
                                        &nodevecout);
    sourceStateWrapper.mesh_get_data<double>(Portage::NODE, "nodedata",
                                        &nodevecin);
    if (numpe == 1 && ntarnodes < 10)
      std::cout << "nodedata vector on target mesh after remapping is:"
                << std::endl;

    for (int c = 0; c < nsrcnodes; ++c) {
      minin = fmin( minin, nodevecin[c] );
      maxin = fmax( maxin, nodevecin[c] );
      double nodevol2 = sourceMeshWrapper.dual_cell_volume(c);
      source_mass += nodevecin[c] * nodevol2;
    }

    Portage::Point<dim> nodexy;
    for (int i = 0; i < ntarnodes; ++i) {
      targetMeshWrapper.node_get_coordinates(i, &nodexy);
      error = source_field(nodexy) - nodevecout[i];
      double dualcellvol = targetMeshWrapper.dual_cell_volume(i);
      minout = fmin( minout, nodevecout[i] );
      maxout = fmax( maxout, nodevecout[i] );
      err_l1 += fabs(error)*dualcellvol;
      //std::printf("Source/Target volumes = %.5e, %.5e \n",
      //      sourceMeshWrapper.dual_cell_volume(i),
      //      targetMeshWrapper.dual_cell_volume(i));
      err_norm  += fabs( nodevecout[i] ) * dualcellvol;
      target_mass += nodevecout[i] * dualcellvol;

      if (!targetMeshWrapper.on_exterior_boundary(Portage::Entity_kind::NODE, i)) {
        totvolume += dualcellvol;
        *L1_error += fabs(error)*dualcellvol;
        *L2_error += error*error*dualcellvol;
      }
      if (ntarnodes < 10) {
        std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf) ", i,
                    nodexy[0], nodexy[1]);
        std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], error);
      }
    }
  }

  if (numpe > 1) {
    std::cout << std::flush << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    double globalerr;
    MPI_Reduce(L1_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    *L1_error = globalerr;

    MPI_Reduce(L2_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    *L2_error = globalerr;
  }
  err_norm = err_l1 / err_norm;
  *L2_error = sqrt(*L2_error);
  if (!remap_back) {
    std::printf("\n\nL1 NORM OF ERROR (excluding boundary) = %lf\n", *L1_error);
    std::printf("L2 NORM OF ERROR (excluding boundary) = %lf\n\n", *L2_error);
    std::printf("===================================================\n");
    std::printf("ON RANK %d\n", rank);
    std::printf("Relative L1 error = %.5e \n", err_norm);
    std::printf("Source min/max    = %.15e %.15e \n", minin, maxin);
    std::printf("Target min/max    = %.15e %.15e \n", minout, maxout);
    std::printf("Source total mass = %.15e \n", source_mass);
    std::printf("Target total mass = %.15e \n", target_mass);
    std::printf("===================================================\n");
  }

  // Write out the meshes if requested
  if (mesh_output) {

    // The current version of MSTK (2.27rc2) has a bug in writing out
    // exodus files with node variables in parallel and so we will avoid
    // the exodus export in this situation. The 'if' statement can be
    // removed once we upgrade to the next version of MSTK

    if (numpe == 1) {
      if (rank == 0)
        std::cout << "Dumping data to Exodus files..." << std::endl;
      if (field_expression.length() > 0) {
        sourceState->export_to_mesh();
        sourceMesh->write_to_exodus_file("input.exo");
      }
      targetState->export_to_mesh();
      targetMesh->write_to_exodus_file("output.exo");
      if (rank == 0)
        std::cout << "...done." << std::endl;
    }
  }  // if (dump_meshes)

  // construct the field file name and open the file

  if (field_filename.length()) {
    std::vector<int> lgid;
    std::vector<double> lvalues;
    std::string entstr;
    if (entityKind == Jali::Entity_kind::CELL) {
      entstr = "cell";
      lgid.resize(ntarcells);
      lvalues.resize(ntarcells);
      for (int i=0; i < ntarcells; i++) {
        lgid[i] = targetMesh->GID(i, Jali::Entity_kind::CELL);
        lvalues[i] = cellvecout[i];
      }
    } else {
      entstr = "node";
      lgid.resize(ntarnodes);
      lvalues.resize(ntarnodes);
      for (int i=0; i < ntarnodes; i++) {
        lgid[i] = targetMesh->GID(i, Jali::Entity_kind::NODE);
        lvalues[i] = nodevecout[i];
      }
    }

    // sort the field values by global ID
    std::vector<int> idx;
    argsort(lgid, idx);   // find sorting indices based on global IDS
    reorder(lgid, idx);   // sort the global ids
    reorder(lvalues, idx);  // sort the values

    if (numpe > 1) {
      int maxwidth = static_cast<long long>(std::ceil(std::log10(numpe)));
      char rankstr[10];
      std::snprintf(rankstr, sizeof(rankstr), "%0*d", maxwidth, rank);
      field_filename = field_filename + "." + std::string(rankstr);
    }
    std::ofstream fout(field_filename);
    fout << std::scientific;
    fout.precision(17);

    // write out the values
    for (int i=0; i < lgid.size(); i++)
      fout << lgid[i] << " " << lvalues[i] << std::endl;
  }
}
