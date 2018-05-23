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


#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/mpi_collate.h"
#include "portage/driver/mmdriver.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/jali/jali_state_wrapper.h"
#include "read_material_data.h"

#ifdef XMOF2D
#endif

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/driver/write_to_gmv.h"
#endif

// For parsing and evaluating user defined expressions in apps

#include "user_field.h"

using Wonton::Jali_Mesh_Wrapper;
using Portage::argsort;
using Portage::reorder;

/*! 
  @file portageapp_jali_multimat.cc 

  @brief A simple application that remaps multi-material fields
  between two meshes - the meshes can be internally generated
  rectangular meshes or externally read unstructured meshes

  This program is used to showcase our capabilities with various types
  of remap operations (e.g. interpolation order) on various types of
  meshes (2d or 3d; node-centered or zone-centered) with multiple
  materials. The material data (material volume fractions, material
  centroids and field data) are specified in a file with a .mat
  extension

*/

//////////////////////////////////////////////////////////////////////


int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp " <<
      "--dim=2|3 --nsourcecells=N --ntargetcells=M --conformal=y|n \n" << 
      "--remap_order=1|2 \n" <<
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

  std::cout << "--mesh_min (default = 0.): " <<
      "coordinates (same in x, y, and z) of the lower corner of a mesh\n" <<
      "ONLY APPLICABLE FOR INTERNALLY GENERATED MESHES\n\n";

  std::cout << "--mesh_max (default = 1.): " <<
      "coordinates (same in x, y, and z) of the upper corner of a mesh\n" <<
      "ONLY APPLICABLE FOR INTERNALLY GENERATED MESHES\n\n";

  std::cout << "--conformal (default = y): 'y' means mesh boundaries match\n" <<
      "If 'n', target mesh is shifted in all directions by half-domain width\n";
  std::cout << " ONLY APPLICABLE TO INTERNALLY GENERATED TARGET MESH\n\n";

  std::cout << "--material_file=filename: file containing volume fractions\n" <<
      "and, optionally, material centroids for cells. The file lists the\n";
  std::cout << "total number of materials and a numeric value indicating\n" <<
      "if only volume fractions are given (1) or volume fractions and\n" <<
      "centroids are given (2). Subsequent lines list the cell ID,\n" <<
      "number of materials in the cell, the material IDs," <<
      " the volume fractions\n" << " and the centroids for the materials.\n\n";

  std::cout << "--material_fields: A comma separated list of quoted math expressions \n" <<
      " expressed in terms of \n" << "x, y and z following the syntax of " <<
      " the expression parser package ExprTk \n" <<
      "(http://www.partow.net/programming/exprtk/)\n";
  std::cout << "The syntax is generally predictable, e.g. " <<
      " \"24*x + y*y\" or \"x + log(abs(y)+1)\"\n";
  std::cout << "There must be as many expressions as there materials in the problem\n\n";

  std::cout << "--remap order (default = 1): " <<
      "order of accuracy of interpolation\n\n";

  std::cout << "--limiter (default = 0): " <<
      "slope limiter for a piecewise linear reconstrution\n\n";

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
                           std::string material_filename,
                           std::vector<std::string> material_field_expressions,
                           std::string field_filename,
                           bool mesh_output, int rank, int numpe,
                           Jali::Entity_kind entityKind,
                           double *L1_error, double *L2_error);

int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  if (argc == 1) print_usage();
  
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
  std::string material_filename;
  std::vector<std::string> material_field_expressions;
  std::string srcfile, trgfile;  // No default

  int interp_order = 1;
  bool mesh_output = true;
  int n_converge = 1;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  Portage::LimiterType limiter = Portage::LimiterType::NOLIMITER;
  double srclo = 0.0, srchi = 1.0;  // bounds of generated mesh in each dir

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
    } else if (keyword == "nsourcecells")
      nsourcecells = stoi(valueword);
    else if (keyword == "ntargetcells")
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
    } else if (keyword == "material_file") {
      material_filename = valueword;
    } else if (keyword == "material_fields") {
      // Expecting comma-separated quoted expressions
      std::string exprlist = valueword;
      std::size_t expr_beg = 0;
      std::size_t expr_end = 0;
      bool exprdone = false;
      std::string expr;
      while (!exprdone) {
        expr_end = exprlist.find_first_of(",", expr_beg);
        if (expr_end == std::string::npos) {
          exprdone = true;
          expr = exprlist.substr(expr_beg);
        } else {
          expr = exprlist.substr(expr_beg, expr_end-expr_beg);
          expr_beg = expr_end+1;
        }
        material_field_expressions.push_back(expr);
      }
    } else if (keyword == "limiter") {
      if (valueword == "barth_jespersen" || valueword == "BARTH_JESPERSEN")
        limiter = Portage::LimiterType::BARTH_JESPERSEN;
    } else if (keyword == "mesh_min") {
      srclo = stof(valueword);
    } else if (keyword == "mesh_max") {
      srchi = stof(valueword);
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
    } else if (keyword == "results_file") {
      field_output_filename = valueword;
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
  if (nsourcecells > 0 && material_field_expressions.size() == 0) {
    std::cout << "No field imposed on internally generated source mesh\n";
    std::cout << "Nothing to remap. Exiting...";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }


  gettimeofday(&begin, 0);


  // The mesh factory and mesh setup
  std::shared_ptr<Jali::Mesh> source_mesh, target_mesh;
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
    if (nsourcecells) {
      if (dim == 2)
        source_mesh = mf(srclo, srclo, srchi, srchi, nsourcecells,
                         nsourcecells);
      else if (dim == 3)
        source_mesh = mf(srclo, srclo, srclo, srchi, srchi, srchi,
                         nsourcecells, nsourcecells, nsourcecells);
    }
    if (ntargetcells) {
      if (dim == 2)
        target_mesh = mf(trglo, trglo, trghi, trghi,
                         ntargetcells, ntargetcells);
      else if (dim == 3)
        target_mesh = mf(trglo, trglo, trglo, trghi, trghi, trghi,
                         ntargetcells, ntargetcells, ntargetcells);
    }
    
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

    double error_L1 = 0.0, error_L2 = 0.0;

    // Now run the remap on the meshes and get back the L2 error
    switch (dim) {
      case 2:
        run<2>(source_mesh, target_mesh, limiter, interp_order,
               material_filename, material_field_expressions,
               field_output_filename, mesh_output,
               rank, numpe, entityKind, &(l1_err[i]), &(l2_err[i]));
        break;
      case 3:
        run<3>(source_mesh, target_mesh, limiter, interp_order,
               material_filename, material_field_expressions,
               field_output_filename, mesh_output,
               rank, numpe, entityKind, &(l1_err[i]), &(l2_err[i]));
        break;
      default:
        std::cerr << "Dimension not 2 or 3" << std::endl;
        return 2;
    }
    std::cout << "L1 norm of error for iteration " << i << " is " <<
        l2_err[i] << std::endl;
    std::cout << "L2 norm of error for iteration " << i << " is " <<
        l2_err[i] << std::endl;

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
// remap only volume fractions (and if specified, centroids)

template<int dim> void run(std::shared_ptr<Jali::Mesh> sourceMesh,
                           std::shared_ptr<Jali::Mesh> targetMesh,
                           Portage::LimiterType limiter,
                           int interp_order,
                           std::string material_filename,
                           std::vector<std::string> material_field_expressions,
                           std::string field_filename, bool mesh_output,
                           int rank, int numpe, Jali::Entity_kind entityKind,
                           double *L1_error, double *L2_error) {

  if (rank == 0)
    std::cout << "starting portageapp_jali_multimat...\n";

  // Wrappers for interfacing with the underlying mesh data structures.
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  const int nsrccells = sourceMeshWrapper.num_owned_cells() +
      sourceMeshWrapper.num_ghost_cells();
  const int ntarcells = targetMeshWrapper.num_owned_cells();

  const int nsrcnodes = sourceMeshWrapper.num_owned_nodes() +
      sourceMeshWrapper.num_ghost_nodes();
  const int ntarnodes = targetMeshWrapper.num_owned_nodes();

  
  // Native jali state managers for source and target
  std::shared_ptr<Jali::State> sourceState(Jali::State::create(sourceMesh));
  std::shared_ptr<Jali::State> targetState(Jali::State::create(targetMesh));

  // Portage wrappers for source and target fields
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  
  // Read volume fraction and centroid data from file
  
  int nmats;
  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;  // flattened 2D array
  std::vector<double> cell_mat_volfracs;  // flattened 2D array
  std::vector<Portage::Point<dim>> cell_mat_centroids;  // flattened 2D array
  
  read_material_data<Wonton::Jali_Mesh_Wrapper, dim>(sourceMeshWrapper,
                                                     material_filename,
                                                     cell_num_mats,
                                                     cell_mat_ids,
                                                     cell_mat_volfracs,
                                                     cell_mat_centroids);
  bool mat_centroids_given = (cell_mat_centroids.size() > 0);

  // Compute offsets into flattened arrays based on cell_num_mats
  
  std::vector<int> offsets(nsrccells);
  offsets[0] = 0;
  for (int i = 1; i < nsrccells; i++)
    offsets[i] = offsets[i-1] + cell_num_mats[i-1];


  // Count the number of materials in the problem and gather their IDs
  
  std::vector<int> mat_ids;
  nmats = 0;
  for (int c = 0; c < nsrccells; c++) {
    int nmats_cell = cell_num_mats[c];
    for (int m = 0; m < nmats_cell; m++) {
      bool found = false;
      for (int m2 = 0; m2 < nmats; m2++)
        if (mat_ids[m2] == cell_mat_ids[offsets[c]+m]) {
          found = true;
          break;
        }
      if (!found) {
        mat_ids.push_back(cell_mat_ids[offsets[c]+m]);
        nmats++;
      }
    }
  }

  // Spit out some information for the user
  
  if (rank == 0) {
    std::cout << "Source mesh has " << nsrccells << " cells\n";
    std::cout << "Target mesh has " << ntarcells << " cells\n";
    if (material_field_expressions.size()) {
      std::cout << " Specified fields for materials are ";
      for (auto const& expr : material_field_expressions)
        std::cout << "    " << expr << ", ";
      std::cout << "\n";
      if (material_field_expressions.size() < nmats)
        std::cout << "Not all material fields are specified. Missing ones will be set to 0\n";
    }
    std::cout << "   Interpolation order is " << interp_order << "\n";
    if (interp_order == 2)
      std::cout << "   Limiter type is " << limiter << "\n";
  }



  // Perform interface reconstruction for pretty pictures (optional)

  if (dim == 2) {  // XMOF2D works only in 2D (I know, shocking!!)
    auto interface_reconstructor =
        std::make_shared<Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
                                         Wonton::Jali_Mesh_Wrapper>>(sourceMeshWrapper);

    // convert from Portage point to Tangram point
    int ncen = cell_mat_centroids.size();
    std::vector<Tangram::Point<2>> Tcell_mat_centroids(ncen);

    // have to do this manually because cell_mat_centroids can be
    // Portage::Point<2> or Portage::Point<3> and due to the lack of
    // static_if condition, the compiler thinks it may have to convert
    // Portage::Point<3> to Tangram::Point<2>
    
    for (int i = 0; i < ncen; i++)
      for (int j = 0; j < dim; j++)
        Tcell_mat_centroids[i][j] = cell_mat_centroids[i][j];
    
    interface_reconstructor->set_volume_fractions(cell_num_mats,
                                                  cell_mat_ids,
                                                  cell_mat_volfracs,
                                                  Tcell_mat_centroids);
    interface_reconstructor->reconstruct();

    std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> const&
        cellmatpoly_list = interface_reconstructor->cell_matpoly_ptrs();

    Tangram::write_to_gmv<Wonton::Jali_Mesh_Wrapper, 2>(sourceMeshWrapper,
                                                          nmats,
                                                          cell_num_mats,
                                                          cell_mat_ids,
                                                          cellmatpoly_list,
                                                          "source_ir.gmv");
  }



  
  // Convert data from cell-centric to material-centric form as we
  // will need it for adding it to the state manager
  
  std::vector<std::vector<int>> matcells(nmats);
  std::vector<std::vector<double>> mat_volfracs(nmats);
  std::vector<std::vector<Portage::Point<dim>>> mat_centroids(nmats);
  for (int c = 0; c < nsrccells; c++) {
    int ibeg = offsets[c];
    for (int j = 0; j < cell_num_mats[c]; j++) {
      int m = cell_mat_ids[ibeg+j];
      matcells[m].push_back(c);
      mat_volfracs[m].push_back(cell_mat_volfracs[ibeg+j]);
      if (mat_centroids_given)
        mat_centroids[m].push_back(cell_mat_centroids[ibeg+j]);
    }
  }

  // Add materials, volume fractions and centroids to source state

  std::vector<std::string> matnames(nmats);
  for (int m = 0; m < nmats; m++) {
    std::stringstream matstr;
    matstr << "mat" << mat_ids[m];
    matnames[m] = matstr.str();
    sourceStateWrapper.add_material(matnames[m], matcells[m]);
  }

  for (int m = 0; m < nmats; m++) {
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &((mat_volfracs[m])[0]));
    if (mat_centroids_given)
      sourceStateWrapper.mat_add_celldata("mat_centroids", m, &((mat_centroids[m])[0]));
  }

  
  // User specified fields on source

  std::vector<user_field_t> mat_fields(nmats);
  if (material_field_expressions.size()) {
    for (int m = 0; m < nmats; m++) {
      if (!mat_fields[m].initialize(dim, material_field_expressions[m]))
        MPI_Abort(MPI_COMM_WORLD, -1);
      
      int nmatcells = matcells[m].size();
      std::vector<double> matData(nmatcells);
      for (int ic = 0; ic < nmatcells; ic++) {
        int c = matcells[m][ic];
        matData[ic] = mat_fields[m](sourceMesh->cell_centroid(c));
      }
      
      sourceStateWrapper.mat_add_celldata("cellmatdata", m, &(matData[0]));
    }
  }




  

  // Add the materials into the target mesh but with empty cell lists
  // The remap algorithm will figure out which cells contain which materials

  std::vector<int> dummylist;
  for (int m = 0; m < nmats; m++)
    targetStateWrapper.add_material(matnames[m], dummylist);

  // Add the volume fractions, centroids and cellmatdata variables
  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Portage::Point<2>>("mat_centroids");

  // Register the variable name and interpolation order with the driver
  std::vector<std::string> remap_fields;

  // If the user specified some material fields, then add a placeholder for
  // them on the target side
  if (material_field_expressions.size()) {
    targetStateWrapper.mat_add_celldata<double>("cellmatdata");
    remap_fields.push_back("cellmatdata");
  }
  
  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);

  if (dim == 2) {
    if (interp_order == 1) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR2D,
        Portage::Interpolate_1stOrder,
        2,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::XMOF2D_Wrapper>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.run(numpe > 1);
    } else if (interp_order == 2) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR2D,
        Portage::Interpolate_2ndOrder,
        2,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::XMOF2D_Wrapper>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields, limiter);
      driver.run(numpe > 1);
    }
  } else {  // 3D
    if (interp_order == 1) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_1stOrder,
        3,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::SLIC>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.run(numpe > 1);
    } else {  // 2nd order & 3D
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_2ndOrder,
        3,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::SLIC>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields, limiter);
      driver.run(numpe > 1);
    }
  }

  // Dump some timing information
  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);


  // Perform interface reconstruction on target mesh for pretty pictures
  // (optional)

  offsets.resize(ntarcells);
  offsets[0] = 0;
  for (int c = 0; c < ntarcells; c++)
    offsets[c+1] = offsets[c] + targetStateWrapper.cell_get_num_mats(c);
  int ntotal = offsets[ntarcells-1] +
      targetStateWrapper.cell_get_num_mats(ntarcells-1);

  std::vector<int> target_cell_num_mats(ntarcells, 0);
  std::vector<int> target_cell_mat_ids(ntotal);
  std::vector<double> target_cell_mat_volfracs(ntotal);
  std::vector<Tangram::Point<2>> target_cell_mat_centroids(ntotal);

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);

    double const *matvf;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf);

    Portage::Point<dim> const *matcen;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen);

    int nmatcells = matcells.size();
    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcells[ic];
      int offset = offsets[c];
      int& ncmats = target_cell_num_mats[c];
      target_cell_mat_ids[offsets[c]+ncmats] = m;
      target_cell_mat_volfracs[offsets[c]+ncmats] = matvf[ic];
      for (int i = 0; i < 2; i++)
        target_cell_mat_centroids[offsets[c]+ncmats][i] = matcen[ic][i];
      ncmats++;
    }
  }

  if (dim == 2) {  // XMOF2D works only in 2D (I know, shocking!!)
    auto interface_reconstructor =
        std::make_shared<Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
                                         Wonton::Jali_Mesh_Wrapper>>(targetMeshWrapper);

    interface_reconstructor->set_volume_fractions(target_cell_num_mats,
                                                  target_cell_mat_ids,
                                                  target_cell_mat_volfracs,
                                                  target_cell_mat_centroids);
    interface_reconstructor->reconstruct();

    std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> const&
        cellmatpoly_list = interface_reconstructor->cell_matpoly_ptrs();

    Tangram::write_to_gmv<Wonton::Jali_Mesh_Wrapper, 2>(targetMeshWrapper,
                                                        nmats,
                                                        target_cell_num_mats,
                                                        target_cell_mat_ids,
                                                        cellmatpoly_list,
                                                        "target_ir.gmv");
  }



  
  // Compute error

  if (material_field_expressions.size()) {
    double error, toterr = 0.0;
    double const * cellmatvals;
    double totvolume = 0.;
    for (int m = 0; m < nmats; m++) {
      targetStateWrapper.mat_get_celldata<double>("cellmatdata", m, &cellmatvals);
      
      std::vector<int> matcells;
      targetStateWrapper.mat_get_cells(m, &matcells);
      
      // Cell error computation
      Portage::Point<dim> ccen;
      int nmatcells = matcells.size();
      for (int ic = 0; ic < nmatcells; ++ic) {
        int c = matcells[ic];
        targetMeshWrapper.cell_centroid(c, &ccen);
        error = mat_fields[m](ccen) - cellmatvals[ic];
        
        if (!targetMeshWrapper.on_exterior_boundary(Portage::Entity_kind::CELL, c)) {
          double cellvol = targetMeshWrapper.cell_volume(c);
          totvolume += cellvol;
          *L1_error += fabs(error)*cellvol;
          *L2_error += error*error*cellvol;
        }
      }
    }
  }

  *L2_error = sqrt(*L2_error);
  if (numpe > 1) {
    std::cout << std::flush << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    double globalerr;
    MPI_Reduce(L1_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    *L1_error = globalerr;

    MPI_Reduce(L2_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    *L2_error = globalerr;
  }
  if (rank == 0) {
    std::printf("\n\nL1 NORM OF ERROR = %lf\n", *L1_error);
    std::printf("L2 NORM OF ERROR = %lf\n\n", *L2_error);
  }

  // Write out the meshes if requested
  if (mesh_output) {

    if (rank == 0)
      std::cout << "Dumping data to Exodus files..." << std::endl;
    if (material_field_expressions.size() > 0) {

      // For now convert each material field into a mesh field with
      // zero values for cells that don't have the material

      for (int m = 0; m < nmats; m++) {
        std::string varname1 = "cellmatdata_" + matnames[m];
        std::string varname2 = "cellmatdata_wtd_" + matnames[m];
        std::vector<double> cellvec(nsrccells, 0.0);
        std::vector<double> cellvec_wtd(nsrccells, 0.0);

        std::vector<int> matcells;
        sourceStateWrapper.mat_get_cells(m, &matcells);
        int nmatcells = matcells.size();

        double *matvec;
        sourceStateWrapper.mat_get_celldata("cellmatdata", m, &matvec);

        double *matvolfracs;
        sourceStateWrapper.mat_get_celldata("mat_volfracs", m, &matvolfracs);

        for (int ic = 0; ic < nmatcells; ic++) {
          int c = matcells[ic];
          cellvec[c] = matvec[ic];
          cellvec_wtd[c] = matvec[ic]*matvolfracs[ic];
        }
        sourceStateWrapper.mesh_add_data(Portage::CELL, varname1, &(cellvec[0]));
        sourceStateWrapper.mesh_add_data(Portage::CELL, varname2, &(cellvec_wtd[0]));
      }
    
      sourceState->export_to_mesh();
      sourceMesh->write_to_exodus_file("input.exo");
    }

    if (material_field_expressions.size() > 0) {
      // For now convert each material field into a mesh field with
      // zero values for cells that don't have the material
      
      for (int m = 0; m < nmats; m++) {
        std::string varname1 = "cellmatdata_" + matnames[m];
        std::string varname2 = "cellmatdata_wtd_" + matnames[m];
        std::vector<double> cellvec(ntarcells, 0.0);
        std::vector<double> cellvec_wtd(ntarcells, 0.0);

        std::vector<int> matcells;
        targetStateWrapper.mat_get_cells(m, &matcells);
        int nmatcells = matcells.size();

        double *matvec;
        targetStateWrapper.mat_get_celldata("cellmatdata", m, &matvec);

        double *matvolfracs;
        targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvolfracs);

        for (int ic = 0; ic < nmatcells; ic++) {
          int c = matcells[ic];
          cellvec[c] = matvec[ic];
          cellvec_wtd[c] = matvec[ic]*matvolfracs[ic];
        }
        targetStateWrapper.mesh_add_data(Portage::CELL, varname1, &(cellvec[0]));
        targetStateWrapper.mesh_add_data(Portage::CELL, varname2, &(cellvec_wtd[0]));
      }
    }

    targetState->export_to_mesh();
    targetMesh->write_to_exodus_file("output.exo");
    if (rank == 0)
      std::cout << "...done." << std::endl;
  }

}
