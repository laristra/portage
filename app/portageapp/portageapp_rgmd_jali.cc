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

#include "portage/support/portage.h"
#include "portage/support/mpi_collate.h"
#include "portage/driver/mmdriver.h"
#include "portage/support/timer.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "tangram/utility/get_material_moments.h"
#include "tangram/utility/rpgtools/cuts.h"
#include "tangram/utility/rpgtools/primitives.h"
#include "tangram/utility/rpgtools/examples/matdata_rotor3d.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/reconstruct/VOF.h"

#ifdef HAVE_XMOF2D
  #include "tangram/reconstruct/xmof2D_wrapper.h"
  #define IR_2D XMOF2D_Wrapper
#else
  #define IR_2D MOF
#endif

#include "portage/driver/write_to_gmv.h"

#if ENABLE_TIMINGS
  #include "portage/support/timer.h"
#endif

// For parsing and evaluating user defined expressions in apps
#include "user_field.h"

using Wonton::Jali_Mesh_Wrapper;

namespace RGMDApp {
  enum Problem_type {
    SRCPURECELLS, TJUNCTION, BALL, ROTOR
  };
}

/*! 
  @file portageapp_rgmd_jali.cc 

  @brief A simple application that remaps multi-material fields
  between two meshes - the meshes can be internally generated
  rectangular meshes or externally read unstructured meshes

  This program is used to showcase our capabilities with various types
  of remap operations (e.g. interpolation order) on various types of
  meshes (2d or 3d; node-centered or zone-centered) with multiple
  materials. The material data is runtime generated (RGMD stands for
  runtime generated material data) based on the chosen problem, which
  defines the material distribution in the computational domain.

  Currently, the following problems are implemented:
  - srcpurecells: domain contains three materials, interfaces between materials
  are linear (two lines in 2D and two planes in 3D) and ARE parallel to
  coordinate lines/planes. Interfaces forms a T-junction in the center of
  the domain. Cells that end up containing more than one material will be reset
  to contain the material with the largest volume fraction.
  - tjunction: domain contains three materials, interfaces between materials
  are linear (two lines in 2D and two planes in 3D) and are NOT parallel to
  coordinate lines/planes. Interfaces forms a T-junction in the center of
  the domain.
  - ball: domain contains two (2D) or three (3D) material, interfaces between
  material are faceted circles/spheres. In 3D, two concentric spheres are used
  to create a relatively thin shell. Both the circle and the spheres are placed
  off-center and in such a way that material interfaces touch the boundary of
  the domain.
  - rotor: only available in 3D, the domain contains 13 materials of complex 
  geometrical shapes, with some of them being thin layers.
*/

//////////////////////////////////////////////////////////////////////

int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp_rgmd_jali " <<
      "--problem=srcpurecells|tjunction|ball|rotor \n" <<  
      "--dim=2|3 --nsourcecells=N --ntargetcells=M \n" << 
      "--source_convex_cells=y|n --target_convex_cells=y|n \n" <<
      "--remap_order=1|2 \n" <<
      "--limiter=barth_jespersen --bnd_limiter=zero_gradient \n"
      "--mesh_min=0. --mesh_max=1. \n" <<
      "--output_meshes=y|n --convergence_study=NREF --only_threads=y|n \n\n";

  std::cout << "--problem (default = tjunction): defines material distribution in the domain\n\n";

  std::cout << "--dim (default = 2): spatial dimension of mesh\n\n";

  std::cout << "--nsourcecells (NO DEFAULT): Internally generated rectangular "
            << "SOURCE mesh with num cells in each coord dir\n\n";
  std::cout << "        -----  OR  ----- \n";
  std::cout << "--source_file=srcfilename: file name of source mesh (Exodus II format only)\n\n";

  std::cout << "--ntargetcells (NO DEFAULT): Internally generated rectangular "
            << "TARGET mesh with num cells in each coord dir\n\n";
  std::cout << "        -----  OR  ----- \n";
  std::cout << "--target_file=trgfilename: file name of target mesh (Exodus II format only)\n\n";

  std::cout << "--source_convex_cells (default = y): " <<
      "specifies whether all the source mesh cells are convex\n" <<
      "if mesh contains non-convex cells, all the cells will be decomposed into simplices\n" <<
      "during material data generation and interface reconstruction\n\n";

  std::cout << "--target_convex_cells (default = y): " <<
      "specifies whether all the target mesh cells are convex\n" <<
      "if mesh contains non-convex cells, all the cells will be decomposed into simplices\n" <<
      "during interface reconstruction\n\n";

  std::cout << "--mesh_min (default = 0.): " <<
      "coordinates (same in x, y, and z) of the lower corner of a mesh\n" <<
      "ONLY APPLICABLE FOR INTERNALLY GENERATED MESHES\n\n";

  std::cout << "--mesh_max (default = 1.): " <<
      "coordinates (same in x, y, and z) of the upper corner of a mesh\n" <<
      "ONLY APPLICABLE FOR INTERNALLY GENERATED MESHES\n\n";

  std::cout << "--material_fields: A comma separated list of quoted math expressions \n" <<
      " expressed in terms of \n" << "x, y and z following the syntax of " <<
      " the expression parser package ExprTk \n" <<
      "(http://www.partow.net/programming/exprtk/)\n";
  std::cout << "The syntax is generally predictable, e.g. " <<
      " \"24*x + y*y\" or \"x + log(abs(y)+1)\"\n";
  std::cout << "There must be as many expressions as there are materials in the problem\n\n";

  std::cout << "--remap order (default = 1): " <<
      "order of accuracy of interpolation\n\n";

  std::cout << "--limiter (default = NOLIMITER): " <<
      "slope limiter for a piecewise linear reconstrution\n\n";

  std::cout << "--bnd_limiter (default = NOLIMITER): " <<
      "slope limiter on the boundary for a piecewise linear reconstruction\n\n";

  std::cout << "--convergence_study (default = 1): provide the number of times "
            << "you want to double source and target mesh sizes \n";
  std::cout << "  ONLY APPLICABLE IF BOTH MESHES ARE INTERNALLY GENERATED\n\n";

  std::cout << "--output_meshes (default = y)\n";
  std::cout << "  If 'y', the source and target meshes are output with the " <<
      "remapped field attached as input.exo and output.exo, \n" <<
      "and interface reconstruction results are output as source_mm.gmv and target_mm.gmv \n\n";
  
#if ENABLE_TIMINGS
  std::cout << "--only_threads (default = n)\n";
  std::cout << " enable if you want to profile only threads scaling\n\n";

  std::cout << "--scaling (default = strong)\n";
  std::cout << " specify the scaling study type [strong|weak]\n\n";
#endif

  return 0;
}
//////////////////////////////////////////////////////////////////////

template<class Mesh_Wrapper>
void srcpurecells_material_data(const Mesh_Wrapper& mesh,
                                std::vector<int>& mesh_material_IDs,
                                std::vector<int>& cell_num_mats,
                                std::vector<int>& cell_mat_ids,
                                std::vector<double>& cell_mat_volfracs,
                                std::vector< Wonton::Point<2> >& cell_mat_centroids,
                                const double vol_tol,
                                const double dst_tol,
                                bool decompose_cells) {
  mesh_material_IDs = {2, 0, 1};
  const std::vector< Wonton::Vector<2> > material_interface_normals = {
    Wonton::Vector<2>(1.0, 0.0), Wonton::Vector<2>(0.0, -1.0)
  };
  const std::vector< Wonton::Point<2> > material_interface_points = {
    Wonton::Point<2>(0.5, 0.5), Wonton::Point<2>(0.5, 0.5)
  };

  int nmesh_materials = static_cast<int>(mesh_material_IDs.size());
  std::vector< Tangram::Plane_t<2> > material_interfaces(nmesh_materials - 1);
  for (int iline = 0; iline < nmesh_materials - 1; iline++) {
    material_interfaces[iline].normal = material_interface_normals[iline];
    material_interfaces[iline].dist2origin =
      -Wonton::dot(material_interface_points[iline].asV(),
                   material_interfaces[iline].normal);
  }  

  get_material_moments<Mesh_Wrapper>(mesh, material_interfaces, mesh_material_IDs, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    vol_tol, dst_tol, decompose_cells);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  assert(ncells == cell_num_mats.size());
  if (ncells != cell_mat_ids.size()) {
    int nextras = cell_mat_ids.size() - ncells;
    for (int icell = 0; icell < ncells; icell++)
      if (cell_num_mats[icell] > 1) {
        int ncellmats = cell_num_mats[icell];
        int idominant_mat = cell_mat_ids[std::distance(cell_mat_volfracs.begin(),
          std::max_element(cell_mat_volfracs.begin() + icell,
                           cell_mat_volfracs.begin() + icell + ncellmats))];

        cell_num_mats[icell] = 1;
        cell_mat_ids[icell] = idominant_mat;
        cell_mat_volfracs[icell] = 1.0;
        mesh.cell_centroid(icell, &cell_mat_centroids[icell]);
        cell_mat_ids.erase(cell_mat_ids.begin() + icell + 1, 
                           cell_mat_ids.begin() + icell + ncellmats);
        cell_mat_volfracs.erase(cell_mat_volfracs.begin() + icell + 1, 
                                cell_mat_volfracs.begin() + icell + ncellmats);
        cell_mat_centroids.erase(cell_mat_centroids.begin() + icell + 1, 
                                 cell_mat_centroids.begin() + icell + ncellmats);

        nextras -= ncellmats - 1;
        if (nextras == 0) break;
      }
  }
}

template<class Mesh_Wrapper>
void tjunction_material_data(const Mesh_Wrapper& mesh,
                             std::vector<int>& mesh_material_IDs,
                             std::vector<int>& cell_num_mats,
                             std::vector<int>& cell_mat_ids,
                             std::vector<double>& cell_mat_volfracs,
                             std::vector< Wonton::Point<2> >& cell_mat_centroids,
                             const double vol_tol,
                             const double dst_tol,
                             bool decompose_cells) {
  mesh_material_IDs = {2, 0, 1};
  const std::vector< Wonton::Vector<2> > material_interface_normals = {
    Wonton::Vector<2>(0.5, 0.5), Wonton::Vector<2>(0.5, -0.375)
  };
  const std::vector< Wonton::Point<2> > material_interface_points = {
    Wonton::Point<2>(0.5, 0.5), Wonton::Point<2>(0.5, 0.5)
  };

  int nmesh_materials = static_cast<int>(mesh_material_IDs.size());
  std::vector< Tangram::Plane_t<2> > material_interfaces(nmesh_materials - 1);
  for (int iline = 0; iline < nmesh_materials - 1; iline++) {
    material_interfaces[iline].normal = material_interface_normals[iline];
    material_interfaces[iline].normal.normalize();
    material_interfaces[iline].dist2origin =
      -Wonton::dot(material_interface_points[iline].asV(),
                   material_interfaces[iline].normal);
  }  

  get_material_moments<Mesh_Wrapper>(mesh, material_interfaces, mesh_material_IDs, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    vol_tol, dst_tol, decompose_cells);  
}

template<class Mesh_Wrapper>
void ball_material_data(const Mesh_Wrapper& mesh,
                        std::vector<int>& mesh_material_IDs,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector< Wonton::Point<2> >& cell_mat_centroids,
                        const double vol_tol,
                        const double dst_tol,
                        bool decompose_cells) {
  mesh_material_IDs = {1, 0};
  Wonton::Point<2> circle_cen(0.75, 0.75);
  double circle_rad = 0.25;
  int nquadrant_samples = 90;

  get_material_moments<Mesh_Wrapper>(mesh, mesh_material_IDs, 
    circle_cen, circle_rad, nquadrant_samples,
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids, 
    vol_tol, dst_tol, decompose_cells);  
}

template<class Mesh_Wrapper>
void rotor_material_data(const Mesh_Wrapper& mesh,
                         std::vector<int>& mesh_material_IDs,
                         std::vector<int>& cell_num_mats,
                         std::vector<int>& cell_mat_ids,
                         std::vector<double>& cell_mat_volfracs,
                         std::vector< Wonton::Point<2> >& cell_mat_centroids,
                         const double vol_tol,
                         const double dst_tol,
                         bool decompose_cells) {
  std::cerr << "Rotor problem is only available in 3D!\n";
  MPI_Abort(MPI_COMM_WORLD, -1);
}

template<class Mesh_Wrapper>
void srcpurecells_material_data(const Mesh_Wrapper& mesh,
                                std::vector<int>& mesh_material_IDs,
                                std::vector<int>& cell_num_mats,
                                std::vector<int>& cell_mat_ids,
                                std::vector<double>& cell_mat_volfracs,
                                std::vector< Wonton::Point<3> >& cell_mat_centroids,
                                const double vol_tol,
                                const double dst_tol,
                                bool decompose_cells) {
  mesh_material_IDs = {2, 0, 1};
  const std::vector< Wonton::Vector<3> > material_interface_normals = {
    Wonton::Vector<3>(1.0, 0.0, 0.0), Wonton::Vector<3>(0.0, 0.0, -1.0)
  };
  const std::vector< Wonton::Point<3> > material_interface_points = {
    Wonton::Point<3>(0.5, 0.5, 0.5), Wonton::Point<3>(0.5, 0.5, 0.5)
  };

  int nmesh_materials = static_cast<int>(mesh_material_IDs.size());
  std::vector< Tangram::Plane_t<3> > material_interfaces(nmesh_materials - 1);
  for (int iplane = 0; iplane < nmesh_materials - 1; iplane++) {
    material_interfaces[iplane].normal = material_interface_normals[iplane];
    material_interfaces[iplane].dist2origin =
      -Wonton::dot(material_interface_points[iplane].asV(),
                   material_interfaces[iplane].normal);
  }

  get_material_moments<Mesh_Wrapper>(mesh, material_interfaces, mesh_material_IDs, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    vol_tol, dst_tol, decompose_cells);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  assert(ncells == cell_num_mats.size());
  if (ncells != cell_mat_ids.size()) {
    int nextras = cell_mat_ids.size() - ncells;
    for (int icell = 0; icell < ncells; icell++)
      if (cell_num_mats[icell] > 1) {
        int ncellmats = cell_num_mats[icell];
        int idominant_mat = cell_mat_ids[std::distance(cell_mat_volfracs.begin(),
          std::max_element(cell_mat_volfracs.begin() + icell,
                           cell_mat_volfracs.begin() + icell + ncellmats))];

        cell_num_mats[icell] = 1;
        cell_mat_ids[icell] = idominant_mat;
        cell_mat_volfracs[icell] = 1.0;
        mesh.cell_centroid(icell, &cell_mat_centroids[icell]);
        cell_mat_ids.erase(cell_mat_ids.begin() + icell + 1, 
                           cell_mat_ids.begin() + icell + ncellmats);
        cell_mat_volfracs.erase(cell_mat_volfracs.begin() + icell + 1, 
                                cell_mat_volfracs.begin() + icell + ncellmats);
        cell_mat_centroids.erase(cell_mat_centroids.begin() + icell + 1, 
                                 cell_mat_centroids.begin() + icell + ncellmats);

        nextras -= ncellmats - 1;
        if (nextras == 0) break;
      }
  }
}

template<class Mesh_Wrapper>
void tjunction_material_data(const Mesh_Wrapper& mesh,
                             std::vector<int>& mesh_material_IDs,
                             std::vector<int>& cell_num_mats,
                             std::vector<int>& cell_mat_ids,
                             std::vector<double>& cell_mat_volfracs,
                             std::vector< Wonton::Point<3> >& cell_mat_centroids,
                             const double vol_tol,
                             const double dst_tol,                             
                             bool decompose_cells) {
  mesh_material_IDs = {2, 0, 1};
  const std::vector< Wonton::Vector<3> > material_interface_normals = {
    Wonton::Vector<3>(0.5, 0.5, 0.0), Wonton::Vector<3>(0.5, 0.025, -0.375)
  };
  const std::vector< Wonton::Point<3> > material_interface_points = {
    Wonton::Point<3>(0.5, 0.5, 0.5), Wonton::Point<3>(0.5, 0.5, 0.5)
  };

  int nmesh_materials = static_cast<int>(mesh_material_IDs.size());
  std::vector< Tangram::Plane_t<3> > material_interfaces(nmesh_materials - 1);
  for (int iplane = 0; iplane < nmesh_materials - 1; iplane++) {
    material_interfaces[iplane].normal = material_interface_normals[iplane];
    material_interfaces[iplane].normal.normalize();
    material_interfaces[iplane].dist2origin =
      -Wonton::dot(material_interface_points[iplane].asV(),
                   material_interfaces[iplane].normal);
  }

  get_material_moments<Mesh_Wrapper>(mesh, material_interfaces, mesh_material_IDs, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    vol_tol, dst_tol, decompose_cells);  
}

template<class Mesh_Wrapper>
void ball_material_data(const Mesh_Wrapper& mesh,
                        std::vector<int>& mesh_material_IDs,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector< Wonton::Point<3> >& cell_mat_centroids,
                        const double vol_tol,
                        const double dst_tol,
                        bool decompose_cells) {
  mesh_material_IDs = {2, 0, 1};
  Wonton::Point<3> spheres_cen(0.75, 0.75, 0.75);
  std::vector<double> sphere_rads = {0.2495, 0.25};
  int nquadrant_samples = 24;

  get_material_moments<Mesh_Wrapper>(mesh, mesh_material_IDs, 
    spheres_cen, sphere_rads, nquadrant_samples,
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids, 
    vol_tol, dst_tol, decompose_cells);  
}

template<class Mesh_Wrapper>
void rotor_material_data(const Mesh_Wrapper& mesh,
                         std::vector<int>& mesh_material_IDs,
                         std::vector<int>& cell_num_mats,
                         std::vector<int>& cell_mat_ids,
                         std::vector<double>& cell_mat_volfracs,
                         std::vector< Wonton::Point<3> >& cell_mat_centroids,
                         const double vol_tol,
                         const double dst_tol,
                         bool decompose_cells) {
  std::vector< std::string > mesh_material_names;

  rotor_material_moments(mesh, mesh_material_IDs, mesh_material_names, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids, 
    vol_tol, dst_tol, decompose_cells);
}

// Generic interface reconstructor factory

template<int dim, class MeshWrapper>
class interface_reconstructor_factory {};
  
// Specializations
template<class MeshWrapper>
class interface_reconstructor_factory<2, MeshWrapper>{
 public:
  interface_reconstructor_factory(MeshWrapper const& mesh,
                                  std::vector<Tangram::IterativeMethodTolerances_t> tols,
                                  bool all_convex) :
      mesh_(mesh), tols_(tols), all_convex_(all_convex) {};

  auto operator()() -> decltype(auto) {
    return std::make_shared<Tangram::Driver<Tangram::IR_2D, 2, MeshWrapper,
                                            Tangram::SplitR2D,
                                            Tangram::ClipR2D>>(mesh_, tols_, all_convex_);
  }

 private:
  MeshWrapper const& mesh_;
  std::vector<Tangram::IterativeMethodTolerances_t> tols_;
  bool all_convex_;
};

template<class MeshWrapper>
class interface_reconstructor_factory<3, MeshWrapper>{
 public:
  interface_reconstructor_factory(MeshWrapper const& mesh,
                                  std::vector<Tangram::IterativeMethodTolerances_t> tols,
                                  bool all_convex) :
      mesh_(mesh), tols_(tols), all_convex_(all_convex) {};

  auto operator()() -> decltype(auto) {
    return std::make_shared<Tangram::Driver<Tangram::MOF, 3, MeshWrapper,
                                            Tangram::SplitR3D,
                                            Tangram::ClipR3D>>(mesh_, tols_, all_convex_);
  }

 private:
  MeshWrapper const& mesh_;
  std::vector<Tangram::IterativeMethodTolerances_t> tols_;
  bool all_convex_;
};

// Forward declaration of function to run remap on two meshes and
// return the L1 and L2 error norm in the remapped field w.r.t. to an
// analytically imposed field. If no field was imposed, the errors are
// returned as 0

template<int dim> void run(std::shared_ptr<Jali::Mesh> sourceMesh,
                           std::shared_ptr<Jali::Mesh> targetMesh,
                           bool source_convex_cells,
                           bool target_convex_cells,
                           Portage::Limiter_type limiter,
                           Portage::Boundary_Limiter_type bnd_limiter,
                           int interp_order,
                           RGMDApp::Problem_type problem,
                           std::vector<std::string> material_field_expressions,
                           std::string field_filename,
                           bool mesh_output, int rank, int numpe,
                           Jali::Entity_kind entityKind,
                           double& L1_error, double& L2_error,
                           std::shared_ptr<Profiler> profiler = nullptr);

// Forward declaration of function to run interface reconstruction and
// write the material polygons and their associated fields to a GMV file

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

  RGMDApp::Problem_type problem = RGMDApp::Problem_type::TJUNCTION;
  int nsourcecells = 0, ntargetcells = 0;  // No default
  int dim = 2;
  bool source_convex_cells = true, target_convex_cells = true;
  std::vector<std::string> material_field_expressions;
  std::string srcfile, trgfile;  // No default
  std::string field_output_filename;  // No default;

  int interp_order = 1;
  bool mesh_output = true;
  int n_converge = 1;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  Portage::Limiter_type limiter = Portage::Limiter_type::NOLIMITER;
  Portage::Boundary_Limiter_type bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;
  double srclo = 0.0, srchi = 1.0;  // bounds of generated mesh in each dir

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

    if (keyword == "problem") {
      if(valueword == "srcpurecells")
        problem = RGMDApp::Problem_type::SRCPURECELLS;
      else if(valueword == "tjunction")
        problem = RGMDApp::Problem_type::TJUNCTION;
      else if (valueword == "ball")
        problem = RGMDApp::Problem_type::BALL;
      else if (valueword == "rotor")
        problem = RGMDApp::Problem_type::ROTOR;
      else {
        std::cerr << "Unknown problem type!\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    } else if (keyword == "dim") {
      dim = stoi(valueword);
      assert(dim == 2 || dim == 3);
    } else if (keyword == "nsourcecells")
      nsourcecells = stoi(valueword);
    else if (keyword == "ntargetcells")
      ntargetcells = stoi(valueword);
    else if (keyword == "source_convex_cells")
      source_convex_cells = (valueword == "y");
    else if (keyword == "target_convex_cells")
      target_convex_cells = (valueword == "y");
    else if (keyword == "source_file")
      srcfile = valueword;
    else if (keyword == "target_file")
      trgfile = valueword;
    else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
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
        limiter = Portage::Limiter_type::BARTH_JESPERSEN;
    } else if (keyword == "bnd_limiter") {
      if (valueword == "zero_gradient" || valueword == "ZERO_GRADIENT")
        bnd_limiter = Portage::Boundary_Limiter_type::BND_ZERO_GRADIENT;
      else if (valueword == "barth_jespersen" || valueword == "BARTH_JESPERSEN")
        bnd_limiter = Portage::Boundary_Limiter_type::BND_BARTH_JESPERSEN;
    } else if (keyword == "mesh_min") {
      srclo = stof(valueword);
    } else if (keyword == "mesh_max") {
      srchi = stof(valueword);
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
    } else if (keyword == "convergence_study") {
      n_converge = stoi(valueword);
      if (n_converge <= 0) {
        std::cerr << "Number of meshes for convergence study should be greater than 0" << std::endl;
        throw std::exception();
      }
#if ENABLE_TIMINGS
      assert(n_converge == 1);
#endif
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
    } else {
      std::cerr << "Unrecognized option " << keyword << " !\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
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
  if (nsourcecells > 0)
    if (material_field_expressions.size() == 0) {
      std::cout << "No field imposed on internally generated source mesh\n";
      std::cout << "Nothing to remap. Exiting...";
      print_usage();
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

#if ENABLE_TIMINGS
  auto profiler = std::make_shared<Profiler>();
  // save params for after
  profiler->params.ranks   = numpe;
  profiler->params.nsource = std::pow(nsourcecells, dim);
  profiler->params.ntarget = std::pow(ntargetcells, dim);
  profiler->params.order   = interp_order;
  profiler->params.nmats   = material_field_expressions.size();
  switch(problem) {
    case RGMDApp::Problem_type::SRCPURECELLS: profiler->params.output  = "srcpurecells_"; break;
    case RGMDApp::Problem_type::TJUNCTION: profiler->params.output  = "t-junction_"; break;
    case RGMDApp::Problem_type::BALL: profiler->params.output  = "ball_"; break;
    case RGMDApp::Problem_type::ROTOR: profiler->params.output  = "rotor_"; break;
    default: std::cerr << "Unknown problem type!\n"; MPI_Abort(MPI_COMM_WORLD, -1);
  }
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
  std::shared_ptr<Jali::Mesh> source_mesh, target_mesh;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});

  double trglo = srclo, trghi = srchi;  // bounds of generated mesh in each dir

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

#if ENABLE_TIMINGS
    profiler->time.mesh_init = timer::elapsed(tic);

    if (rank == 0) {
      if (n_converge == 1)
        std::cout << "Mesh Initialization Time: " << profiler->time.mesh_init << std::endl;
      else
        std::cout << "Mesh Initialization Time (Iteration i): " <<
            profiler->time.mesh_init << std::endl;
    }
#endif

    // Make sure we have the right dimension and that source and
    // target mesh dimensions match (important when one or both of the
    // meshes are read in)
    assert(source_mesh->space_dimension() == target_mesh->space_dimension());
    dim = source_mesh->space_dimension();

    double error_L1 = 0.0, error_L2 = 0.0;

    // Now run the remap on the meshes and get back the L2 error
    switch (dim) {
      case 2:
        run<2>(source_mesh, target_mesh, source_convex_cells, target_convex_cells,
               limiter, bnd_limiter, interp_order,
               problem, material_field_expressions,
               field_output_filename, mesh_output,
               rank, numpe, entityKind, l1_err[i], l2_err[i], profiler);
        break;
      case 3:
        run<3>(source_mesh, target_mesh, source_convex_cells, target_convex_cells,
               limiter, bnd_limiter, interp_order,
               problem, material_field_expressions,
               field_output_filename, mesh_output,
               rank, numpe, entityKind, l1_err[i], l2_err[i], profiler);
        break;
      default:
        std::cerr << "Dimension not 2 or 3" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    std::cout << "L1 norm of error for iteration " << i << " is " <<
        l1_err[i] << std::endl;
    std::cout << "L2 norm of error for iteration " << i << " is " <<
        l2_err[i] << std::endl;

    // if convergence study, double the mesh resolution
    nsourcecells *= 2;
    ntargetcells *= 2;
  }

  for (int i = 1; i < n_converge; i++) {
    std::cout << "Error ratio L1(" << i - 1 << ")/L1(" << i << ") is " <<
        l1_err[i - 1]/l1_err[i] << std::endl;
    std::cout << "Error ratio L2(" << i - 1 << ")/L2(" << i << ") is " <<
        l2_err[i - 1]/l2_err[i] << std::endl;
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

// Run a remap between two meshes and return the L1 and L2 error norms
// with respect to the specified field. If a field was not specified and
// remap only volume fractions (and if specified, centroids)

template<int dim> void run(std::shared_ptr<Jali::Mesh> sourceMesh,
                           std::shared_ptr<Jali::Mesh> targetMesh,
                           bool source_convex_cells,
                           bool target_convex_cells,
                           Portage::Limiter_type limiter,
                           Portage::Boundary_Limiter_type bnd_limiter,
                           int interp_order,
                           RGMDApp::Problem_type problem,
                           std::vector<std::string> material_field_expressions,
                           std::string field_filename, bool mesh_output,
                           int rank, int numpe, Jali::Entity_kind entityKind,
                           double& L1_error, double& L2_error,
                           std::shared_ptr<Profiler> profiler) {
  if (rank == 0)
    std::cout << "starting portageapp_rgmd_jali...\n";

  // Wrappers for interfacing with the underlying mesh data structures.
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  const int nsrccells = sourceMeshWrapper.num_owned_cells() +
    sourceMeshWrapper.num_ghost_cells();
  const int ntarcells = targetMeshWrapper.num_owned_cells() +
    targetMeshWrapper.num_ghost_cells();
  const int ntar_owned_cells = targetMeshWrapper.num_owned_cells();

  // Native jali state managers for source and target
  std::shared_ptr<Jali::State> sourceState(Jali::State::create(sourceMesh));
  std::shared_ptr<Jali::State> targetState(Jali::State::create(targetMesh));

  // Portage wrappers for source and target fields
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  // Get volume fraction and centroid data for the specified problem

  double dst_tol = sqrt(dim)*std::numeric_limits<double>::epsilon();
  double vol_tol = std::numeric_limits<double>::epsilon();

  std::vector<int> global_material_IDs;
  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;  // flattened 2D array
  std::vector<double> cell_mat_volfracs;  // flattened 2D array
  std::vector<Wonton::Point<dim>> cell_mat_centroids;  // flattened 2D array

  switch(problem) {
    case RGMDApp::Problem_type::SRCPURECELLS:
      srcpurecells_material_data<Wonton::Jali_Mesh_Wrapper>(sourceMeshWrapper,
                                                            global_material_IDs,
                                                            cell_num_mats,
                                                            cell_mat_ids,
                                                            cell_mat_volfracs,
                                                            cell_mat_centroids,
                                                            vol_tol,
                                                            dst_tol,
                                                            !source_convex_cells);
      break;
    case RGMDApp::Problem_type::TJUNCTION:
      tjunction_material_data<Wonton::Jali_Mesh_Wrapper>(sourceMeshWrapper,
                                                         global_material_IDs,
                                                         cell_num_mats,
                                                         cell_mat_ids,
                                                         cell_mat_volfracs,
                                                         cell_mat_centroids,
                                                         vol_tol,
                                                         dst_tol,
                                                         !source_convex_cells);
      break;
    case RGMDApp::Problem_type::BALL:
      ball_material_data<Wonton::Jali_Mesh_Wrapper>(sourceMeshWrapper,
                                                    global_material_IDs,
                                                    cell_num_mats,
                                                    cell_mat_ids,
                                                    cell_mat_volfracs,
                                                    cell_mat_centroids,
                                                    vol_tol,
                                                    dst_tol,
                                                    !source_convex_cells);
      break;
    case RGMDApp::Problem_type::ROTOR:
      rotor_material_data<Wonton::Jali_Mesh_Wrapper>(sourceMeshWrapper,
                                                     global_material_IDs,
                                                     cell_num_mats,
                                                     cell_mat_ids,
                                                     cell_mat_volfracs,
                                                     cell_mat_centroids,
                                                     vol_tol,
                                                     dst_tol,
                                                     !source_convex_cells);
      break;
    default:
      std::cerr << "Unknown problem type!\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int nglobal_mats = static_cast<int>(global_material_IDs.size());
  if (std::set<int>(global_material_IDs.begin(), global_material_IDs.end()).size() !=
      nglobal_mats) {
    std::cerr << "Generated materials had repeated indices, " << 
                 "check the implementation of the problem!\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if ( (*std::min_element(global_material_IDs.begin(), 
                          global_material_IDs.end()) != 0) ||
       (*std::max_element(global_material_IDs.begin(), 
                          global_material_IDs.end()) != nglobal_mats - 1) ) {
    std::cerr << "Generated materials didn't have contiguous indices, " << 
                 "check if there are material subdomains not covered by the mesh!\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (material_field_expressions.size() != nglobal_mats) {
    std::cerr << "Number of imposed fields (" << material_field_expressions.size() <<
                 ") is not equal to the number of materials in the problem (" <<
                 nglobal_mats << ")!\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // Count the number of local materials and gather their IDs
  std::set<int> mat_ids;
  for (auto id: cell_mat_ids) mat_ids.insert(id);
  std::vector<int> local_material_IDs(mat_ids.begin(), mat_ids.end());
  int nlocal_mats = static_cast<int>(local_material_IDs.size());
 
  // Compute offsets into flattened arrays based on cell_num_mats
  std::vector<int> offsets(nsrccells, 0);
  for (int i = 0; i < nsrccells - 1; i++)
    offsets[i + 1] = offsets[i] + cell_num_mats[i];

  // Output some information for the user
  if (rank == 0) {
    std::cout << "Source mesh has " << nsrccells << " cells\n";
    std::cout << "Target mesh has " << ntarcells << " cells\n";
    if (material_field_expressions.size()) {
      std::cout << " Specified fields for materials are ";
      for (auto const& expr : material_field_expressions)
        std::cout << "    " << expr << ", ";
      std::cout << "\n";
    }
    std::cout << "   Interpolation order is " << interp_order << "\n";
    if (interp_order == 2) {
      std::cout << "   Limiter type is " << limiter << "\n";
      std::cout << "   Boundary limiter type is " << bnd_limiter << "\n";
    }
  }

#if ENABLE_TIMINGS
  auto tic = timer::now();
#endif

  // Convert data from cell-centric to material-centric form as we
  // will need it for adding it to the state manager

  std::unordered_map<int, std::vector<int>> matcells;
  std::unordered_map<int, std::vector<double>> mat_volfracs;
  std::unordered_map<int, std::vector<Wonton::Point<dim>>> mat_centroids;
  for (int c = 0; c < nsrccells; c++) {
    int ibeg = offsets[c];
    for (int j = 0; j < cell_num_mats[c]; j++) {
      int mat_id = cell_mat_ids[ibeg + j];
      matcells[mat_id].push_back(c);
      mat_volfracs[mat_id].push_back(cell_mat_volfracs[ibeg + j]);
      mat_centroids[mat_id].push_back(cell_mat_centroids[ibeg + j]);
    }
  }

  // Add materials, volume fractions, and centroids to the source state
  // If the material isn't found in the partition, add empty data to the source state
  // Note names are dimensioned on all materials whether in the partition or not
  std::vector<std::string> matnames(nglobal_mats);
  for (int m = 0; m < nglobal_mats; m++) {
    matnames[m] = "mat" + std::to_string(m);
    if (matcells.find(m) != matcells.end()) {
      sourceStateWrapper.add_material(matnames[m], matcells[m]);
      sourceStateWrapper.mat_add_celldata("mat_volfracs", m, mat_volfracs[m].data());
      sourceStateWrapper.mat_add_celldata("mat_centroids", m, mat_centroids[m].data());
    }
    else {
      sourceStateWrapper.add_material(matnames[m], {});
      // need to cast the empty data so that type_id works
      sourceStateWrapper.mat_add_celldata("mat_volfracs", m, static_cast<double*>(nullptr));
      sourceStateWrapper.mat_add_celldata("mat_centroids", m, static_cast<Wonton::Point<dim>*>(nullptr));
    }
  }

#if ENABLE_TIMINGS
  profiler->time.remap = timer::elapsed(tic, true);
#endif

  // Executor
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  Wonton::Executor_type *executor = (numpe > 1) ? &mpiexecutor : nullptr;
  
  // Perform interface reconstruction

  std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {1000, dst_tol, vol_tol};
  ims_tols[1] = {100, 1.0e-15, 1.0e-15};

  interface_reconstructor_factory<dim, Wonton::Jali_Mesh_Wrapper> 
    source_IRFactory(sourceMeshWrapper, ims_tols, source_convex_cells);
  auto source_interface_reconstructor = source_IRFactory();

  source_interface_reconstructor->set_volume_fractions(cell_num_mats,
                                                       cell_mat_ids,
                                                       cell_mat_volfracs,
                                                       cell_mat_centroids);
  source_interface_reconstructor->reconstruct();

#if ENABLE_TIMINGS
  profiler->time.interface = timer::elapsed(tic, true);
#endif

  // Material fields are evaluated at material polygon centroids.
  // This allows us to test for near zero error in cases where an
  // interface reconstruction method can reconstruct the interface
  // exactly (e.g. MOF with linear interfaces) and the remapping
  // method can reproduce a field exactly (linear fields with a 2nd
  // order accurate method)

  // User specified fields on source
  std::vector<user_field_t> mat_fields(nglobal_mats);  

  if (material_field_expressions.size()) {
    for (int m = 0; m < nglobal_mats; m++) {
      if (!mat_fields[m].initialize(dim, material_field_expressions[m]))
        MPI_Abort(MPI_COMM_WORLD, -1);

      int nmatcells = matcells[m].size();
      std::vector<double> matData(nmatcells);
      for (int ic = 0; ic < nmatcells; ic++) {
        int c = matcells[m][ic];
        if (cell_num_mats[c] == 1) {
          Wonton::Point<dim> ccen;
          sourceMeshWrapper.cell_centroid(c, &ccen);
          matData[ic] = mat_fields[m](ccen);
        } else {
          Tangram::CellMatPoly<dim> const& cellmatpoly =
            source_interface_reconstructor->cell_matpoly_data(c);
          int nmp = cellmatpoly.num_matpolys();
          for (int i = 0; i < nmp; i++) {
            if (cellmatpoly.matpoly_matid(i) == m) {
              Wonton::Point<dim> mcen = cellmatpoly.matpoly_centroid(i);
              matData[ic] = mat_fields[m](mcen);
            }
          }
        }
      }

      sourceStateWrapper.mat_add_celldata("cellmatdata", m, &(matData[0]));
    }
  }

  if (rank == 0) {
    std::cout << "***Registered fieldnames:\n";
    for (auto & field_name: sourceStateWrapper.names()) 
      std::cout << " registered fieldname: " << field_name << std::endl;
  }

  std::vector<std::string> fieldnames;
  fieldnames.push_back("cellmatdata");
  
  // Add the materials into the target mesh but with empty cell lists
  // The remap algorithm will figure out which cells contain which materials

  for (int m = 0; m < nglobal_mats; m++)
    targetStateWrapper.add_material(matnames[m], {});
    
  // Add the volume fractions, centroids and cellmatdata variables
  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Wonton::Point<dim>>("mat_centroids");

  // Register the variable name and interpolation order with the driver
  std::vector<std::string> remap_fields;

  // If the user specified some material fields, then add a placeholder for
  // them on the target side
  if (material_field_expressions.size()) {
    targetStateWrapper.mat_add_celldata<double>("cellmatdata");
    remap_fields.push_back("cellmatdata");
  }

  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);


  if (mesh_output) {  
    std::string filename = "source_mm_" + std::to_string(rank) + ".gmv";
    Portage::write_to_gmv<dim>(sourceMeshWrapper, sourceStateWrapper,
                               source_interface_reconstructor, fieldnames,
                               filename);
  }

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
        Tangram::IR_2D, 
        Tangram::SplitR2D,
        Tangram::ClipR2D>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_reconstructor_options(ims_tols, source_convex_cells);
      driver.run(executor);
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
        Tangram::IR_2D,
        Tangram::SplitR2D,
        Tangram::ClipR2D>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_limiter(limiter);
      driver.set_bnd_limiter(bnd_limiter);
      driver.set_reconstructor_options(ims_tols, source_convex_cells);
      driver.run(executor);
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
        Tangram::MOF,
        Tangram::SplitR3D,
        Tangram::ClipR3D>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_reconstructor_options(ims_tols, source_convex_cells);
      driver.run(executor);
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
        Tangram::MOF,
        Tangram::SplitR3D,
        Tangram::ClipR3D>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_limiter(limiter);
      driver.set_bnd_limiter(bnd_limiter);
      driver.set_reconstructor_options(ims_tols, source_convex_cells);
      driver.run(executor);
    }
  }

  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);

  // Convert data from material-centric to cell-centric form

  offsets.resize(ntarcells, 0);
  std::vector<int> owned2all(ntar_owned_cells, -1);
  for (int i = 0, iowned = 0; i < ntarcells - 1; i++) {
    // Target state has no information on ghost cells
    if (targetMeshWrapper.cell_get_type(i) == Wonton::PARALLEL_GHOST)
      offsets[i + 1] = offsets[i];
    else {
      offsets[i + 1] = offsets[i] + targetStateWrapper.cell_get_num_mats(iowned);
      owned2all[iowned] = i;
      iowned++;
    }
  }

  int mat_vec_len = offsets[ntarcells - 1];
  if (targetMeshWrapper.cell_get_type(ntarcells - 1) == Wonton::PARALLEL_OWNED) {
    mat_vec_len += targetStateWrapper.cell_get_num_mats(ntar_owned_cells - 1);
    owned2all[ntar_owned_cells - 1] = ntarcells - 1;
  }

  std::vector<int> target_cell_num_mats(ntarcells, 0);
  std::vector<int> target_cell_mat_ids(mat_vec_len);
  std::vector<double> target_cell_mat_volfracs(mat_vec_len);
  std::vector<Tangram::Point<dim>> target_cell_mat_centroids(mat_vec_len);

  // We only set the material data for owned cells
  for (int m = 0; m < nglobal_mats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);

    double const *matvf;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf);

    Wonton::Point<dim> const *matcen;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen);

    int nmatcells = matcells.size();
    for (int ic = 0; ic < nmatcells; ic++) {
      int state_cell_id = matcells[ic];
      int mesh_cell_id = owned2all[state_cell_id];
      int offset = offsets[mesh_cell_id];
      int& ncmats = target_cell_num_mats[mesh_cell_id];
      target_cell_mat_ids[offsets[mesh_cell_id] + ncmats] = m;
      target_cell_mat_volfracs[offsets[mesh_cell_id] + ncmats] = matvf[ic];
      for (int i = 0; i < dim; i++)
        target_cell_mat_centroids[offsets[mesh_cell_id] + ncmats][i] = matcen[ic][i];
      ncmats++;
    }
  }

#if !ENABLE_TIMINGS
  // Output some information for the user
  for (int m = 0; m < nglobal_mats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);

    double const *matvf;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf);

    Wonton::Point<dim> const *matcen;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen);

    double const *cellmatdata;
    targetStateWrapper.mat_get_celldata("cellmatdata", m, &cellmatdata);

    int nmatcells = matcells.size();
    
    std::cout << "\n----target owned cell global indices on rank " << rank << " for material " << m << ": ";
    for (int ic = 0; ic < nmatcells; ic++) std::cout <<
      targetMeshWrapper.get_global_id(owned2all[matcells[ic]], Wonton::Entity_kind::CELL) << " ";
    std::cout << std::endl; 
    
    std::cout << "----mat_volfracs on rank " << rank << " for material " << m << ": ";
    for (int ic = 0; ic < nmatcells; ic++) std::cout << matvf[ic] << " ";
    std::cout << std::endl; 
    
    std::cout << "----mat_centroids on rank " << rank<< " for material " << m << ": ";
    for (int ic = 0; ic < nmatcells; ic++) std::cout << "(" << matcen[ic][0]<< ", " << matcen[ic][1] << ") ";
    std::cout << std::endl; 
    
    std::cout << "----cellmatdata on rank " << rank << " for material " << m << ": ";
    for (int ic = 0; ic < nmatcells; ic++) std::cout << cellmatdata[ic] << " ";
    std::cout << std::endl << std::endl; 
  }
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

    if (material_field_expressions.size() > 0) {
      // For now convert each material field into a mesh field with
      // zero values for cells that don't have the material

      for (int m = 0; m < nlocal_mats; m++) {
        int mat_id = local_material_IDs[m];

        std::string varname1 = "cellmatdata_" + matnames[mat_id];
        std::string varname2 = "cellmatdata_wtd_" + matnames[mat_id];
        std::vector<double> cellvec(nsrccells, 0.0);
        std::vector<double> cellvec_wtd(nsrccells, 0.0);

        std::vector<int> matcells;
        sourceStateWrapper.mat_get_cells(mat_id, &matcells);
        int nmatcells = matcells.size();

        double *matvec;
        sourceStateWrapper.mat_get_celldata("cellmatdata", mat_id, &matvec);

        double *matvolfracs;
        sourceStateWrapper.mat_get_celldata("mat_volfracs", mat_id, &matvolfracs);

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

      for (int m = 0; m < nlocal_mats; m++) {
        int mat_id = local_material_IDs[m];

        std::string varname1 = "cellmatdata_" + matnames[mat_id];
        std::string varname2 = "cellmatdata_wtd_" + matnames[mat_id];
        std::vector<double> cellvec(ntarcells, 0.0);
        std::vector<double> cellvec_wtd(ntarcells, 0.0);

        std::vector<int> matcells;
        targetStateWrapper.mat_get_cells(mat_id, &matcells);
        int nmatcells = matcells.size();

        double *matvec;
        targetStateWrapper.mat_get_celldata("cellmatdata", mat_id, &matvec);

        double *matvolfracs;
        targetStateWrapper.mat_get_celldata("mat_volfracs", mat_id, &matvolfracs);

        for (int ic = 0; ic < nmatcells; ic++) {
          int state_cell_id = matcells[ic];
          int mesh_cell_id = owned2all[state_cell_id];
          cellvec[mesh_cell_id] = matvec[ic];
          cellvec_wtd[mesh_cell_id] = matvec[ic]*matvolfracs[ic];
        }
        targetStateWrapper.mesh_add_data(Portage::CELL, varname1, &(cellvec[0]));
        targetStateWrapper.mesh_add_data(Portage::CELL, varname2, &(cellvec_wtd[0]));
      }

      targetState->export_to_mesh();
      targetMesh->write_to_exodus_file("output.exo");
    }
    
    if (rank == 0)
      std::cout << "...done." << std::endl;
  }

  // Perform interface reconstruction on target mesh for pretty pictures
  // and error computation of material fields

  // Because we use the one and only MOF, we do NOT need the material data
  // on the ghost cells

  interface_reconstructor_factory<dim, Wonton::Jali_Mesh_Wrapper>
    target_IRFactory(targetMeshWrapper, ims_tols, target_convex_cells);
  auto target_interface_reconstructor = target_IRFactory();

  target_interface_reconstructor->set_volume_fractions(target_cell_num_mats,
                                                       target_cell_mat_ids,
                                                       target_cell_mat_volfracs,
                                                       target_cell_mat_centroids);
  target_interface_reconstructor->reconstruct(executor);

  if (mesh_output) {  
    std::string filename = "target_mm_" + std::to_string(rank) + ".gmv";
    Portage::write_to_gmv<dim>(targetMeshWrapper, targetStateWrapper,
                               target_interface_reconstructor, fieldnames,
                               filename);
  }

  // Compute error
  L1_error = 0.0; L2_error = 0.0;
  if (material_field_expressions.size()) {
    double error, toterr = 0.0;
    double const * cellmatvals;
    double totvolume = 0.;
    for (int m = 0; m < nlocal_mats; m++) {
      int mat_id = local_material_IDs[m];

      targetStateWrapper.mat_get_celldata<double>("cellmatdata", mat_id, &cellmatvals);

      std::vector<int> matcells;
      targetStateWrapper.mat_get_cells(mat_id, &matcells);

      // Cell error computation
      int nmatcells = matcells.size();
      for (int ic = 0; ic < nmatcells; ++ic) {
        int state_cell_id = matcells[ic];
        int mesh_cell_id = owned2all[state_cell_id];
        if (target_cell_num_mats[mesh_cell_id] == 1) {
          if (target_cell_mat_ids[offsets[mesh_cell_id]] == mat_id) {
            Wonton::Point<dim> ccen;
            targetMeshWrapper.cell_centroid(mesh_cell_id, &ccen);
            error = mat_fields[mat_id](ccen) - cellmatvals[ic];

            double cellvol = targetMeshWrapper.cell_volume(mesh_cell_id);
            totvolume += cellvol;
            L1_error += fabs(error)*cellvol;
            L2_error += error*error*cellvol;
          }
        } else {
          Tangram::CellMatPoly<dim> const& cellmatpoly =
              target_interface_reconstructor->cell_matpoly_data(mesh_cell_id);
          int nmp = cellmatpoly.num_matpolys();
          for (int i = 0; i < nmp; i++) {
            if (cellmatpoly.matpoly_matid(i) == mat_id) {
              Wonton::Point<dim> mcen = cellmatpoly.matpoly_centroid(i);
              error = mat_fields[mat_id](mcen) - cellmatvals[ic];

              double matpolyvol = cellmatpoly.matpoly_volume(i);
              totvolume += matpolyvol;
              L1_error += fabs(error)*matpolyvol;
              L2_error += error*error*matpolyvol;
            }
          }
        }
      }
    }
  }

  L2_error = sqrt(L2_error);
  if (numpe > 1) {
    std::cout << std::flush << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    double globalerr;
    MPI_Reduce(&L1_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    L1_error = globalerr;

    MPI_Reduce(&L2_error, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    L2_error = globalerr;
  }

  if (rank == 0) {
    std::cout << "L1 NORM OF ERROR IS " << L1_error << std::endl;
    std::cout << "L2 NORM OF ERROR IS " << L2_error << std::endl;    
  }
}
