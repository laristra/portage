/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "portage/support/portage.h"
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

#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/support/mpi_collate.h"
#include "portage/driver/mmdriver.h"
#include "portage/support/timer.h"
#include "portage/search/search_swept_face.h"
#include "portage/intersect/intersect_swept_face.h"

#include "tangram/utility/get_material_moments.h"
#include "tangram/utility/rpgtools/cuts.h"
#include "tangram/utility/rpgtools/primitives.h"
#include "tangram/utility/rpgtools/examples/matdata_rotor3d.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/reconstruct/VOF.h"

#define IR_2D MOF

#include "portage/driver/write_to_gmv.h"

#if ENABLE_TIMINGS
  #include "portage/support/timer.h"
#endif

//#define DEBUG_PRINT

// For parsing and evaluating user defined expressions in apps
#include "user_field.h"

using Wonton::Jali_Mesh_Wrapper;

namespace RGMDApp {
  enum Problem_type {
    SRCPURECELLS, TJUNCTION, BALL, ROTOR, CHECKERBOARD
  };
  enum Mesh_perturb_type {
    NO, SHIFT, PSEUDORANDOM, ROTATE
  };
  struct Recon_params_type {
    std::vector<double> ball_radii;
    std::vector<double> ball_center;
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
  - checkerboard: domain contains two alternating materials in 2D or 8 in 3D.
  There are only single-material cells in this problem.
*/

//////////////////////////////////////////////////////////////////////

int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp_rgmd_jali " <<
      "--problem=srcpurecells|tjunction|ball|rotor|checkerboard \n" <<  
      "--dim=2|3 --nsourcecells=N --ntargetcells=M \n" << 
      "--source_convex_cells=y|n --target_convex_cells=y|n \n" <<
      "--remap_order=1|2 \n" <<
      "--limiter=barth_jespersen --bnd_limiter=zero_gradient \n"
      "--mesh_min=0. --mesh_max=1. --perturb_source=n|shift|pseudorandom|rotate \n" <<
      "--output_meshes=y|n --convergence_study=NREF --only_threads=y|n " <<
      "--field_filename=string --intersect=y|n \n\n";

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

  std::cout << "--perturb_source (default = n): " <<
      "add a shift, pseudorandom perturbation or rotation to coordinates of a source mesh\n\n";

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

  std::cout << "--intersect (default = y)\n";
  std::cout << "  If 'y', intersection-based remap is used. Set to 'n' to enable  " <<
      "swept-face remap \n\n";
  
#if ENABLE_TIMINGS
  std::cout << "--only_threads (default = n)\n";
  std::cout << " enable if you want to profile only threads scaling\n\n";

  std::cout << "--scaling (default = strong)\n";
  std::cout << " specify the scaling study type [strong|weak]\n\n";
#endif

  std::cout << "--field_filename\n\n";
  std::cout << "  If defined, the field output filename. Rank is appended if parallel\n\n";

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

  get_material_moments<Mesh_Wrapper, Wonton::CartesianCoordinates>(
      mesh, material_interfaces, mesh_material_IDs, 
      cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
      vol_tol, dst_tol, decompose_cells);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  assert(unsigned(ncells) == cell_num_mats.size());
  if (unsigned(ncells) != cell_mat_ids.size()) {
    int nextras = static_cast<int>(cell_mat_ids.size() - ncells);
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
void checkerboard_material_data(const Mesh_Wrapper& mesh,
                                std::vector<int>& mesh_material_IDs,
                                std::vector<int>& cell_num_mats,
                                std::vector<int>& cell_mat_ids,
                                std::vector<double>& cell_mat_volfracs,
                                std::vector< Wonton::Point<2> >& cell_mat_centroids,
                                const double vol_tol,
                                const double dst_tol,
                                bool decompose_cells) {
  mesh_material_IDs = {0, 1};
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  cell_num_mats.resize(ncells);
  cell_mat_ids.resize(ncells);
  cell_mat_volfracs.resize(ncells);
  cell_mat_centroids.resize(ncells);
  std::fill(cell_num_mats.begin(), cell_num_mats.end(), 1);
  for (int icell = 0; icell < ncells; icell++)
    cell_mat_ids[icell] = icell%2;
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

  get_material_moments<Mesh_Wrapper, Wonton::CartesianCoordinates>(
      mesh, material_interfaces, mesh_material_IDs, 
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
                        RGMDApp::Recon_params_type recon_params,
                        bool decompose_cells) {
  mesh_material_IDs = {1, 0};
  Wonton::Point<2> circle_cen(0.75, 0.75);
  if(!recon_params.ball_center.empty()) {
    circle_cen[0] = recon_params.ball_center[0];
    circle_cen[1] = recon_params.ball_center[1];
  }
  double circle_rad = 0.25;
  if(!recon_params.ball_radii.empty())
    circle_rad = recon_params.ball_radii[0];
  int nquadrant_samples = 90;

  get_material_moments<Mesh_Wrapper, Wonton::CartesianCoordinates>(
      mesh, mesh_material_IDs, 
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

  get_material_moments<Mesh_Wrapper, Wonton::CartesianCoordinates>(
      mesh, material_interfaces, mesh_material_IDs, 
      cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
      vol_tol, dst_tol, decompose_cells);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  assert(unsigned(ncells) == cell_num_mats.size());
  if (unsigned(ncells) != cell_mat_ids.size()) {
    int nextras = static_cast<int>(cell_mat_ids.size() - ncells);
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
void checkerboard_material_data(const Mesh_Wrapper& mesh,
                                std::vector<int>& mesh_material_IDs,
                                std::vector<int>& cell_num_mats,
                                std::vector<int>& cell_mat_ids,
                                std::vector<double>& cell_mat_volfracs,
                                std::vector< Wonton::Point<3> >& cell_mat_centroids,
                                const double vol_tol,
                                const double dst_tol,
                                bool decompose_cells) {
  mesh_material_IDs = {0, 1, 2, 3, 4, 5, 6, 7};
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  cell_num_mats.resize(ncells);
  cell_mat_ids.resize(ncells);
  cell_mat_volfracs.resize(ncells);
  cell_mat_centroids.resize(ncells);
  std::fill(cell_num_mats.begin(), cell_num_mats.end(), 1);
  for (int icell = 0; icell < ncells; icell++) {
    int nc = int(std::cbrt(ncells));
    int i = icell/(nc*nc);
    int j = (icell%(nc*nc))/nc;
    int k = icell%nc;
    cell_mat_ids[icell] = i%2+2*(j%2)+4*(k%2);
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

  get_material_moments<Mesh_Wrapper, Wonton::CartesianCoordinates>(
      mesh, material_interfaces, mesh_material_IDs, 
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
                        RGMDApp::Recon_params_type recon_params,
                        bool decompose_cells) {
  mesh_material_IDs = {2, 0, 1};
  Wonton::Point<3> spheres_cen(0.75, 0.75, 0.75);
  if(!recon_params.ball_center.empty()) {
    spheres_cen[0] = recon_params.ball_center[0];
    spheres_cen[1] = recon_params.ball_center[1];
    spheres_cen[2] = recon_params.ball_center[2];
  }
  std::vector<double> sphere_rads = {0.2495, 0.25};
  if(!recon_params.ball_radii.empty()) {
    sphere_rads[0] = recon_params.ball_radii[0];
    sphere_rads[1] = recon_params.ball_radii[1];
  }
  int nquadrant_samples = 24;

  get_material_moments<Mesh_Wrapper, Wonton::CartesianCoordinates>(
      mesh, mesh_material_IDs, 
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
                           int interp_order, bool intersect,
                           double srclo, double srchi,
                           RGMDApp::Mesh_perturb_type mesh_perturb,
                           RGMDApp::Problem_type problem,
                           RGMDApp::Recon_params_type recon_params,
                           std::vector<std::string> material_field_expressions,
                           std::string field_filename,
                           bool mesh_output, int rank, int numpe,
                           Jali::Entity_kind entityKind,
                           double& L1_error, double& L2_error,
			   int iteration,
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
  RGMDApp::Mesh_perturb_type mesh_perturb = RGMDApp::Mesh_perturb_type::NO;
  int nsourcecells = 0, ntargetcells = 0;  // No default
  int dim = 2;
  bool source_convex_cells = true, target_convex_cells = true;
  std::vector<std::string> material_field_expressions;
  std::string srcfile, trgfile;  // No default
  std::string field_filename;  // No default;

  int interp_order = 1;
  bool intersect = true;
  // since write_to_gmv segfaults in parallel, default to false and force the
  // user to output in serial
  bool mesh_output = false;
  int n_converge = 1;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  Portage::Limiter_type limiter = Portage::Limiter_type::NOLIMITER;
  Portage::Boundary_Limiter_type bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;
  double srclo = 0.0, srchi = 1.0;  // bounds of generated mesh in each dir
  RGMDApp::Recon_params_type recon_params;

#if ENABLE_TIMINGS
  bool only_threads = false;
  std::string scaling_type = "strong";
#endif

  // Parse the input
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of('=');
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
      else if (valueword == "checkerboard")
        problem = RGMDApp::Problem_type::CHECKERBOARD;
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
    } else if (keyword == "material_fields" or keyword == "ball_radii" or
               keyword == "ball_center") {
      // Expecting comma-separated quoted expressions
      std::string const& exprlist = valueword;
      std::size_t expr_beg = 0;
      std::size_t expr_end = 0;
      bool exprdone = false;
      std::string expr;
      while (!exprdone) {
        expr_end = exprlist.find_first_of(',', expr_beg);
        if (expr_end == std::string::npos) {
          exprdone = true;
          expr = exprlist.substr(expr_beg);
        } else {
          expr = exprlist.substr(expr_beg, expr_end-expr_beg);
          expr_beg = expr_end+1;
        }
        if (keyword == "material_fields")
          material_field_expressions.push_back(expr);
        else if (keyword == "ball_radii")
          recon_params.ball_radii.push_back(std::stod(expr));
        else if (keyword == "ball_center")
          recon_params.ball_center.push_back(std::stod(expr));
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
      srclo = stod(valueword);
    } else if (keyword == "mesh_max") {
      srchi = stod(valueword);
    } else if (keyword == "perturb_source") {
      if(valueword == "n")
        mesh_perturb = RGMDApp::Mesh_perturb_type::NO;
      else if(valueword == "shift")
        mesh_perturb = RGMDApp::Mesh_perturb_type::SHIFT;
      else if (valueword == "pseudorandom")
        mesh_perturb = RGMDApp::Mesh_perturb_type::PSEUDORANDOM;
      else if (valueword == "rotate")
        mesh_perturb = RGMDApp::Mesh_perturb_type::ROTATE;
      else {
        std::cerr << "Unknown mesh perturbation type!\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
    } else if (keyword == "intersect") {
      intersect = (valueword == "y");
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
    } else if (keyword == "field_filename") {
      field_filename=valueword;
    } else {
      std::cerr << "Unrecognized option " << keyword << " !\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Some input error checking

  if (nsourcecells > 0 && srcfile.length() > 0) {
    std::cerr << "Cannot request internally generated source mesh "
              << "(--nsourcecells) and external file read (--source_file)\n\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if (!nsourcecells && srcfile.length() == 0) {
    std::cerr << "Must specify one of the two options --nsourcecells "
              << "or --source_file\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (ntargetcells > 0 && trgfile.length() > 0) {
    std::cerr << "Cannot request internally generated target mesh "
              << "(--ntargetcells) and external file read (--target_file)\n\n";
    print_usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if (!ntargetcells && trgfile.length() == 0) {
    std::cerr << "Must specify one of the two options --ntargetcells "
              << "or --target_file\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if ((srcfile.length() || trgfile.length()) && n_converge > 1) {
    std::cerr <<
        "Convergence study possible only for internally generated meshes\n";
    std::cerr << "Will do single remap and exit\n";
    n_converge = 1;
  }
  if (nsourcecells > 0)
    if (material_field_expressions.empty()) {
      std::cerr << "No field imposed on internally generated source mesh\n";
      std::cerr << "Nothing to remap. Exiting...";
      print_usage();
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  if (not intersect) {
    if (nsourcecells != ntargetcells) {
      std::cerr << "Need the same mesh topology for a swept-face remap\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (source_convex_cells)
      std::cerr << "source_convex_cells should be false with a swept remap\n";
  }
  if (dim == 2) {
    if (!recon_params.ball_center.empty() && recon_params.ball_center.size() != 2) {
      std::cerr << "ball_center should contain 2 components separated be a comma\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (!recon_params.ball_radii.empty() && recon_params.ball_radii.size() != 1) {
      std::cerr << "ball_radii should contain one number in 2D\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  if (dim == 3) {
    if (!recon_params.ball_center.empty() && recon_params.ball_center.size() != 3) {
      std::cerr << "ball_center should contain 3 components separated be commas\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (!recon_params.ball_radii.empty() && recon_params.ball_radii.size() != 2) {
      std::cerr << "ball_radii should contain 2 radii for two concentric spheres\n";
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
  profiler->params.nmats   = material_field_expressions.size();
  switch(problem) {
    case RGMDApp::Problem_type::SRCPURECELLS: profiler->params.output  = "srcpurecells_"; break;
    case RGMDApp::Problem_type::TJUNCTION: profiler->params.output  = "t-junction_"; break;
    case RGMDApp::Problem_type::BALL: profiler->params.output  = "ball_"; break;
    case RGMDApp::Problem_type::ROTOR: profiler->params.output  = "rotor_"; break;
    case RGMDApp::Problem_type::CHECKERBOARD: profiler->params.output  = "checkerboard_"; break;
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
  mf.included_entities(Jali::Entity_kind::ALL_KIND);

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

    // Now run the remap on the meshes and get back the L2 error
    switch (dim) {
      case 2:
        run<2>(source_mesh, target_mesh, source_convex_cells, target_convex_cells,
               limiter, bnd_limiter, interp_order, intersect,
               srclo, srchi, mesh_perturb,
               problem, recon_params, material_field_expressions,
               field_filename, mesh_output,
               rank, numpe, entityKind, l1_err[i], l2_err[i], i, profiler);
        break;
      case 3:
        run<3>(source_mesh, target_mesh, source_convex_cells, target_convex_cells,
               limiter, bnd_limiter, interp_order, intersect,
               srclo, srchi, mesh_perturb,
               problem, recon_params, material_field_expressions,
               field_filename, mesh_output,
               rank, numpe, entityKind, l1_err[i], l2_err[i], i, profiler);
        break;
      default:
        std::cerr << "Dimension not 2 or 3" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
#ifndef NDEBUG
    std::cout << "L1 norm of error for iteration " << i << " is " <<
        l1_err[i] << std::endl;
    std::cout << "L2 norm of error for iteration " << i << " is " <<
        l2_err[i] << std::endl;
#endif
    // if convergence study, double the mesh resolution
    nsourcecells *= 2;
    ntargetcells *= 2;
  }
#ifndef NDEBUG 
  for (int i = 1; i < n_converge; i++) {
    std::cout << "Error ratio L1(" << i - 1 << ")/L1(" << i << ") is " <<
        l1_err[i - 1]/l1_err[i] << std::endl;
    std::cout << "Error ratio L2(" << i - 1 << ")/L2(" << i << ") is " <<
        l2_err[i - 1]/l2_err[i] << std::endl;
  }
#endif
  MPI_Finalize();

#if ENABLE_TIMINGS
  profiler->time.total = timer::elapsed(start);

  // dump timing data
  if (rank == 0) {
    profiler->dump();
  }
#endif
}


// write a field a file in the format needed by distributed_cmp
// only works for scalar fields at the moment
void write_field(std::string filename, 
    const Wonton::Jali_Mesh_Wrapper& meshWrapper, 
    const Wonton::Jali_State_Wrapper& stateWrapper,
    std::string field_name="cellmatdata"){
  
  // open the stream
  std::ofstream f(filename);
  
  // get the number of materials in the problem frm the state manager
  int nmats=stateWrapper.num_materials();
  
  // loop over materials since we are material dominant
  for (int m = 0; m < nmats;  ++m) {
  
    // get the material cells
    std::vector<int> matcells;
    stateWrapper.mat_get_cells(m, &matcells);

    // get the field values
    double const *data;
    stateWrapper.mat_get_celldata(field_name, m, &data);
    
    // write a line for each cell (using global id) for each material
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      f << meshWrapper.get_global_id(matcells[ic], Wonton::Entity_kind::CELL)
      << " " << m << " " 
      << std::fixed
      << std::setprecision(16) 
      << data[ic]<< "\n";
    }    
  }
  
  // close the stream
  f.close();
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
                   int interp_order, bool intersect,
                   double srclo, double srchi,
                   RGMDApp::Mesh_perturb_type mesh_perturb,
                   RGMDApp::Problem_type problem,
                   RGMDApp::Recon_params_type recon_params,
                   std::vector<std::string> material_field_expressions,
                   std::string field_filename, bool mesh_output,
                   int rank, int numpe, Jali::Entity_kind entityKind,
                   double& L1_error, double& L2_error,
	           int iteration,
                   std::shared_ptr<Profiler> profiler) {
  if (rank == 0)
    std::cout << "starting portageapp_rgmd_jali...\n";

  // Add perturbations to a source mesh
  if(mesh_perturb != RGMDApp::Mesh_perturb_type::NO) {
    int const numsrcnodes = sourceMesh->num_nodes<Jali::Entity_type::ALL>();
    int const numsrccells = sourceMesh->num_cells<Jali::Entity_type::ALL>();
    double dx = (srchi-srclo)/pow(numsrccells,1.0/dim);
    std::array<double, dim> pnt;
    std::array<double, dim> new_pnt;
    for (int i = 0; i < numsrcnodes; i++) {
      sourceMesh->node_get_coordinates(i, &pnt);
      if (mesh_perturb == RGMDApp::Mesh_perturb_type::PSEUDORANDOM) {
        for (int dir = 0; dir < 2; dir++)
          new_pnt[dir] = pnt[dir] + dx/20.0 * sin( 2.0*M_PI / (srchi-srclo) *
                  (pnt[dir]-srclo)*(pnt[(dir+1)%2]-srclo) *
                  (srchi-pnt[dir])*(srchi-pnt[(dir+1)%2]) * 9001 );
        if (dim==3) new_pnt[2] = pnt[2];
      }
      else if (mesh_perturb == RGMDApp::Mesh_perturb_type::SHIFT) { // perturbation ala Misha
        double alpha = 0.1*dx;
        new_pnt[0] = (1. - alpha) * pnt[0] + alpha * pnt[0]*pnt[0]*pnt[0];
        new_pnt[1] = (1. - alpha) * pnt[1] + alpha * pnt[1]*pnt[1];
        if (dim==3) new_pnt[2] = pnt[2];
      }
      else if(mesh_perturb == RGMDApp::Mesh_perturb_type::ROTATE) {
        if (dim==3) { // scaling by 4 and rotations by pi/4 around z and by arccos(1/sqrt(3)) around y
          new_pnt[0] = 4.0 * (  1.0 / sqrt(6.0)* pnt[0] + 1.0 / sqrt(6.0)* pnt[1] - sqrt(2.0/3.0)* pnt[2] );
          new_pnt[1] = 4.0 * ( -1.0 / sqrt(2.0)* pnt[0] + 1.0 / sqrt(2.0)* pnt[1] );
          new_pnt[2] = 4.0 * (  1.0 / sqrt(3.0)* pnt[0] + 1.0 / sqrt(3.0)* pnt[1] + 1.0 / sqrt(3.0)* pnt[2] );
        }
        else {
          new_pnt[0] = 2.0 * sqrt(2.0) * ( pnt[0] - pnt[1] );
          new_pnt[1] = 2.0 * sqrt(2.0) * ( pnt[0] + pnt[1] );
        }
      }
      sourceMesh->node_set_coordinates(i, new_pnt.data());
    }
    if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);
  }

  {//print out some info
     std::cout<<"\n\nITERATION "<<iteration<<"---->"<<std::endl;
     std::cout<<"Using intersect-based remap ? "<<intersect<<std::endl;
     std::cout<<"source_convex_cells ? "<<source_convex_cells<<std::endl;
     std::cout<<"Problem Type :"<<problem<<std::endl;
     std::cout<<"AFTER PERTURBATION "<<std::endl;
     int const numsrccells = sourceMesh->num_cells<Jali::Entity_type::ALL>();
       
     int nedges = 0; 
     double hmax = -std::numeric_limits<double>::max();
     double hmin = std::numeric_limits<double>::max();
     double havg = 0.0 ; 

     for (auto c = 0; c < numsrccells; c++){ 
      std::vector<int> edges, dirs, nodes;
      sourceMesh->cell_get_faces(c, &edges);
      int const nb_edges = edges.size();
      for (int i = 0; i < nb_edges; ++i) {
        sourceMesh->face_get_nodes(edges[i], &nodes);

        std::array<double,dim> coord1, coord2;
        sourceMesh->node_get_coordinates(nodes[0], &coord1);
        sourceMesh->node_get_coordinates(nodes[1], &coord2);
    
        double dist = (dim == 2 ? 
		       (pow(coord1[0]-coord2[0],2) + 
    		        pow(coord1[1]-coord2[1],2)) : 
		       (pow(coord1[0]-coord2[0],2) + 
		        pow(coord1[1]-coord2[1],2) +
		        pow(coord1[2]-coord2[2],2)));

        dist = sqrt(dist);
        if (dist > hmax) hmax=dist;  
        if (dist < hmin) hmin=dist;  
        havg += dist;
        nedges ++;        
       } //nedges
     }//ncells

     havg = havg/nedges; 

     double dh = (srchi-srclo)/pow(numsrccells,1.0/dim);
     std::cout<<"NCELLS = "<<numsrccells<<std::endl;
     std::cout<<"SRC: HMIN = "<<hmin<<", HMAX = "<<hmax<<", HAVG = "<< havg<<std::endl;
     std::cout<<"TRG: DH = "<<dh<<"\n\n"<<std::endl;
  }//print info

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
    case RGMDApp::Problem_type::CHECKERBOARD:
      checkerboard_material_data<Wonton::Jali_Mesh_Wrapper>(sourceMeshWrapper,
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
                                                    recon_params,
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
  int const num_unique_ids = std::set<int>(global_material_IDs.begin(), global_material_IDs.end()).size();
  if (num_unique_ids != nglobal_mats) {
    std::cerr << "Generated materials had repeated indices, " << 
                 "check the implementation of the problem!\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int const min_id = *std::min_element(global_material_IDs.begin(),
                                       global_material_IDs.end());
  int const max_id = *std::max_element(global_material_IDs.begin(),
                                       global_material_IDs.end());

  if ((min_id != 0) || (max_id != nglobal_mats - 1) ) {
    std::cerr << "Generated materials didn't have contiguous indices, " << 
                 "check if there are material subdomains not covered by the mesh!\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int const num_field_expressions = material_field_expressions.size();

  if (num_field_expressions != nglobal_mats) {
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
#ifndef NDEBUG
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
#endif

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

  auto source_interface_reconstructor =
      std::make_shared<Tangram::Driver<Tangram::MOF,
                                       dim,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Tangram::SplitRnD<dim>,
                                       Tangram::ClipRnD<dim>
                                       >
                       >(sourceMeshWrapper, ims_tols, source_convex_cells);


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
          std::vector<double> cell_mat_moments = cellmatpoly.material_moments(m);
          Wonton::Point<dim> mcen;
          for (int dir = 1; dir <= dim; dir++)
            mcen[dir-1] = cell_mat_moments[dir] / cell_mat_moments[0];
          matData[ic] = mat_fields[m](mcen);
        }
      }

      sourceStateWrapper.mat_add_celldata("cellmatdata", m, &(matData[0]));
    }
  }

#ifndef NDEBUG
  if (rank == 0) {
    std::cout << "***Registered fieldnames:\n";
    for (auto & field_name: sourceStateWrapper.names()) 
      std::cout << " registered fieldname: " << field_name << std::endl;
  }
#endif
  std::vector<std::string> fieldnames;
  fieldnames.emplace_back("cellmatdata");
  
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
  if (!material_field_expressions.empty()) {
    targetStateWrapper.mat_add_celldata<double>("cellmatdata");
    remap_fields.emplace_back("cellmatdata");
  }

  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);


  if (mesh_output) {  
    std::string filename = "source_mm_" + std::to_string(rank) + 
                           "_iteration_" + std::to_string(iteration) + ".gmv";
    Portage::write_to_gmv<dim>(sourceMeshWrapper, sourceStateWrapper,
                               source_interface_reconstructor, fieldnames,
                               filename);
  }

  if(intersect) {
    if (interp_order == 1) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectRnD,
        Portage::Interpolate_1stOrder,
        dim,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::MOF,
        Tangram::SplitRnD<dim>,
        Tangram::ClipRnD<dim>>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_reconstructor_options(source_convex_cells, ims_tols);
      driver.run(executor);
    } else if (interp_order == 2) {
      Portage::MMDriver<
        Portage::SearchKDTree,
        Portage::IntersectRnD,
        Portage::Interpolate_2ndOrder,
        dim,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::MOF,
        Tangram::SplitRnD<dim>,
        Tangram::ClipRnD<dim>>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_limiter(limiter);
      driver.set_bnd_limiter(bnd_limiter);
      driver.set_reconstructor_options(source_convex_cells, ims_tols);
      driver.run(executor);
    }
  } else { // swept face
    if (interp_order == 1) {
      Portage::MMDriver<
        Portage::SearchSweptFace,
        Portage::IntersectSweptFace,
        Portage::Interpolate_1stOrder,
        dim,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::MOF,
        Tangram::SplitRnD<dim>,
        Tangram::ClipRnD<dim>>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_reconstructor_options(source_convex_cells, ims_tols);
      driver.run(executor);
    } else if (interp_order == 2) {
      Portage::MMDriver<
        Portage::SearchSweptFace,
        Portage::IntersectSweptFace,
        Portage::Interpolate_2ndOrder,
        dim,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Wonton::Jali_Mesh_Wrapper,
        Wonton::Jali_State_Wrapper,
        Tangram::MOF,
        Tangram::SplitRnD<dim>,
        Tangram::ClipRnD<dim>>
          driver(sourceMeshWrapper, sourceStateWrapper,
                 targetMeshWrapper, targetStateWrapper);
      driver.set_remap_var_names(remap_fields);
      driver.set_limiter(limiter);
      driver.set_bnd_limiter(bnd_limiter);
      driver.set_reconstructor_options(source_convex_cells, ims_tols);
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
    std::vector<int> current_matcells;
    targetStateWrapper.mat_get_cells(m, &current_matcells);

    double const *matvf;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf);

    Wonton::Point<dim> const *matcen;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen);

    int nmatcells = current_matcells.size();
    for (int ic = 0; ic < nmatcells; ic++) {
      int state_cell_id = current_matcells[ic];
      int mesh_cell_id = owned2all[state_cell_id];
      int& ncmats = target_cell_num_mats[mesh_cell_id];
      target_cell_mat_ids[offsets[mesh_cell_id] + ncmats] = m;
      target_cell_mat_volfracs[offsets[mesh_cell_id] + ncmats] = matvf[ic];
      for (int i = 0; i < dim; i++)
        target_cell_mat_centroids[offsets[mesh_cell_id] + ncmats][i] = matcen[ic][i];
      ncmats++;
    }
  }

#ifndef NDEBUG
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

    std::cout << "\n----target owned cell global id's on rank " << rank << " for material " << m << ": ";
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
#endif
#if ENABLE_TIMINGS
  profiler->time.remap += timer::elapsed(tic);

  if (rank == 0) {
    std::cout << "Remap Time: " << profiler->time.remap << std::endl;
  }
#endif  

  // Write out the meshes if requested
  if (mesh_output) {
    if (rank == 0)
      std::cout << "Dumping data to Exodus files..." << std::endl;

    if (!material_field_expressions.empty()) {
      // For now convert each material field into a mesh field with
      // zero values for cells that don't have the material

      for (int m = 0; m < nlocal_mats; m++) {
        int mat_id = local_material_IDs[m];

        std::string varname1 = "cellmatdata_" + matnames[mat_id];
        std::string varname2 = "cellmatdata_wtd_" + matnames[mat_id];
        std::vector<double> cellvec(nsrccells, 0.0);
        std::vector<double> cellvec_wtd(nsrccells, 0.0);

        std::vector<int> current_matcells;
        sourceStateWrapper.mat_get_cells(mat_id, &current_matcells);
        int nmatcells = current_matcells.size();

        double *matvec;
        sourceStateWrapper.mat_get_celldata("cellmatdata", mat_id, &matvec);

        double *matvolfracs;
        sourceStateWrapper.mat_get_celldata("mat_volfracs", mat_id, &matvolfracs);

        for (int ic = 0; ic < nmatcells; ic++) {
          int c = current_matcells[ic];
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

        std::vector<int> current_matcells;
        targetStateWrapper.mat_get_cells(mat_id, &current_matcells);
        int nmatcells = current_matcells.size();

        double *matvec;
        targetStateWrapper.mat_get_celldata("cellmatdata", mat_id, &matvec);

        double *matvolfracs;
        targetStateWrapper.mat_get_celldata("mat_volfracs", mat_id, &matvolfracs);

        for (int ic = 0; ic < nmatcells; ic++) {
          int state_cell_id = current_matcells[ic];
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

  auto target_interface_reconstructor =
      std::make_shared<Tangram::Driver<Tangram::MOF,
                                       dim,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Tangram::SplitRnD<dim>,
                                       Tangram::ClipRnD<dim>
                                       >
                       >(targetMeshWrapper, ims_tols, target_convex_cells);

  target_interface_reconstructor->set_volume_fractions(target_cell_num_mats,
                                                       target_cell_mat_ids,
                                                       target_cell_mat_volfracs,
                                                       target_cell_mat_centroids);
  target_interface_reconstructor->reconstruct(executor);

  if (mesh_output) {
    std::string filename_prefix;
    if (!field_filename.empty())
      filename_prefix = field_filename;
    else
      filename_prefix = "target_mm_";
    std::string filename = filename_prefix + std::to_string(rank) +
                           "_iteration_" + std::to_string(iteration) + ".gmv";
    Portage::write_to_gmv<dim>(targetMeshWrapper, targetStateWrapper,
                               target_interface_reconstructor, fieldnames,
                               filename);
  }

  // Compute error
  L1_error = 0.0; L2_error = 0.0;
  if (!material_field_expressions.empty()) {
    double error;
    double const * cellmatvals;
    double totvolume = 0.;
    for (int m = 0; m < nlocal_mats; m++) {
      int mat_id = local_material_IDs[m];

      targetStateWrapper.mat_get_celldata<double>("cellmatdata", mat_id, &cellmatvals);

      std::vector<int> current_matcells;
      targetStateWrapper.mat_get_cells(mat_id, &current_matcells);

      // Cell error computation
      int nmatcells = current_matcells.size();
      for (int ic = 0; ic < nmatcells; ++ic) {
        int state_cell_id = current_matcells[ic];
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
  
  // print data for distributed compare
  if (!field_filename.empty()) {   
  
    // name mangle the filename if more than one partition
    std::string field_filename_ = field_filename + (numpe==1 ? "" : "." + std::to_string(rank));
    
    std::cout << "*************writing to: " << field_filename_ <<"\n";
    
    // write the date file
    write_field(field_filename_, targetMeshWrapper, targetStateWrapper);
  } 

#ifndef NDEBUG
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

  ///////////////////////////////////////////////////////////////////////////// 
  //New error computation 
    
  //Step 1: Compute the material layout on the target mesh. 
  std::vector<int> trgex_global_material_IDs;
  std::vector<int> trgex_cell_num_mats;
  std::vector<int> trgex_cell_mat_ids;  // flattened 2D array
  std::vector<double> trgex_cell_mat_volfracs;  // flattened 2D array
  std::vector<Wonton::Point<dim>> trgex_cell_mat_centroids;  // flattened 2D array

  switch(problem) {
    case RGMDApp::Problem_type::SRCPURECELLS:
      srcpurecells_material_data<Wonton::Jali_Mesh_Wrapper>(targetMeshWrapper,
                                                            trgex_global_material_IDs,
                                                            trgex_cell_num_mats,
                                                            trgex_cell_mat_ids,
                                                            trgex_cell_mat_volfracs,
                                                            trgex_cell_mat_centroids,
                                                            vol_tol,
                                                            dst_tol,
                                                            !target_convex_cells);
      break;
    case RGMDApp::Problem_type::CHECKERBOARD:
      checkerboard_material_data<Wonton::Jali_Mesh_Wrapper>(targetMeshWrapper,
                                                            trgex_global_material_IDs,
                                                            trgex_cell_num_mats,
                                                            trgex_cell_mat_ids,
                                                            trgex_cell_mat_volfracs,
                                                            trgex_cell_mat_centroids,
                                                            vol_tol,
                                                            dst_tol,
                                                            !target_convex_cells);
      break;
    case RGMDApp::Problem_type::TJUNCTION:
      tjunction_material_data<Wonton::Jali_Mesh_Wrapper>(targetMeshWrapper,
                                                         trgex_global_material_IDs,
                                                         trgex_cell_num_mats,
                                                         trgex_cell_mat_ids,
                                                         trgex_cell_mat_volfracs,
                                                         trgex_cell_mat_centroids,
                                                         vol_tol,
                                                         dst_tol,
                                                         !target_convex_cells);
      break;
    case RGMDApp::Problem_type::BALL:
      ball_material_data<Wonton::Jali_Mesh_Wrapper>(targetMeshWrapper,
                                                    trgex_global_material_IDs,
                                                    trgex_cell_num_mats,
                                                    trgex_cell_mat_ids,
                                                    trgex_cell_mat_volfracs,
                                                    trgex_cell_mat_centroids,
                                                    vol_tol,
                                                    dst_tol,
                                                    recon_params,
                                                    !target_convex_cells);
      break;
    case RGMDApp::Problem_type::ROTOR:
      rotor_material_data<Wonton::Jali_Mesh_Wrapper>(targetMeshWrapper,
                                                     trgex_global_material_IDs,
                                                     trgex_cell_num_mats,
                                                     trgex_cell_mat_ids,
                                                     trgex_cell_mat_volfracs,
                                                     trgex_cell_mat_centroids,
                                                     vol_tol,
                                                     dst_tol,
                                                     !target_convex_cells);
      break;
    default:
      std::cerr << "Unknown problem type!\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
  }
 
 //Step 2: Reconstruct interface using the exact material layout on target
  auto target_interface_reconstructor_exact =
      std::make_shared<Tangram::Driver<Tangram::MOF,
                                       dim,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Tangram::SplitRnD<dim>,
                                       Tangram::ClipRnD<dim>
                                       >
                       >(targetMeshWrapper, ims_tols, target_convex_cells);

  target_interface_reconstructor_exact->set_volume_fractions(trgex_cell_num_mats,
                                                       trgex_cell_mat_ids,
                                                       trgex_cell_mat_volfracs,
                                                       trgex_cell_mat_centroids);
  target_interface_reconstructor_exact->reconstruct(executor);

  //Step 3: Compute offsets into flattened arrays based on trgex_cell_num_mats
  std::vector<int> offsets_exact(ntarcells, 0);
  for (int i = 0; i < ntarcells - 1; i++)
    offsets_exact[i + 1] = offsets_exact[i] + trgex_cell_num_mats[i];

  //Step 4: Convert data from cell-centric to material-centric form as we
  // will need it for comparision: vf, centroids, data
  std::vector<int> matcells_exact[nlocal_mats];
  std::vector<double> mat_volfracs_exact[nlocal_mats];
  std::vector<Wonton::Point<dim>> mat_centroids_exact[nlocal_mats];
  
  for (int c = 0; c < ntarcells; c++) {
    int ibeg = offsets_exact[c];
    for (int j = 0; j < trgex_cell_num_mats[c]; j++) {
      int mat_id = trgex_cell_mat_ids[ibeg + j];
      matcells_exact[mat_id].push_back(c);
      mat_volfracs_exact[mat_id].push_back(trgex_cell_mat_volfracs[ibeg + j]);
      mat_centroids_exact[mat_id].push_back(trgex_cell_mat_centroids[ibeg + j]);
    }
  }

  // User specified fields on target
  std::vector<double> mat_data_exact[nlocal_mats];
  for (int m = 0; m < nlocal_mats; m++) {
    int nmatcells_ex = matcells_exact[m].size();
    for (int ic = 0; ic < nmatcells_ex; ic++) {
      int c = matcells_exact[m][ic];
      if (trgex_cell_num_mats[c] == 1) {
        Wonton::Point<dim> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
        mat_data_exact[m].push_back(mat_fields[m](ccen));
      } else {
        // obtain matpoly's for this material
        Tangram::CellMatPoly<dim> const& cellmatpoly =
          target_interface_reconstructor_exact->cell_matpoly_data(c);
        
        std::vector<double> cell_mat_moments = cellmatpoly.material_moments(m);
        Wonton::Point<dim> mcen;
        for (int dir = 1; dir <= dim; dir++)
          mcen[dir-1] = cell_mat_moments[dir] / cell_mat_moments[0];

        mat_data_exact[m].push_back(mat_fields[m](mcen));
      }        
    }
  } 


 // Step 5: Add mesh fields to target state to store volume and mass errors
 // per material 
 for (int m = 0; m < nlocal_mats; m++) {
    std::string field_vol = "error_volume_" + std::to_string(local_material_IDs[m]);
    targetStateWrapper.mesh_add_data<double>( Portage::CELL, field_vol, 0.0);

    std::string field_mass = "error_mass_" + std::to_string(local_material_IDs[m]);
    targetStateWrapper.mesh_add_data<double>( Portage::CELL, field_mass, 0.0);
 } //material loop
 
 //Step 6: Compute error between remapped and material layout over target mesh
  double L1_error_vol[nlocal_mats], L2_error_vol[nlocal_mats];
  double L1_error_mass[nlocal_mats], L2_error_mass[nlocal_mats];
  double totalvol_exact[nlocal_mats], totalmass_exact[nlocal_mats];
  double totalvol[nlocal_mats], totalmass[nlocal_mats];

  if (!material_field_expressions.empty()) {
    double const * cellmatvals, *vf;
    Wonton::Point<dim> *cen;
    
    double *error_field_volume, *error_field_mass;

    for (int m = 0; m < nlocal_mats; m++) {
      totalvol_exact[m] = totalvol[m] = 0.0; 
      totalmass_exact[m] = totalmass[m] = 0.0; 

      //Get error fields from target state 
      std::string field_vol = "error_volume_" + std::to_string(local_material_IDs[m]);
      targetStateWrapper.mesh_get_data<double>( Portage::CELL, field_vol, &error_field_volume);

      std::string field_mass = "error_mass_" + std::to_string(local_material_IDs[m]);
      targetStateWrapper.mesh_get_data<double>( Portage::CELL, field_mass, &error_field_mass);

      //Get the remapped data from target state
      int mat_id = local_material_IDs[m];
      std::vector<int> current_matcells;
      targetStateWrapper.mat_get_cells(mat_id, &current_matcells);
      targetStateWrapper.mat_get_celldata("mat_volfracs", mat_id, &vf);
      targetStateWrapper.mat_get_celldata("mat_centroids", mat_id, &cen);
      targetStateWrapper.mat_get_celldata("cellmatdata", mat_id, &cellmatvals);

      //Obtain total volume and mass from exact/close approximation
      //to  material layout
      int nmatcells_ex = matcells_exact[m].size();
      for (int ic = 0; ic <nmatcells_ex; ++ic) {
        int mesh_cell_id = matcells_exact[m][ic];

        if (trgex_cell_num_mats[mesh_cell_id] == 1) {
          if (trgex_cell_mat_ids[offsets_exact[mesh_cell_id]] == mat_id) {
            Wonton::Point<dim> ccen;
            targetMeshWrapper.cell_centroid(mesh_cell_id, &ccen);
            double cellvol = targetMeshWrapper.cell_volume(mesh_cell_id);
            totalvol_exact[m] += cellvol;
            totalmass_exact[m] += mat_data_exact[m][ic]*cellvol; 
          }
        } else {
          Tangram::CellMatPoly<dim> const& cellmatpoly =
              target_interface_reconstructor_exact->cell_matpoly_data(mesh_cell_id);
          int nmp = cellmatpoly.num_matpolys();
          double vol_ex = 0.0; 
          for (int i = 0; i < nmp; i++) {
            if (cellmatpoly.matpoly_matid(i) == mat_id) { 
             vol_ex += cellmatpoly.matpoly_volume(i);
            }
          }
          totalvol_exact[m] += vol_ex;
          totalmass_exact[m] += mat_data_exact[m][ic]*vol_ex;       
        }
      }
 
      //Obtain total volume and mass from remapped values
      int nmatcells_cur = current_matcells.size();
      for (int ic = 0; ic < nmatcells_cur; ++ic) {
        int state_cell_id = current_matcells[ic];
        int mesh_cell_id = owned2all[state_cell_id];

        if (target_cell_num_mats[mesh_cell_id] == 1) {
          if (target_cell_mat_ids[offsets[mesh_cell_id]] == mat_id) {
            double cellvol = targetMeshWrapper.cell_volume(mesh_cell_id);
            totalvol[m] += cellvol;
            totalmass[m] += cellmatvals[ic]*cellvol; 
          }
        } else {
          Tangram::CellMatPoly<dim> const& cellmatpoly =
              target_interface_reconstructor->cell_matpoly_data(mesh_cell_id);
          int nmp = cellmatpoly.num_matpolys();
          double vol = 0.0; 
          for (int i = 0; i < nmp; i++) {
            if (cellmatpoly.matpoly_matid(i) == mat_id) { 
             vol += cellmatpoly.matpoly_volume(i);
            }
          }
          totalvol[m] += vol;
          totalmass[m] += cellmatvals[ic]*vol;    
        }
      }

      // Fill out error fields 
      std::cout<<"mat "<<m<<" nmatcells :"<<current_matcells.size()<<std::endl;
      for (int ic = 0; ic < nmatcells_cur; ++ic) {
        int state_cell_id = current_matcells[ic];
        int mesh_cell_id = owned2all[state_cell_id];

        auto cell_in_mat = std::find(matcells_exact[m].begin(), 
        matcells_exact[m].end(), mesh_cell_id);

        if (cell_in_mat == matcells_exact[m].end()) {
          //this target cell after remap doesn't have the current
          //material when compared to the exact or a close approximation
          //of the material layout on the target mesh.
       
          if (target_cell_num_mats[mesh_cell_id] == 1) {
            if (target_cell_mat_ids[offsets[mesh_cell_id]] == mat_id) {
              double cellvol = targetMeshWrapper.cell_volume(mesh_cell_id);
              error_field_volume[mesh_cell_id] = cellvol; 
              error_field_mass[mesh_cell_id] = cellvol*cellmatvals[ic]; 
            }
          } else {
             Tangram::CellMatPoly<dim> const& cellmatpoly =
              target_interface_reconstructor->cell_matpoly_data(mesh_cell_id);
             int nmp = cellmatpoly.num_matpolys();
             double matpolyvol=0.0;
             for (int i = 0; i < nmp; i++) {
              if (cellmatpoly.matpoly_matid(i) == mat_id) {
                matpolyvol += cellmatpoly.matpoly_volume(i);  
              }
             }
             error_field_volume[mesh_cell_id] = matpolyvol; 
             error_field_mass[mesh_cell_id] = cellmatvals[ic]*matpolyvol; 
          }
        } else {
          auto index = std::distance(matcells_exact[m].begin(), cell_in_mat);

          if (target_cell_num_mats[mesh_cell_id] == 1) { 
            if (target_cell_mat_ids[offsets[mesh_cell_id]] == mat_id) {
              double cellvol = targetMeshWrapper.cell_volume(mesh_cell_id);

              double vol_ex = 0.0; 
              if (trgex_cell_num_mats[*cell_in_mat]==1) {
                 if (trgex_cell_mat_ids[offsets_exact[*cell_in_mat]] == mat_id) {
                   vol_ex = targetMeshWrapper.cell_volume(*cell_in_mat);
                 }
              } else {
                Tangram::CellMatPoly<dim> const& cellmatpoly =
                 target_interface_reconstructor_exact->cell_matpoly_data(*cell_in_mat);
                int nmp = cellmatpoly.num_matpolys();
                for (int i = 0; i < nmp; i++) {
                 if (cellmatpoly.matpoly_matid(i) == mat_id) {
                  vol_ex += cellmatpoly.matpoly_volume(i);  
                 }
                }
              }

              error_field_volume[mesh_cell_id] = cellvol - vol_ex; 
              error_field_mass[mesh_cell_id] = cellmatvals[ic]*cellvol 
                                               - mat_data_exact[m][index]*vol_ex; 
            }
          } else {
             Tangram::CellMatPoly<dim> const& cellmatpoly =
              target_interface_reconstructor->cell_matpoly_data(mesh_cell_id);
             int nmp = cellmatpoly.num_matpolys();
             double matpolyvol=0.0;
             for (int i = 0; i < nmp; i++) {
              if (cellmatpoly.matpoly_matid(i) == mat_id) {
                matpolyvol += cellmatpoly.matpoly_volume(i);  
              }
             }

             double vol_ex = 0.0; 
              if (trgex_cell_num_mats[*cell_in_mat]==1) {
                 if (trgex_cell_mat_ids[offsets_exact[*cell_in_mat]] == mat_id) {
                   vol_ex = targetMeshWrapper.cell_volume(*cell_in_mat);
                 }
              } else {
                Tangram::CellMatPoly<dim> const& cellmatpoly =
                 target_interface_reconstructor_exact->cell_matpoly_data(*cell_in_mat);
                int nmp = cellmatpoly.num_matpolys();
                for (int i = 0; i < nmp; i++) {
                 if (cellmatpoly.matpoly_matid(i) == mat_id) {
                  vol_ex += cellmatpoly.matpoly_volume(i);  
                 }
                }
              }

            error_field_volume[mesh_cell_id] = matpolyvol - vol_ex; 
            error_field_mass[mesh_cell_id] = cellmatvals[ic]*matpolyvol 
                                             - mat_data_exact[m][index]*vol_ex; 
          }
        }
      } //cells in the material
    
    //Compute L1 and L2 error
    L1_error_vol[m] = L2_error_vol[m] = 0.0;
    L1_error_mass[m] = L2_error_mass[m] = 0.0;

    for (auto ic = 0; ic < ntarcells; ic++)
    {
       //vol
       L1_error_vol[m] += fabs(error_field_volume[ic]);
       L2_error_vol[m] += error_field_volume[ic]*error_field_volume[ic];

       //mass
       L1_error_mass[m] += fabs(error_field_mass[ic]);
       L2_error_mass[m] += error_field_mass[ic]*error_field_mass[ic];

     } 

     L2_error_vol[m] = sqrt(L2_error_vol[m]);
     L2_error_mass[m] = sqrt(L2_error_mass[m]);

     
     L1_error_vol[m] /= totalvol_exact[m];
     L1_error_mass[m] /= totalmass_exact[m];
     
    } //material loop
  } //non-empty list of materials

 
  if (rank == 0) {
    std::cout.precision(8);
    for (int m = 0; m < nlocal_mats; m++) {
      std::cout<<"\n\nMATERIAL ID "<<m<<std::endl; 
      std::cout<<" -----> ERRORS: L1 NORM "<<std::endl; 
      std::cout << "VOLUME  "<< L1_error_vol[m]<< std::endl;
      std::cout << "MASS  "<< L1_error_mass[m]<< std::endl;
      std::cout<<" -----> CONSERVATION "<<std::endl; 
      std::cout<<"Total Volume Exact  "<<totalvol_exact[m]<<std::endl;
      std::cout<<"Total Volume Remapped  "<<totalvol[m]<<std::endl;

      std::cout<<"Total Mass Exact  "<<totalmass_exact[m]<<std::endl;
      std::cout<<"Total Mass Remapped  "<<totalmass[m]<<std::endl;
      std::cout<<std::endl;
    }
  }

  if (mesh_output) {  
    std::string filename = "target_errors_" + std::to_string(rank) + "_iteration_" 
                           + std::to_string(iteration) + ".exo";
    targetState->export_to_mesh();
    targetMesh->write_to_exodus_file(filename);
  }

  ///////////////////////////////////////////////////////////////////////////// 

}
