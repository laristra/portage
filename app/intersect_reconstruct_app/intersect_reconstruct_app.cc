/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <stdlib.h>
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/driver/write_to_gmv.h"
#include "portage/search/search_simple.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"
#include "portage/intersect/intersect_r2d.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

extern "C" {
#include "portage/intersect/r2d.h"
}

/*!
 @file intersect_reconstruct_app.cc
 @brief Implements multi-material intersect using Tangram, XMOF2D, and Jali.
 The domain is a unit square. Source and target meshes are read from Exodus II files
 (Jali), the material data for the source mesh is read from a provided file. 
 Tangram driver is used to perform interface reconstuction on the source mesh
 and obtain a vector of CellMatPoly objects for multi-material cells. The reconstructor
 used by Tangram is XMOF2D.
 For target cells intersecting with multi-material source cells volume fractions and
 centroids of their intersections with single-material polygons in respective CellMatPoly
 objects are computed, and for both materials their overall volume fractions and centroids
 are found in every target cell.
 The obtained volume fractions and centroids are used to perform interface reconstruction
 on the target mesh using Tangram driver with XMOF2D as a reconstructor.
 The material polygons produced by Tangram are written out as GMV meshes both for the source
 and the target meshes.
 Requires Tangram, XMOF2D, and Jali libraries to be linked.
*/

const double seps = 1.0e-14;  //Size(area) tolerance;

/*!
 @brief Reads material data for a mesh from a provided binary file.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] mesh_data_fname Name of the file with material data.
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
*/
template <class Mesh_Wrapper>
void read_material_data(const Mesh_Wrapper& Mesh,
                        const std::string& mesh_data_fname,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Tangram::Point2>& cell_mat_centroids);

/*!
 @brief Intersects a material polygon in a source cell with a target cell and accumulates
 material moment data.
 @param[in] source_points Points corresonding to a material polygon in a source cell
 (coincident with a source cell for single-material cells).
 @param[in] target_points Points corresonding to a cell in the target mesh.
 @param[in] mat_id Material ID corresponding to a material polygon in a source cell.
 @param[in/out] mat_moments 2D vector containing moment data for a target cell by material
 (indexed first by mat_id and then by moment).
 */
void add_intersect_moments(const std::vector<Portage::Point<2>>& source_points,
                           const std::vector<Portage::Point<2>>& target_points,
                           const int& mat_id,
                           std::vector<std::vector<double>>& mat_moments);

/*!
 @brief Finds material data for a target mesh by intersecting its cells with material
 polygons of a source mesh
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.
 @param[in] SourceMesh Source mesh wrapper.
 @param[in] TargetMesh Target mesh wrapper.
 @param[in] cellmatpoly_list Vector of pointers to CellMatPoly objects corresponding to
 multi-material source cells
 @param[in] source_offsets Offsets into source_mat_ids vector
 @param[in] source_mat_ids Indices of materials in each source cell, a flat vector
 @param[out] cell_num_mats Number of material in each target cell
 @param[out] cell_mat_ids Indices of materials in each target cell, a flat vector
 @param[out] cell_mat_volfracs Volume fractions of materials in each target cell, a flat
 vector
 @param[out] cell_mat_centroids Centroids of materials in each target cell, a flat vector
*/
template <class Mesh_Wrapper>
void get_target_material_data(const Mesh_Wrapper& SourceMesh,
                              const Mesh_Wrapper& TargetMesh,
                              const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
                              cellmatpoly_list,
                              const std::vector<int>& source_offsets,
                              const std::vector<int>& source_mat_ids,
                              std::vector<int>& cell_num_mats,
                              std::vector<int>& cell_mat_ids,
                              std::vector<double>& cell_mat_volfracs,
                              std::vector<Tangram::Point2>& cell_mat_centroids);

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank = 0;
  int world_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &world_size);
  if ((comm_rank == 0) && (world_size > 1)) {
    std::cerr << std::endl << "This app is not designed to run in distributed mode!" <<
      std::endl << std::endl;
    exit(EXIT_FAILURE);
  }

  if (argc != 4) {
    std::cerr << std::endl << "Correct usage: intersect_reconstruct_app <material_data_file> <source_mesh_file> <target_mesh_file>" <<
      std::endl << std::endl;
    exit(EXIT_FAILURE);
  }

  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> source_mesh = mesh_factory(argv[2]),
                              target_mesh = mesh_factory(argv[3]);
  
  assert(source_mesh != nullptr);
  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh, true, false, false);
  assert(target_mesh != nullptr);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh, true, false, false);

  std::vector<int> scell_num_mats;
  std::vector<int> scell_mat_ids;
  std::vector<double> scell_mat_volfracs;
  std::vector<Tangram::Point2> scell_mat_centroids;
  
  std::cout << "Reading material data for the source mesh..." << std::endl;
  std::string in_data_fname = argv[1];
  read_material_data(source_mesh_wrapper, in_data_fname, scell_num_mats, scell_mat_ids,
                     scell_mat_volfracs, scell_mat_centroids);
  
  int nsource_cells = (int) scell_num_mats.size();
  std::vector<int> source_offsets(nsource_cells, 0);
  for (int isc = 0; isc < nsource_cells - 1; isc++)
    source_offsets[isc + 1] = source_offsets[isc] + scell_num_mats[isc];
  
  std::vector<int> ncells_with_nmats;
  for (int icell = 0; icell < nsource_cells; icell++) {
    int ncell_mats = scell_num_mats[icell];
    if (ncells_with_nmats.size() < ncell_mats)
      ncells_with_nmats.resize(ncell_mats, 0);
    ncells_with_nmats[ncell_mats - 1]++;
  }
  std::cout << "Source mesh has " << nsource_cells << " cells, of which:" << std::endl <<
    ncells_with_nmats[0] << " are SMC's";
  for (int inm = 0; inm < ncells_with_nmats.size() - 1; inm++)
    std::cout << ", " << ncells_with_nmats[inm + 1] << " are " << inm + 2 << "MC's";
  std::cout << std::endl;
  
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Wonton::Jali_Mesh_Wrapper> source_xmof_driver(source_mesh_wrapper);
  
  source_xmof_driver.set_volume_fractions(scell_num_mats, scell_mat_ids,
                                          scell_mat_volfracs, scell_mat_centroids);
  std::cout << "Performing interface reconstruction on the source mesh..." << std::endl;
  source_xmof_driver.reconstruct();

  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    scellmatpoly_list = source_xmof_driver.cell_matpoly_ptrs();
 
  in_data_fname.resize(in_data_fname.size() - 4);
  std::size_t path_end = in_data_fname.rfind('/');
  if (path_end != std::string::npos)
    in_data_fname = in_data_fname.substr(path_end + 1);
  std::string target_mesh_fname = argv[3];
  target_mesh_fname.resize(target_mesh_fname.size() - 4);
  path_end = target_mesh_fname.rfind('/');
  if (path_end != std::string::npos)
    target_mesh_fname = target_mesh_fname.substr(path_end + 1);
  
  std::string source_gmv_fname = in_data_fname + "_on_source_mesh.gmv";
  std::string target_gmv_fname = in_data_fname + "_mapped_to_" + target_mesh_fname +
    "_target_mesh.gmv";
  
  int max_gmv_mat_id = *std::max_element(scell_mat_ids.begin(), scell_mat_ids.end()) + 1;
  std::cout << "Writing source mesh material polygons to " << source_gmv_fname << std::endl;
  Tangram::write_to_gmv(source_mesh_wrapper, max_gmv_mat_id, scell_num_mats, scell_mat_ids,
                        scellmatpoly_list, source_gmv_fname);

  std::vector<int> tcell_num_mats;
  std::vector<int> tcell_mat_ids;
  std::vector<double> tcell_mat_volfracs;
  std::vector<Tangram::Point2> tcell_mat_centroids;
  
  std::cout << "Intersecting target mesh with source material polygons to get material data..." << std::endl;
  get_target_material_data(source_mesh_wrapper, target_mesh_wrapper, scellmatpoly_list,
                           source_offsets, scell_mat_ids, tcell_num_mats,
                           tcell_mat_ids, tcell_mat_volfracs, tcell_mat_centroids);
  
  ncells_with_nmats.clear();
  int ntarget_cells = (int) tcell_num_mats.size();
  for (int icell = 0; icell < ntarget_cells; icell++) {
    int ncell_mats = tcell_num_mats[icell];
    if (ncells_with_nmats.size() < ncell_mats)
      ncells_with_nmats.resize(ncell_mats, 0);
    ncells_with_nmats[ncell_mats - 1]++;
  }
  std::cout << "Target mesh has " << ntarget_cells << " cells, of which:" << std::endl <<
  ncells_with_nmats[0] << " are SMC's";
  for (int inm = 0; inm < ncells_with_nmats.size() - 1; inm++)
    std::cout << ", " << ncells_with_nmats[inm + 1] << " are " << inm + 2 << "MC's";
  std::cout << std::endl;
  
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Wonton::Jali_Mesh_Wrapper> target_xmof_driver(target_mesh_wrapper);
  target_xmof_driver.set_volume_fractions(tcell_num_mats, tcell_mat_ids,
                                          tcell_mat_volfracs, tcell_mat_centroids);
  std::cout << "Performing interface reconstruction on the target mesh..." << std::endl;
  target_xmof_driver.reconstruct();
  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    tcellmatpoly_list = target_xmof_driver.cell_matpoly_ptrs();

  std::cout << "Writing target mesh material polygons to " << target_gmv_fname << std::endl;
  Tangram::write_to_gmv(target_mesh_wrapper, max_gmv_mat_id, tcell_num_mats, tcell_mat_ids,
                        tcellmatpoly_list, target_gmv_fname);
  std::cout << "All done!" << std::endl;
  MPI_Finalize();
}

/*!
 @brief Reads material data for a mesh from a provided binary file.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] mesh_data_fname Name of the file with material data.
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
 computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
 vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
 requires computations of offsets
*/
template <class Mesh_Wrapper>
void read_material_data(const Mesh_Wrapper& Mesh,
                        const std::string& mesh_data_fname,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Tangram::Point2>& cell_mat_centroids) {
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();
  
  std::ifstream os(mesh_data_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl << "Cannot open " << mesh_data_fname <<
      " for binary input" << std::endl;
    throw std::runtime_error(os.str());
  }
  
  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  assert(data_dim == 2);
  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));
  assert(ncells = Mesh.num_owned_cells());
  cell_num_mats.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    os.read(reinterpret_cast<char *>(&cell_num_mats[icell]), sizeof(int));
    for (int im = 0; im < cell_num_mats[icell]; im++) {
      int imat;
      os.read(reinterpret_cast<char *>(&imat), sizeof(int));
      cell_mat_ids.push_back(imat);
    }
  }
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      cell_mat_volfracs.push_back(1.0);
      continue;
    }
    for (int im = 0; im < cell_num_mats[icell]; im++) {
      double vfrac;
      os.read(reinterpret_cast<char *>(&vfrac), sizeof(double));
      cell_mat_volfracs.push_back(vfrac);
    }
  }
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      Portage::Point2 cur_cell_cen;
      Mesh.cell_centroid(icell, &cur_cell_cen);
      cell_mat_centroids.push_back(Tangram::Point2(cur_cell_cen));
      continue;
    }
    for (int im = 0; im < cell_num_mats[icell]; im++) {
      double cen_x, cen_y;
      os.read(reinterpret_cast<char *>(&cen_x), sizeof(double));
      os.read(reinterpret_cast<char *>(&cen_y), sizeof(double));
      cell_mat_centroids.push_back(Tangram::Point2(cen_x, cen_y));
    }
  }
  os.close();
}

/*!
 @brief Intersects a material polygon in a source cell with a target cell and accumulates
 material moment data.
 @param[in] source_points Points corresonding to a material polygon in a source cell
 (coincident with a source cell for single-material cells).
 @param[in] target_points Points corresonding to a cell in the target mesh.
 @param[in] mat_id Material ID corresponding to a material polygon in a source cell.
 @param[in/out] mat_moments 2D vector containing moment data for a target cell by material
 (indexed first by mat_id and then by moment).
*/
void add_intersect_moments(const std::vector<Portage::Point<2>>& source_points,
                           const std::vector<Portage::Point<2>>& target_points,
                           const int& mat_id,
                           std::vector<std::vector<double>>& mat_moments) {
  // Intersect source candidate matpoly with target cell
  std::vector<double> moments =
    Portage::intersect_polys_r2d(source_points, target_points);
  // Accumulate moments (if any) from the intersection
  if (moments[0] > seps) {
    //Check if new mat_id is the max of all previously added
    if (mat_id >= mat_moments.size())
      mat_moments.resize(mat_id + 1);
    //Check if this is the first time this material is added
    if (mat_moments[mat_id].empty())
      mat_moments[mat_id].resize(3, 0.0);
    for (int im = 0; im < 3; im++)
      mat_moments[mat_id][im] += moments[im];
  }
}

template <class Mesh_Wrapper>
void get_target_material_data(const Mesh_Wrapper& SourceMesh,
                              const Mesh_Wrapper& TargetMesh,
                              const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
                              cellmatpoly_list,
                              const std::vector<int>& source_offsets,
                              const std::vector<int>& source_mat_ids,
                              std::vector<int>& cell_num_mats,
                              std::vector<int>& cell_mat_ids,
                              std::vector<double>& cell_mat_volfracs,
                              std::vector<Tangram::Point2>& cell_mat_centroids) {
  int ncells = TargetMesh.num_owned_cells();
  cell_num_mats = std::vector<int>(ncells, 0);
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();
  
  Portage::SearchSimple<Wonton::Jali_Mesh_Wrapper, Wonton::Jali_Mesh_Wrapper>
    search(SourceMesh, TargetMesh);
  
  // Target cells loop
  for (int icell = 0; icell < ncells; icell++) {
    //Target cell data
    std::vector<Portage::Point<2>> target_points;
    TargetMesh.cell_get_coordinates(icell, &target_points);
    double cell_size = TargetMesh.cell_volume(icell);

    // Accumulator for target cell material moments
    std::vector<std::vector<double>> mat_moments;
    // Search for candidate source cells to intersect with the target cell
    std::vector<int> sc_candidates;
    search(icell, &sc_candidates);

    //Accumulate moments of intersections of the target cell and material 
    //poly's in canditate source cells
    for (int isc : sc_candidates) {
      // Source candidate cellmatpoly
      auto cellmatpoly_ptr = cellmatpoly_list[isc];

      if (cellmatpoly_ptr != nullptr) { // Mixed material cell
        // Access all matpolys from source candidate cellmatpoly
        for (int ipoly = 0; ipoly < cellmatpoly_ptr->num_matpolys(); ipoly++) {
          // Get the matpoly material
          int mat_id = cellmatpoly_ptr->matpoly_matid(ipoly);

          // Get Tangram points for matpoly, convert to Portage points
          std::vector<Portage::Point<2>> source_points;
          for (auto p : cellmatpoly_ptr->matpoly_points(ipoly))
            source_points.push_back(Portage::Point<2>(p));
          
          add_intersect_moments(source_points, target_points, mat_id, mat_moments);
        }
      } else { // Single material cell
        std::vector<Portage::Point<2>> source_points;
        SourceMesh.cell_get_coordinates(isc, &source_points);
        int mat_id = source_mat_ids[source_offsets[isc]];

        add_intersect_moments(source_points, target_points, mat_id, mat_moments);
      }
    } // Finish gathering material moments for target cell
    
    // Populate material data for target cell
    for (int imat = 0; imat < mat_moments.size(); imat++)
      if (!mat_moments[imat].empty()) {
        cell_num_mats[icell]++;
        cell_mat_ids.push_back(imat);

        cell_mat_volfracs.push_back(mat_moments[imat][0]/cell_size);

        cell_mat_centroids.push_back(Portage::Point2(
          mat_moments[imat][1]/mat_moments[imat][0],
          mat_moments[imat][2]/mat_moments[imat][0]));
      }
  }
}
