/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SRC_DRIVER_DRIVER_H_
#define SRC_DRIVER_DRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersectClipper.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"
#include "portage/distributed/mpi_bounding_boxes.h"

/*!
  @file driver.h
  @brief Example driver for mapping between two Jali meshes.

  This should serve as a good example for how to write your own driver routine
  and datastructures.
*/

namespace Portage {

  /*!
    @class MeshWrapperDual "driver.h"
    @brief Wrapper for dual mesh.

    Utilizes a Jali_Mesh_Wrapper to the original mesh, but treats
    the nodes of the original mesh as the centroids of the dual mesh.
  */
template<class Mesh_Wrapper_Type>
class MeshWrapperDual {  // cellid is the dual cell (i.e. node) id
 public:
  /*!
    @brief Constructor of a wrapper to a 2d mesh.
    @param[in] w Jali_Mesh_Wrapper to original mesh.
  */
  explicit MeshWrapperDual(const Mesh_Wrapper_Type &w) : w_(w) {}

  /*!
    @brief Get the volume of the dual cell.
    @param[in] dualCellId The index of the dual cell.
    @return The volume of the specified dual cell.
  */
  double cell_volume(const int dualCellId) const {
    return w_.dual_cell_volume(dualCellId);
  }

  /*!
    @brief Get the spatial dimensions of the mesh.
    @return The spatial dimension of the mesh.
  */
  int space_dimension() const { return w_.space_dimension(); }

  /*!
    @brief Get the number of cells on this processor.
    @return The number of cells owned on this processor.
  */
  int num_owned_cells() const { return w_.num_owned_nodes(); }

  /*!
    @brief Get the number of ghost (at domain boundaries @e and processor
    boundaries) cells for this processor.
    @return The number of ghost cells on this processor.
  */
  int num_ghost_cells() const { return w_.num_ghost_nodes(); }

  /*!
    @brief Gets the coordinates of the cell centroid in the dual mesh.
    @param[in] dualcellid The dual cell id (i.e. the node id in the original
    mesh).
    @param[in,out] pplist The list of coordinate points for each node of
    the dual cell given by @c dualcellid
  */

  template<long D>
  void cell_get_coordinates(int const dualcellid,
                            std::vector<Portage::Point<D>> *pplist) const {
    w_.dual_cell_get_coordinates(dualcellid, pplist);
  }

  /*!
    @brief Get an iterator to the start of the vector of @e entity-type
    objects in the dual mesh.
    @param[in] entity Which type of data do you want to iterate over (e.g.
    @c CELL, @c NODE, etc.)
    @return An iterator pointing to the beginning of the list of @c entity data
    for this mesh.
  */
  counting_iterator begin(Entity_kind const entity) const {
    if (entity == NODE) return w_.begin(CELL);
    return w_.begin(NODE);
  }

  /*!
    @brief Get an iterator to the end of the vector of @e entity-type objects
    in the dual mesh.
    @param[in] entity Which type of data do you want to iterate over (e.g.
    @c CELL, @c NODE, etc.)
    @return An iterator pointing to the end of the list of @c entity data for
    this mesh.
  */
  counting_iterator end(Entity_kind const entity) const {
    if (entity == NODE) return w_.end(CELL);
    return w_.end(NODE);
  }

  /*!
    @brief Gets the coordinates of the cell centroid in the dual mesh.
    @param[in] dualcellid The dual cell id (i.e. the node id in the original
    mesh).
    @return The list of (x,y) coordinate pairs for each node of the
    dual cell given by @c dualcellid
    @todo Remove this in favor of @c cell_get_coordinates() ?
  */
  std::vector<Portage::Point<2>> cellToXY(int const dualcellID) const {
    std::vector<Portage::Point<2>> cellPoints;
    cell_get_coordinates(dualcellID, &cellPoints);
    return cellPoints;
  }

  /*!
    @brief Get the IDs of all cells that share a node with the specified cell
    <em> of the original mesh </em>.

    Sharing of cells is determined from the Entity_type (e.g. @c OWNED,
    @c GHOST, @c ALL ).

    @param[in] dualcellID The cell ID for which you would like to find the
    neighbors.
    @param[in] type The type of data you want (@c OWNED, @c GHOST, @c ALL)
    @param[in,out] adjcells List of IDs of adjacent cells.
  */
  void cell_get_node_adj_cells(int const dualcellID,
                               Entity_type const type,
                               std::vector<int> *adjcells) const {
    w_.node_get_cell_adj_nodes(dualcellID, type, adjcells);
  }

  /*!
    @brief Get the IDs of all cells that share a node with the specified cell
    of the <em> dual mesh </em>.

    Sharing of cells is determined from the Entity_type (e.g. @c OWNED,
    @c GHOST, @c ALL).

    @param[in] dualnodeID The cell ID for which you would like to find the
    neighbors.
    @param[in] type The type of data you want (@c OWNED, @c GHOST, @c ALL)
    @param[in,out] adjnodes List of IDs of adjacent cells.

    @todo Clarify this wrt to @c MeshWrapperDual::cell_get_node_adj_cells()
  */
  void dual_cell_get_node_adj_cells(int const dualnodeID,
                                    Entity_type const type,
                                    std::vector<int> *adjnodes) const {
    w_.cell_get_node_adj_cells(dualnodeID, type, adjnodes);
  }

  /*!
    @brief Get the coordinates of the centroid of a given dual mesh cell ID.
    @param[in] dualcellID ID of the cell in the dual mesh.
    @param[in,out] centroid (x,y,z) coordinates of the cell center (for 3d).
  */
  void cell_centroid(int const dualcellID, std::vector<double> *centroid)
      const {
    w_.dual_cell_centroid(dualcellID, centroid);
  }

  /*!
    @brief Get the coordinates of the centroid of a given dual mesh cell ID.
    @param[in] dualcellID ID of the cell in the dual mesh.
    @param[in,out] centroid (x,y,z) coordinates of the cell center (for 3d).
    @todo Clarify this wrt to @c MeshWrapperDual::cell_centroid().
  */
  void dual_cell_centroid(int const dualnodeID, std::vector<double> *centroid)
      const {
    w_.cell_centroid(dualnodeID, centroid);
  }

  typedef std::array<Portage::Point<3>, 4> wedgeCoords;
  /*!
    @brief Get the coordinates of the points that make up the wedge.

    A wedge corresponds to a cell center, a face center, and two points that
    share an edge.

    @param[in] dualcellid ID of the cell in the dual mesh.
    @param[in,out] wcoords (x,y,z) coordinates of each of the four points that
    comprise the tetrahedron that is the wedge.
  */
  void wedges_get_coordinates(int const dualcellid,
                              std::vector<wedgeCoords> *wcoords) const {
    w_.dual_wedges_get_coordinates(dualcellid, wcoords);
  }

  // Get the simplest possible decomposition of a 3D cell into tets.
  // For a dual mesh, that means returning a list of wedges.
  void decompose_cell_into_tets(int const dualcellid,
                                std::vector<wedgeCoords> *tcoords) const {
    wedges_get_coordinates(dualcellid, tcoords);
  }


 private:
  const Mesh_Wrapper_Type &w_;
};


// Forward definition
template <typename SearchType, typename IsectType, typename InterpType>
    struct RemapFunctor;

/*!
  @class Driver "driver.h"
  @brief Driver provides the API to mapping from one mesh to another.
  @tparam SourceMesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.  See Jali_Mesh_Wrapper
  for an example.
  @tparam SourceState_Wrapper A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.  See
  Jali_State_Wrapper for an example.
  @tparam TargetMesh_Wrapper A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.  See Jali_Mesh_Wrapper
  for an example.
  @tparam TargetState_Wrapper A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.  See
  Jali_State_Wrapper for an example.

*/
template <class SourceMesh_Wrapper, class SourceState_Wrapper,
    class TargetMesh_Wrapper = SourceMesh_Wrapper,
    class TargetState_Wrapper = SourceState_Wrapper>
class Driver {
 public:
 /*!
   @brief Constructor for running the interpolation driver.
   @param[in] sourceMesh A @c SourceMesh_Wrapper to the source mesh.
   @param[in] sourceState A @c SourceState_Wrapperfor the data that lives on the
   source mesh.
   @param[in] targetMesh A @c TargetMesh_Wrapper to the target mesh.
   @param[in,out] targetState A @c TargetState_Wrapper for the data that will
   be mapped to the target mesh.
 */
Driver(SourceMesh_Wrapper const & sourceMesh,
       SourceState_Wrapper const & sourceState,
       TargetMesh_Wrapper const & targetMesh,
       TargetState_Wrapper const & targetState)
 : source_mesh_(sourceMesh), source_state_(sourceState),
 target_mesh_(targetMesh), target_state_(targetState),
 interp_order_(1), dim_(sourceMesh.space_dimension()) {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
  }

/// Copy constructor (disabled)
Driver(const Driver &) = delete;

/// Assignment operator (disabled)
Driver & operator = (const Driver &) = delete;

/// Destructor
~Driver() {}

 /*!
   @brief Specify the names of the variables to be interpolated
   @param[in] remap_var_names A list of variable names of the variables to
   interpolate from the source mesh to the target mesh.  This variable must
   exist in both meshes' state manager
 */
void set_remap_var_names(
    std::vector<std::string> const &remap_var_names) {
  source_remap_var_names_ = remap_var_names;
  target_remap_var_names_ = remap_var_names;
}

 /*!
   @brief Specify the names of the variables to be interpolated
   @param[in] source_remap_var_names A list of the variables names of the
   variables to interpolate from the source mesh.
   @param[in] target_remap_var_names  A list of the variables names of the
   variables to interpolate to the target mesh.
 */
void set_remap_var_names(
    std::vector<std::string> const &source_remap_var_names,
    std::vector<std::string> const &target_remap_var_names) {
  assert(source_remap_var_names.size() == target_remap_var_names.size());

  int nvars = source_remap_var_names.size();
  for (int i = 0; i < nvars; ++i)
    assert(source_state_.get_entity(source_remap_var_names[i]) ==
           target_state_.get_entity(target_remap_var_names[i]));

  source_remap_var_names_ = source_remap_var_names;
  target_remap_var_names_ = target_remap_var_names;
}

 /*!
   @brief Get the names of the variables to be remapped from the
   source mesh.
   @return A vector of variable names to be remapped.
 */
std::vector<std::string> source_remap_var_names() const {
  return source_remap_var_names_;
}

 /*!
   @brief Get the names of the variables to be remapped to the
   target mesh.
   @return A vector of variable names to be remapped.
 */
std::vector<std::string> target_remap_var_names() const {
  return target_remap_var_names_;
}

/// Set the order of accuracy of interpolation

void set_interpolation_order(unsigned int const order) {
  interp_order_ = order;
}

 /*!
   @brief Get the order of accuracy of interpolation
   @return The order of accuracy for the interpolation.
 */
unsigned int interpolation_order() const {
  return interp_order_;
}

 /*!
   @brief Get the dimensionality of the meshes.
   @return The dimensionality of the meshes.
 */
unsigned int dim() const {
  return dim_;
}

 /*!
   @brief This method calls specialized functions to do the remapping
   based on the dimensionality of the mesh, the type of data and the
   order of interpolation.

   The individual routines run specialized search, intersect, and
   interpolation routines needed to map one mesh to another. Most of the
   heavy lifting in these routines is via a @c Portage::transform()
   over the cells in the target mesh, applying a custom @c
   RemapFunctor() (defined below) that specifies how the search,
   intersect, and interpolation calculations should be performed.
 */

void run() {
  std::printf("in Driver::run()...\n");

  int numTargetCells = target_mesh_.num_owned_cells();
  std::cout << "Number of target cells in target mesh "
  << numTargetCells << std::endl;

  int nvars = source_remap_var_names_.size();

  // Collect all cell based variables and remap them
  {
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    std::vector<std::string> source_cellvar_names;
    std::vector<std::string> target_cellvar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == CELL) {
        source_cellvar_names.emplace_back(source_remap_var_names_[i]);
        target_cellvar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_cellvar_names.size() > 0) {
      switch (dim_) {
        case 1:
          std::cerr << "Remapping not implemented for 1D" << std::endl;
          exit(-1);
        case 2: {
          (interp_order_ == 1) ?
              run_2D_CELL_order1(source_cellvar_names, target_cellvar_names) :
              run_2D_CELL_order2(source_cellvar_names, target_cellvar_names);
          break;
        }
        case 3: {
          (interp_order_ == 1) ?
              (comm_size > 1) ? 
                 run_3D_CELL_order1_distributed(source_cellvar_names, target_cellvar_names) :
                 run_3D_CELL_order1(source_cellvar_names, target_cellvar_names) :
              run_3D_CELL_order2(source_cellvar_names, target_cellvar_names);
          break;
        }
        default:
          std::cerr << "Invalid dimension" << std::endl;
          exit(-1);
      }
    }
  }


  // Collect all node based variables and remap them
  {
    std::vector<std::string> source_nodevar_names;
    std::vector<std::string> target_nodevar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == NODE) {
        source_nodevar_names.emplace_back(source_remap_var_names_[i]);
        target_nodevar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_nodevar_names.size() > 0) {
      switch (dim_) {
        case 1:
        std::cerr << "Remapping not implemented for 1D" << std::endl;
        exit(-1);
        case 2: {
          (interp_order_ == 1) ?
              run_2D_NODE_order1(source_nodevar_names, target_nodevar_names) :
              run_2D_NODE_order2(source_nodevar_names, target_nodevar_names);
          break;
        }
        case 3: {
          (interp_order_ == 1) ?
          run_3D_NODE_order1(source_nodevar_names, target_nodevar_names) :
          run_3D_NODE_order2(source_nodevar_names, target_nodevar_names);
          break;
        }
        default:
        std::cerr << "Invalid dimension" << std::endl;
        exit(-1);
      }
    }
  }
}


/// @brief 1st order remapping of cell centered data on 2D meshes
void run_2D_CELL_order1(std::vector<std::string> source_cellvar_names,
                        std::vector<std::string> target_cellvar_names);
/// @brief 2nd order remapping of cell centered data on 2D meshes
void run_2D_CELL_order2(std::vector<std::string> source_cellvar_names,
                        std::vector<std::string> target_cellvar_names);
/// @brief 1st order remapping of cell centered data on 3D meshes
void run_3D_CELL_order1(std::vector<std::string> source_cellvar_names,
                        std::vector<std::string> target_cellvar_names);
/// @brief Disributed 1st order remapping of cell centered data on 3D meshes
void run_3D_CELL_order1_distributed(std::vector<std::string> source_cellvar_names,
                                    std::vector<std::string> target_cellvar_names);
/// @brief 2nd order remapping of cell centered data on 3D meshes
void run_3D_CELL_order2(std::vector<std::string> source_cellvar_names,
                        std::vector<std::string> target_cellvar_names);
/// @brief 1st order remapping of node centered data on 2D meshes
void run_2D_NODE_order1(std::vector<std::string> source_nodevar_names,
                        std::vector<std::string> target_nodevar_names);
/// @brief 2nd order remapping of node centered data on 2D meshes
void run_2D_NODE_order2(std::vector<std::string> source_nodevar_names,
                        std::vector<std::string> target_nodevar_names);
/// @brief 1st order remapping of node centered data on 3D meshes
void run_3D_NODE_order1(std::vector<std::string> source_nodevar_names,
                        std::vector<std::string> target_nodevar_names);
/// @brief 2nd order remapping of node centered data on 3D meshes
void run_3D_NODE_order2(std::vector<std::string> source_nodevar_names,
                        std::vector<std::string> target_nodevar_names);


 private:
SourceMesh_Wrapper  const & source_mesh_;
TargetMesh_Wrapper const & target_mesh_;
SourceState_Wrapper const & source_state_;
TargetState_Wrapper const & target_state_;
std::vector<std::string> source_remap_var_names_;
std::vector<std::string> target_remap_var_names_;
unsigned int interp_order_;
unsigned int dim_;
};  // class Driver


  // 1st order remapping of cell centered data on 2D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_2D_CELL_order1(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    // Get an instance of the desired search algorithm type
    const SearchKDTree<2, SourceMesh_Wrapper, TargetMesh_Wrapper>
        search(source_mesh_, target_mesh_);

    // Get an instance of the desired intersect algorithm type
    const IntersectClipper<SourceMesh_Wrapper, TargetMesh_Wrapper>
        intersect{source_mesh_, target_mesh_};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 1st order accurate algorithm" << std::endl;

      const Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i]);

      // Make the remapper instance
      RemapFunctor<SearchKDTree<2, SourceMesh_Wrapper, TargetMesh_Wrapper>,
          IntersectClipper<SourceMesh_Wrapper, TargetMesh_Wrapper>,
          Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper,
          CELL>>
          remapper(&search, &intersect, &interpolater);


#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(CELL, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                         (counting_iterator)(target_mesh_.end(CELL)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }


    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }

  // 2nd order remapping of cell centered data on 2D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_2D_CELL_order2(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    // Get an instance of the desired search algorithm type
    const SearchKDTree<2, SourceMesh_Wrapper, TargetMesh_Wrapper>
        search(source_mesh_, target_mesh_);

    // Get an instance of the desired intersect algorithm type
    const IntersectClipper<SourceMesh_Wrapper, TargetMesh_Wrapper>
        intersect{source_mesh_, target_mesh_};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 2nd order accurate algorithm" << std::endl;

      // Get an instance of the 2nd order remapping algorithm
      const Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i],
                       NOLIMITER);

      // Make the remapper instance
      RemapFunctor<SearchKDTree<2, SourceMesh_Wrapper, TargetMesh_Wrapper>,
          IntersectClipper<SourceMesh_Wrapper, TargetMesh_Wrapper>,
          Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL> >
          remapper(&search, &intersect, &interpolater);

#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(CELL, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                         (counting_iterator)(target_mesh_.end(CELL)),
                         target_field, remapper);
      /* UNCOMMENT WHEN WE RESTORE get_type
         } else {
         std::cerr << "Cannot remap " << source_var_names[i] <<
         " because it is not a scalar double variable\n";
         continue;
         }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }


  // 1st order remapping of cell centered data on 3D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_3D_CELL_order1(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    // Get an instance of the desired search algorithm type
    const SearchKDTree<3, SourceMesh_Wrapper, TargetMesh_Wrapper>
        search(source_mesh_, target_mesh_);

    // Get an instance of the desired intersect algorithm type
    const IntersectR3D<SourceMesh_Wrapper, TargetMesh_Wrapper>
        intersect{source_mesh_, target_mesh_};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 1st order accurate algorithm" << std::endl;

      // Get an instance of the 1st order algorithm
      const Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i]);

      // Make the remapper instance
      RemapFunctor<SearchKDTree<3, SourceMesh_Wrapper, TargetMesh_Wrapper>,
          IntersectR3D<SourceMesh_Wrapper, TargetMesh_Wrapper>,
          Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL> >
          remapper(&search, &intersect, &interpolater);


#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(CELL, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                         (counting_iterator)(target_mesh_.end(CELL)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }

  // Distributed 1st order remapping of cell centered data on 3D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_3D_CELL_order1_distributed(std::vector<std::string>
                                                           source_var_names,
                                                           std::vector<std::string>
                                                           target_var_names) {
    float tot_seconds = 0.0;

    // Get the rank for this process
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    // Convert the source mesh and state to a flat representation;
    // Since we are not sending any target mesh data over MPI, we don't need
    // to convert the target mesh or state to a flat representation
    Flat_Mesh_Wrapper<> source_mesh_flat(8, source_mesh_);
    Flat_State_Wrapper<> source_state_flat(source_state_, source_remap_var_names_);

    // Use a bounding box distributor to send the source cells to the target
    // paritions where they are needed
    MPI_Bounding_Boxes distributor;
    distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_, target_state_);

    // Get an instance of the desired search algorithm type
    const SearchKDTree<3, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>
        search(source_mesh_flat, target_mesh_);

    // Get an instance of the desired intersect algorithm type
    const IntersectR3D<Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>
        intersect{source_mesh_flat, target_mesh_};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      if (comm_rank == 0)
        std::cout << "Remapping variable " << source_var_names[i]
                  << " to variable " << target_var_names[i]
                  << " using a 1st order accurate algorithm" << std::endl;

      // Get an instance of the 1st order algorithm
      const Interpolate_1stOrder<Flat_Mesh_Wrapper<>, TargetMesh_Wrapper,
          Flat_State_Wrapper<>, CELL>
          interpolater(source_mesh_flat, target_mesh_,
                       source_state_flat, source_var_names[i]);

      // Make the remapper instance
      RemapFunctor<SearchKDTree<3, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>,
          IntersectR3D<Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>,
          Interpolate_1stOrder<Flat_Mesh_Wrapper<>, TargetMesh_Wrapper,
          Flat_State_Wrapper<>, CELL> >
          remapper(&search, &intersect, &interpolater);

#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

       // This populates targetField with the values returned by the
       // remapper operator

       /* UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {*/
      
      double *target_field_raw = nullptr;
      target_state_.get_data(CELL, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                         (counting_iterator)(target_mesh_.end(CELL)),
                         target_field, remapper);
      
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }*/
      
#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_rank == 0) std::cout << "Transform Time: " << tot_seconds << std::endl;
  }

  // 2nd order remapping of cell centered data on 3D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_3D_CELL_order2(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    // Get an instance of the desired search algorithm type
    const SearchKDTree<3, SourceMesh_Wrapper, TargetMesh_Wrapper>
        search(source_mesh_, target_mesh_);

    // Get an instance of the desired intersect algorithm type
    const IntersectR3D<SourceMesh_Wrapper, TargetMesh_Wrapper>
        intersect{source_mesh_, target_mesh_};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 2nd order accurate algorithm" << std::endl;

      // Get an instance of the 2nd order algorithm
      const Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i],
                       NOLIMITER);


      RemapFunctor<SearchKDTree<3, SourceMesh_Wrapper, TargetMesh_Wrapper>,
          IntersectR3D<SourceMesh_Wrapper, TargetMesh_Wrapper>,
          Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, CELL> >
          remapper(&search, &intersect, &interpolater);


#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(CELL, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                         (counting_iterator)(target_mesh_.end(CELL)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif
      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }

  // 1st order remapping of node centered data on 2D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_2D_NODE_order1(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    MeshWrapperDual<SourceMesh_Wrapper> source_mesh_dual(source_mesh_);
    MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

    // Get an instance of the desired search algorithm type
    const SearchKDTree<2, MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        search(source_mesh_dual, target_mesh_dual);

    // Get an instance of the desired intersect algorithm type
    const IntersectClipper<MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        intersect{source_mesh_dual, target_mesh_dual};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 1st order accurate algorithm" << std::endl;

      const Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i]);

      RemapFunctor<SearchKDTree<2, MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          IntersectClipper<MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE> >
          remapper(&search, &intersect, &interpolater);


#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(NODE, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                         (counting_iterator)(target_mesh_.end(NODE)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }

  // 2nd order remapping of cell centered data on 2D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_2D_NODE_order2(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    MeshWrapperDual<SourceMesh_Wrapper> source_mesh_dual(source_mesh_);
    MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

    // Get an instance of the desired search algorithm type
    const SearchKDTree<2, MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        search(source_mesh_dual, target_mesh_dual);

    // Get an instance of the desired intersect algorithm type
    const IntersectClipper<MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        intersect{source_mesh_dual, target_mesh_dual};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 2nd order accurate algorithm" << std::endl;

      // Get an instance of the 2nd order interpolate algorithm
      const Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i],
                       NOLIMITER);

      // Make the remapper instance
      RemapFunctor<SearchKDTree<2, MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          IntersectClipper<MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE>>
          remapper(&search, &intersect, &interpolater);

#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(NODE, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                         (counting_iterator)(target_mesh_.end(NODE)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }


  // 1st order remapping of cell centered data on 3D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_3D_NODE_order1(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    MeshWrapperDual<SourceMesh_Wrapper> source_mesh_dual(source_mesh_);
    MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

    // Get an instance of the desired search algorithm type
    const SearchKDTree<3, MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        search(source_mesh_dual, target_mesh_dual);

    // Get an instance of the desired intersect algorithm type
    const IntersectR3D<MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        intersect{source_mesh_dual, target_mesh_dual};

    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 1st order accurate algorithm" << std::endl;

      // Get an instance of the 1st order algorithm
      const Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i]);

      // Make the remapper instance
      RemapFunctor<SearchKDTree<3, MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          IntersectR3D<MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          Interpolate_1stOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE>>
          remapper(&search, &intersect, &interpolater);


#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(NODE, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                         (counting_iterator)(target_mesh_.end(NODE)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }

  // 2nd order remapping of cell centered data on 3D meshes
  template<class SourceMesh_Wrapper, class SourceState_Wrapper,
      class TargetMesh_Wrapper, class TargetState_Wrapper>
      void
      Driver<SourceMesh_Wrapper,
      SourceState_Wrapper,
      TargetMesh_Wrapper,
      TargetState_Wrapper>::run_3D_NODE_order2(std::vector<std::string>
                                               source_var_names,
                                               std::vector<std::string>
                                               target_var_names) {
    float tot_seconds = 0.0;

    MeshWrapperDual<SourceMesh_Wrapper> source_mesh_dual(source_mesh_);
    MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

    // Get an instance of the desired search algorithm type
    const SearchKDTree<3, MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        search(source_mesh_dual, target_mesh_dual);

    // Get an instance of the desired intersect algorithm type
    const IntersectR3D<MeshWrapperDual<SourceMesh_Wrapper>,
        MeshWrapperDual<TargetMesh_Wrapper>>
        intersect{source_mesh_dual, target_mesh_dual};


    int nvars = source_var_names.size();
    for (int i = 0; i < nvars; ++i) {
      std::cout << "Remapping variable " << source_var_names[i]
                << " to variable " << target_var_names[i]
                << " using a 2nd order accurate algorithm" << std::endl;

      // Get an instance of the 2nd order algorithm
      const Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE>
          interpolater(source_mesh_, target_mesh_,
                       source_state_, source_var_names[i],
                       NOLIMITER);


      RemapFunctor<SearchKDTree<3, MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          IntersectR3D<MeshWrapperDual<SourceMesh_Wrapper>,
          MeshWrapperDual<TargetMesh_Wrapper>>,
          Interpolate_2ndOrder<SourceMesh_Wrapper, TargetMesh_Wrapper,
          SourceState_Wrapper, NODE> >
          remapper(&search, &intersect, &interpolater);


#ifdef ENABLE_PROFILE
      __itt_resume();
#endif

      struct timeval begin, end, diff;
      gettimeofday(&begin, 0);

      // This populates targetField with the values returned by the
      // remapper operator

      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          if (typeid(source_state_.get_type(source_var_names[i])) ==
          typeid(double)) {
      */
      double *target_field_raw = nullptr;
      target_state_.get_data(NODE, target_var_names[i], &target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                         (counting_iterator)(target_mesh_.end(NODE)),
                         target_field, remapper);
      /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
          } else {
          std::cerr << "Cannot remap " << source_var_names[i] <<
          " because it is not a scalar double variable\n";
          continue;
          }
      */

#ifdef ENABLE_PROFILE
      __itt_pause();
#endif

      gettimeofday(&end, 0);
      timersub(&end, &begin, &diff);
      float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
      tot_seconds += seconds;
    }

    std::cout << "Transform Time: " << tot_seconds << std::endl;
  }


  /*!
    @struct RemapFunctor "driver.h"
    @brief This functor is used inside a Portage::transform() inside
    Driver::run() to actually do the search, intersect, and interpolation
    calculations.
    @tparam SearchType The type of search method (e.g. SearchSimple or
    SearchKDTree).
    @tparam IsectType The type of intersect method (e.g. IntersectClipper).
    @tparam InterpType The type of interpolation method (e.g. Interpolate_1stOrder
    or Interpolate_2ndOrder).
  */
  template <typename SearchType, typename IsectType, typename InterpType>
      struct RemapFunctor {
        const SearchType* search_;       ///< search method (e.g. SearchSimple)
        const IsectType* intersect_;     ///< intersect method (e.g. IntersectClipper)
        const InterpType* interpolater_;  ///< interpolation method (e.g. Interpolate_2ndOrder)

        /*!
          @brief Constructor.
          @param[in] searcher The search method to use (e.g. SearchSimple)
          @param[in] intersecter The intersect method to use (e.g. IntersectClipper)
          @param[in] interpolater The interpolation method to use (e.g.
          Interpolate_2ndOrder)
        */
       RemapFunctor(const SearchType* searcher,
                    const IsectType* intersecter,
                    const InterpType* interpolater)
       : search_(searcher), intersect_(intersecter),
          interpolater_(interpolater)
        {}

        /*!
          @brief Operator for making this struct a functor

          This is called from within a Portage::transform() operation that iterates
          over the cells in a target mesh.

          @param[in] targetCellindex The cell ID in the target mesh that this functor
          is currently operating on.

          @return Value of the field (stored in @c interpolater ) in the target mesh
          cell with ID @c targetCellIndex.
        */
        double operator()(int const targetCellIndex) {
          // Search for candidates and return their cells indices
          std::vector<int> candidates;
          (*search_)(targetCellIndex, &candidates);

          // Intersect target cell with cells of source mesh and return the
          // moments of intersection
          std::vector<std::vector<std::vector<double>>> moments(candidates.size());
          for (int i = 0; i < candidates.size(); ++i)
            moments[i] = (*intersect_)(candidates[i], targetCellIndex);

          // Compute new value on target cell based on source mesh
          // values and intersection moments

          // Each cell-cell intersection can result in multiple
          // disjointed pieces if one of the cells in non-convex.
          // therefore, there may be more than one set of moments per
          // cell pair. Transform the 3 nested std::vector form to 2
          // nested std::vector form with duplicate candidate entries if
          // need be
          int nalloc = 0;
          for (const auto &moment : moments)
            nalloc += moment.size();

          std::vector<int> candidates_dup(nalloc);
          std::vector< std::vector<double> > interp_moments(nalloc);

          int ninserted = 0;
          for (int i = 0; i < candidates.size(); ++i) {
            std::vector< std::vector<double> > & candidate_moments = moments[i];
            int num_moment_sets = candidate_moments.size();
            for (const auto & candidate_moment : candidate_moments) {
              candidates_dup[ninserted] = candidates[i];  // repeated as needed
              interp_moments[ninserted] = candidate_moment;
              ++ninserted;
            }
          }

          std::pair<std::vector<int> const &,
              std::vector<std::vector<double>> const &>
              source_cells_and_weights(candidates_dup, interp_moments);

          double interpolatedValue = (*interpolater_)(source_cells_and_weights,
                                                      targetCellIndex);

          return interpolatedValue;
        }
      };  // struct RemapFunctor

}  // namespace Portage

#endif  // SRC_DRIVER_DRIVER_H_
