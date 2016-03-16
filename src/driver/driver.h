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
#include "portage/search/search_kdtree2.h"
#include "portage/search/search_kdtree3.h"
#include "portage/intersect/intersectClipper.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"

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

    Sharing of cells is determined from the Parallel_type (e.g. @c OWNED,
    @c GHOST, @c ALL ).

    @param[in] dualcellID The cell ID for which you would like to find the
    neighbors.
    @param[in] ptype The type of data you want (@c OWNED, @c GHOST, @c ALL)
    @param[in,out] adjcells List of IDs of adjacent cells.
  */
  void cell_get_node_adj_cells(int const dualcellID,
                               Parallel_type const ptype,
                               std::vector<int> *adjcells) const {
    w_.node_get_cell_adj_nodes(dualcellID, ptype, adjcells);
  }

  /*!
    @brief Get the IDs of all cells that share a node with the specified cell
    of the <em> dual mesh </em>.

    Sharing of cells is determined from the Parallel_type (e.g. @c OWNED,
    @c GHOST, @c ALL).

    @param[in] dualnodeID The cell ID for which you would like to find the
    neighbors.
    @param[in] ptype The type of data you want (@c OWNED, @c GHOST, @c ALL)
    @param[in,out] adjnodes List of IDs of adjacent cells.

    @todo Clarify this wrt to @c MeshWrapperDual::cell_get_node_adj_cells()
  */
  void dual_cell_get_node_adj_cells(int const dualnodeID,
                                    Parallel_type const ptype,
                                    std::vector<int> *adjnodes) const {
    w_.cell_get_node_adj_cells(dualnodeID, ptype, adjnodes);
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

 private:
  const Mesh_Wrapper_Type &w_;
};


// Forward definition
template <typename SearchType, typename IsectType, typename InterpType>
    struct RemapFunctor;

/*!
  @class Driver "driver.h"
  @brief Driver provides the API to mapping from one mesh to another.
  @tparam InputMesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.  See Jali_Mesh_Wrapper
  for an example.
  @tparam InputState_Wrapper A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.  See
  Jali_State_Wrapper for an example.
  @tparam TargetMesh_Wrapper A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.  See Jali_Mesh_Wrapper
  for an example.
  @tparam TargetState_Wrapper A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.  See
  Jali_State_Wrapper for an example.

*/
template <class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper=InputMesh_Wrapper,
    class TargetState_Wrapper=InputState_Wrapper>
class Driver {
 public:
  /*!
    @brief Constructor for running the interpolation driver.
    @param[in] remapEntity  The type of entity the remapping is on (CELL, NODE)
    @param[in] sourceMesh A @c InputMesh_Wrapper to the source mesh.
    @param[in] sourceState A @c InputState_Wrapperfor the data that lives on the
    source mesh.
    @param[in] targetMesh A @c TargetMesh_Wrapper to the target mesh.
    @param[in,out] targetState A @c TargetState_Wrapper for the data that will
    be mapped to the target mesh.
  */
  Driver(Entity_kind const remapEntity,
         InputMesh_Wrapper const & sourceMesh,
         InputState_Wrapper const & sourceState,
         TargetMesh_Wrapper const & targetMesh,
         TargetState_Wrapper const & targetState)
      : remap_entity_(remapEntity),
      source_mesh_(sourceMesh), source_state_(sourceState),
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
    src_remap_var_names_ = remap_var_names;
    tar_remap_var_names_ = remap_var_names;
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] src_remap_var_names A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] tar_remap_var_names  A list of the variables names of the
    variables to interpolate to the target mesh.
   */
  void set_remap_var_names(
      std::vector<std::string> const &src_remap_var_names,
      std::vector<std::string> const &tar_remap_var_names) {
    assert(src_remap_var_names.size() == tar_remap_var_names.size());
    src_remap_var_names_ = src_remap_var_names;
    tar_remap_var_names_ = tar_remap_var_names;
  }

  /*!
    @brief Get the names of the variables to be remapped from the
    source mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> src_remap_var_names() const {
    return src_remap_var_names_;
  }

  /*!
    @brief Get the names of the variables to be remapped to the
    target mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> tar_remap_var_names() const {
    return tar_remap_var_names_;
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

    // Get raw pointer to the data from the state manager
    double *target_field_raw = nullptr;
    target_state_.get_data(remap_entity_, tar_remap_var_names_[0],
                           &target_field_raw);
    Portage::pointer<double> target_field(target_field_raw);

    // @todo We might be able to make this less verbose using boost::mpl

    switch (dim_) {
      case 1:
        std::cerr << "Remapping not implemented for 1D" << std::endl;
        exit(-1);
      case 2: {
        switch (remap_entity_) {
          case CELL: {
            (interp_order_ == 1) ?
                run_2D_CELL_order1(target_field) :
                run_2D_CELL_order2(target_field);
            break;
          }
          case NODE: {
            (interp_order_ == 1) ?
                run_2D_NODE_order1(target_field) :
                run_2D_NODE_order2(target_field);
            break;
          }
          default:
            std::cerr << "Remapping only implemented for CELLS and NODES"
                      << std::endl;
            exit(-1);
        }
        break;
      }
      case 3: {
        switch (remap_entity_) {
          case CELL: {
            (interp_order_ == 1) ?
                run_3D_CELL_order1(target_field) :
                run_3D_CELL_order2(target_field);
            break;
          }
          case NODE: {
            (interp_order_ == 1) ?
                run_3D_NODE_order1(target_field) :
                run_3D_NODE_order2(target_field);
            break;
          }
          default: {
            std::cerr << "Remapping only implemented for CELLS and NODES"
                      << std::endl;
            exit(-1);
          }
        }
        break;
      }
      default:
        std::cerr << "Invalid dimension" << std::endl;
        exit(-1);
    }
  }


  /// @brief 1st order remapping of cell centered data on 2D meshes
  void run_2D_CELL_order1(Portage::pointer<double> target_field);
  /// @brief 2nd order remapping of cell centered data on 2D meshes
  void run_2D_CELL_order2(Portage::pointer<double> target_field);
  /// @brief 1st order remapping of cell centered data on 3D meshes
  void run_3D_CELL_order1(Portage::pointer<double> target_field);
  /// @brief 2nd order remapping of cell centered data on 3D meshes
  void run_3D_CELL_order2(Portage::pointer<double> target_field);
  /// @brief 1st order remapping of node centered data on 2D meshes
  void run_2D_NODE_order1(Portage::pointer<double> target_field);
  /// @brief 2nd order remapping of node centered data on 2D meshes
  void run_2D_NODE_order2(Portage::pointer<double> target_field);
  /// @brief 1st order remapping of node centered data on 3D meshes
  void run_3D_NODE_order1(Portage::pointer<double> target_field);
  /// @brief 2nd order remapping of node centered data on 3D meshes
  void run_3D_NODE_order2(Portage::pointer<double> target_field);


 private:
  InputMesh_Wrapper  const & source_mesh_;
  TargetMesh_Wrapper const & target_mesh_;
  InputState_Wrapper const & source_state_;
  TargetState_Wrapper const & target_state_;
  std::vector<std::string> src_remap_var_names_;
  std::vector<std::string> tar_remap_var_names_;
  Entity_kind const remap_entity_;
  unsigned int interp_order_;
  unsigned int dim_;
};  // class Driver


// 1st order remapping of cell centered data on 2D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
    Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_2D_CELL_order1(Portage::pointer<double>
                                             target_field) {
  // Get an instance of the desired search algorithm type
  const SearchKDTree2<InputMesh_Wrapper, TargetMesh_Wrapper>
      search(source_mesh_, target_mesh_);

  // Get an instance of the desired intersect algorithm type
  const IntersectClipper<InputMesh_Wrapper, TargetMesh_Wrapper>
      intersect{source_mesh_, target_mesh_};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 1st order accurate algorithm" << std::endl;

  // Eventually put this in a loop over remapping variable names as well

  const Interpolate_1stOrder<InputMesh_Wrapper, InputState_Wrapper, CELL>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0]);

  RemapFunctor<SearchKDTree2<InputMesh_Wrapper, TargetMesh_Wrapper>,
               IntersectClipper<InputMesh_Wrapper, TargetMesh_Wrapper>,
               Interpolate_1stOrder<InputMesh_Wrapper,
                                    InputState_Wrapper, CELL> >
  remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                     (counting_iterator)(target_mesh_.end(CELL)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}

// 2nd order remapping of cell centered data on 2D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_2D_CELL_order2(Portage::pointer<double>
                                             target_field) {
  // Get an instance of the desired search algorithm type
  const SearchKDTree2<InputMesh_Wrapper, TargetMesh_Wrapper>
      search(source_mesh_, target_mesh_);

  // Get an instance of the desired intersect algorithm type
  const IntersectClipper<InputMesh_Wrapper, TargetMesh_Wrapper>
      intersect{source_mesh_, target_mesh_};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 2nd order accurate algorithm" << std::endl;

  // Get an instance of the 2nd order remapping algorithm
  const Interpolate_2ndOrder<InputMesh_Wrapper, InputState_Wrapper, CELL>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0],
            NOLIMITER);

  // Make the remapper instance
  RemapFunctor<SearchKDTree2<InputMesh_Wrapper, TargetMesh_Wrapper>,
               IntersectClipper<InputMesh_Wrapper, TargetMesh_Wrapper>,
               Interpolate_2ndOrder<InputMesh_Wrapper,
                                    InputState_Wrapper, CELL> >
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                     (counting_iterator)(target_mesh_.end(CELL)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}

// 1st order remapping of cell centered data on 3D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_3D_CELL_order1(Portage::pointer<double>
                                             target_field) {
  // Get an instance of the desired search algorithm type
  const SearchKDTree3<InputMesh_Wrapper, TargetMesh_Wrapper>
      search(source_mesh_, target_mesh_);

  // Get an instance of the desired intersect algorithm type
  const IntersectR3D<InputMesh_Wrapper, TargetMesh_Wrapper>
      intersect{source_mesh_, target_mesh_};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 1st order accurate algorithm" << std::endl;

  // Get an instance of the 1st order algorithm
  const Interpolate_1stOrder<InputMesh_Wrapper, InputState_Wrapper, CELL>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0]);

  // Make the remapper instance
  RemapFunctor<SearchKDTree3<InputMesh_Wrapper, TargetMesh_Wrapper>,
               IntersectR3D<InputMesh_Wrapper, TargetMesh_Wrapper>,
               Interpolate_1stOrder<InputMesh_Wrapper,
                                    InputState_Wrapper, CELL> >
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                     (counting_iterator)(target_mesh_.end(CELL)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}

// 2nd order remapping of cell centered data on 3D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_3D_CELL_order2(Portage::pointer<double>
                                             target_field) {
  // Get an instance of the desired search algorithm type
  const SearchKDTree3<InputMesh_Wrapper, TargetMesh_Wrapper>
      search(source_mesh_, target_mesh_);

  // Get an instance of the desired intersect algorithm type
  const IntersectR3D<InputMesh_Wrapper, TargetMesh_Wrapper>
      intersect{source_mesh_, target_mesh_};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 2nd order accurate algorithm" << std::endl;

  // Get an instance of the 2nd order algorithm
  const Interpolate_2ndOrder<InputMesh_Wrapper, InputState_Wrapper, CELL>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0],
            NOLIMITER);


  RemapFunctor<SearchKDTree3<InputMesh_Wrapper, TargetMesh_Wrapper>,
               IntersectR3D<InputMesh_Wrapper, TargetMesh_Wrapper>,
               Interpolate_2ndOrder<InputMesh_Wrapper,
                                    InputState_Wrapper, CELL> >
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
                     (counting_iterator)(target_mesh_.end(CELL)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}

// 1st order remapping of node centered data on 2D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_2D_NODE_order1(Portage::pointer<double>
                                             target_field) {
  MeshWrapperDual<InputMesh_Wrapper> source_mesh_dual(source_mesh_);
  MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

  // Get an instance of the desired search algorithm type
  const SearchKDTree2<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      search(source_mesh_dual, target_mesh_dual);

  // Get an instance of the desired intersect algorithm type
  const IntersectClipper<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      intersect{source_mesh_dual, target_mesh_dual};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 1st order accurate algorithm" << std::endl;

  // Eventually put this in a loop over remapped variable names as well

  const Interpolate_1stOrder<InputMesh_Wrapper, InputState_Wrapper, NODE>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0]);

  RemapFunctor<SearchKDTree2<MeshWrapperDual<InputMesh_Wrapper>,
                             MeshWrapperDual<TargetMesh_Wrapper>>,
      IntersectClipper<MeshWrapperDual<InputMesh_Wrapper>,
                       MeshWrapperDual<TargetMesh_Wrapper>>,
      Interpolate_1stOrder<InputMesh_Wrapper, InputState_Wrapper, NODE> >
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                     (counting_iterator)(target_mesh_.end(NODE)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}

// 2nd order remapping of cell centered data on 2D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_2D_NODE_order2(Portage::pointer<double>
                                             target_field) {
  MeshWrapperDual<InputMesh_Wrapper> source_mesh_dual(source_mesh_);
  MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

  // Get an instance of the desired search algorithm type
  const SearchKDTree2<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      search(source_mesh_dual, target_mesh_dual);

  // Get an instance of the desired intersect algorithm type
  const IntersectClipper<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      intersect{source_mesh_dual, target_mesh_dual};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 2nd order accurate algorithm" << std::endl;

  // Get an instance of the 2nd order interpolate algorithm
  const Interpolate_2ndOrder<InputMesh_Wrapper, InputState_Wrapper, NODE>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0],
            NOLIMITER);

  // Make the remapper instance
  RemapFunctor<SearchKDTree2<MeshWrapperDual<InputMesh_Wrapper>,
                             MeshWrapperDual<TargetMesh_Wrapper>>,
      IntersectClipper<MeshWrapperDual<InputMesh_Wrapper>,
                       MeshWrapperDual<TargetMesh_Wrapper>>,
      Interpolate_2ndOrder<InputMesh_Wrapper, InputState_Wrapper, NODE> >
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                     (counting_iterator)(target_mesh_.end(NODE)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}


// 1st order remapping of cell centered data on 3D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_3D_NODE_order1(Portage::pointer<double>
                                             target_field) {
  MeshWrapperDual<InputMesh_Wrapper> source_mesh_dual(source_mesh_);
  MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

  // Get an instance of the desired search algorithm type
  const SearchKDTree3<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      search(source_mesh_dual, target_mesh_dual);

  // Get an instance of the desired intersect algorithm type
  const IntersectR3D<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      intersect{source_mesh_dual, target_mesh_dual};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 1st order accurate algorithm" << std::endl;

  // Get an instance of the 1st order algorithm
  const Interpolate_1stOrder<InputMesh_Wrapper, InputState_Wrapper, NODE>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0]);

  // Make the remapper instance
  RemapFunctor<SearchKDTree3<MeshWrapperDual<InputMesh_Wrapper>,
                             MeshWrapperDual<TargetMesh_Wrapper>>,
      IntersectR3D<MeshWrapperDual<InputMesh_Wrapper>,
                   MeshWrapperDual<TargetMesh_Wrapper>>,
      Interpolate_1stOrder<InputMesh_Wrapper, InputState_Wrapper, NODE>>
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                     (counting_iterator)(target_mesh_.end(NODE)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}

// 2nd order remapping of cell centered data on 3D meshes
template<class InputMesh_Wrapper, class InputState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
void
Driver<InputMesh_Wrapper,
    InputState_Wrapper,
    TargetMesh_Wrapper,
    TargetState_Wrapper>::run_3D_NODE_order2(Portage::pointer<double>
                                             target_field) {
  MeshWrapperDual<InputMesh_Wrapper> source_mesh_dual(source_mesh_);
  MeshWrapperDual<TargetMesh_Wrapper> target_mesh_dual(target_mesh_);

  // Get an instance of the desired search algorithm type
  const SearchKDTree3<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      search(source_mesh_dual, target_mesh_dual);

  // Get an instance of the desired intersect algorithm type
  const IntersectR3D<MeshWrapperDual<InputMesh_Wrapper>,
      MeshWrapperDual<TargetMesh_Wrapper>>
      intersect{source_mesh_dual, target_mesh_dual};

  std::cout << "Remapping variable " << src_remap_var_names_[0]
            << " to variable " << tar_remap_var_names_[0]
            << " using a 2nd order accurate algorithm" << std::endl;

  // Get an instance of the 2nd order algorithm
  const Interpolate_2ndOrder<InputMesh_Wrapper, InputState_Wrapper, NODE>
      interpolater(source_mesh_, source_state_, src_remap_var_names_[0],
            NOLIMITER);


  RemapFunctor<SearchKDTree3<MeshWrapperDual<InputMesh_Wrapper>,
                             MeshWrapperDual<TargetMesh_Wrapper>>,
      IntersectR3D<MeshWrapperDual<InputMesh_Wrapper>,
                   MeshWrapperDual<TargetMesh_Wrapper>>,
      Interpolate_2ndOrder<InputMesh_Wrapper, InputState_Wrapper, NODE> >
      remapper(&search, &intersect, &interpolater);

  // This populates targetField with the doubles returned from
  // the final remapping

#ifdef ENABLE_PROFILE
  __itt_resume();
#endif

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  Portage::transform((counting_iterator)(target_mesh_.begin(NODE)),
                     (counting_iterator)(target_mesh_.end(NODE)),
                     target_field, remapper);

#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Transform Time: " << seconds << std::endl;
}


/*!
  @struct RemapFunctor "driver.h"
  @brief This functor is used inside a Portage::transform() inside
  Driver::run() to actually do the search, intersect, and interpolation
  calculations.
  @tparam SearchType The type of search method (e.g. SearchSimple or
  SearchKDTree3).
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
  : search_(searcher), intersect_(intersecter), interpolater_(interpolater) { }

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

    std::pair< std::vector<int> const &,
        std::vector< std::vector<double> > const & >
        source_cells_and_weights(candidates_dup, interp_moments);

    double interpolatedValue = (*interpolater_)(source_cells_and_weights);

    return interpolatedValue;
  }
};  // struct RemapFunctor

}  // namespace Portage

#endif  // SRC_DRIVER_DRIVER_H_

