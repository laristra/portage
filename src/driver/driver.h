/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



#ifndef SRC_DRIVER_DRIVER_H_
#define SRC_DRIVER_DRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"

#ifdef ENABLE_MPI
#include "portage/distributed/mpi_bounding_boxes.h"
#endif

/*!
  @file driver.h
  @brief Example driver for mapping between two meshes.

  This should serve as a good example for how to write your own driver routine
  and datastructures.
*/

namespace Portage {

  /*!
    @class MeshWrapperDual "driver.h"
    @brief Wrapper for dual mesh.

    Utilizes a Mesh_Wrapper to the original mesh, but treats
    the nodes of the original mesh as the centroids of the dual mesh.
  */
template<class Mesh_Wrapper_Type>
class MeshWrapperDual {  // cellid is the dual cell (i.e. node) id
 public:
  /*!
    @brief Constructor of a wrapper to a 2d mesh.
    @param[in] w Mesh_Wrapper to original mesh.
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

    Sharing of cells is determined from the Entity_type (e.g. @c PARALLEL_OWNED,
    @c GHOST, @c ALL ).

    @param[in] dualcellID The cell ID for which you would like to find the
    neighbors.
    @param[in] type The type of data you want (@c PARALLEL_OWNED, @c GHOST, @c ALL)
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

    Sharing of cells is determined from the Entity_type (e.g. @c PARALLEL_OWNED,
    @c GHOST, @c ALL).

    @param[in] dualnodeID The cell ID for which you would like to find the
    neighbors.
    @param[in] type The type of data you want (@c PARALLEL_OWNED, @c GHOST, @c ALL)
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
  template <long D>
  void cell_centroid(int const dualcellID, Point<D> *centroid)
      const {
    w_.dual_cell_centroid(dualcellID, centroid);
  }

  /*!
    @brief Get the coordinates of the centroid of a given dual mesh cell ID.
    @param[in] dualcellID ID of the cell in the dual mesh.
    @param[in,out] centroid (x,y,z) coordinates of the cell center (for 3d).
    @todo Clarify this wrt to @c MeshWrapperDual::cell_centroid().
  */
  template<long D>
  void dual_cell_centroid(int const dualnodeID, Point<D> *centroid)
      const {
    w_.cell_centroid(dualnodeID, centroid);
  }

  /*!
    @brief Get the coordinates of the points that make up the wedge.

    A wedge corresponds to a cell center, a face center, and two points that
    share an edge.

    @param[in] dualcellid ID of the cell in the dual mesh.
    @param[in,out] wcoords (x,y,z) coordinates of each of the four points that
    comprise the tetrahedron that is the wedge.
  */
  template <int D>
  void wedges_get_coordinates(int const dualcellid,
              std::vector<std::array<Portage::Point<D>, D+1>> *wcoords) const {
    w_.dual_wedges_get_coordinates(dualcellid, wcoords);
  }

  //! Get a decomposition of a 2D dual cell into simplices
  void decompose_cell_into_simplices(int dualcellid,
               std::vector<std::array<Portage::Point<2>, 3>> *tcoords) const {
    decompose_cell_into_tris(dualcellid, tcoords);
  }

  //! Get a decomposition of a 3D dual cell into simplices
  void decompose_cell_into_simplices(int dualcellid,
               std::vector<std::array<Portage::Point<3>, 4>> *tcoords) const {
    // needs us to say if it is a planar hex or not - in general, its not
    decompose_cell_into_tets(dualcellid, tcoords, false);
  }

  //! Get a decomposition of a 2D cell into tris
  // For a dual mesh, that means returning a list of wedges
  void decompose_cell_into_tris(int const dualcellid,
               std::vector<std::array<Portage::Point<2>, 3>> *tcoords) const {
    wedges_get_coordinates<2>(dualcellid, tcoords);
  }

  // Get the simplest possible decomposition of a 3D cell into tets.
  // For a dual mesh, that means returning a list of wedges.
  void decompose_cell_into_tets(int const dualcellid,
               std::vector<std::array<Portage::Point<3>, 4>> *tcoords,
                                const bool planar_hex) const {
    wedges_get_coordinates<3>(dualcellid, tcoords);
  }

 private:
  const Mesh_Wrapper_Type &w_;
};

// Forward definitions
template <typename SearchType> struct SearchFunctor;
template <int D, typename IntersectType,
          typename SourceMeshWrapper,
          typename TargetMeshWrapper> struct IntersectFunctor;

/*!
  @class Driver "driver.h"
  @brief Driver provides the API to mapping from one mesh to another.
  @tparam SourceMesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.
  @tparam SourceState_Wrapper A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.
  @tparam TargetMesh_Wrapper A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.
  @tparam TargetState_Wrapper A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.
*/
template <template <int, class, class> class Search,
          class Intersect,
          template<class, class, class, Entity_kind, long> class Interpolate,
          int Dim,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper = SourceMesh_Wrapper,
          class TargetState_Wrapper = SourceState_Wrapper>
class Driver {

  // Something like this would be very helpful to users
  // static_assert(
  //   Dim == Interpolate::Dim,
  //   "The dimension of Driver and Interpolate do not match!"
  // );


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
  Driver(SourceMesh_Wrapper const& sourceMesh,
         SourceState_Wrapper const& sourceState,
         TargetMesh_Wrapper const& targetMesh,
         TargetState_Wrapper& targetState)
      : source_mesh_(sourceMesh), source_state_(sourceState),
        target_mesh_(targetMesh), target_state_(targetState),
        dim_(sourceMesh.space_dimension()) {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
  }

  /// Copy constructor (disabled)
  Driver(const Driver &) = delete;

  /// Assignment operator (disabled)
  Driver & operator = (const Driver &) = delete;

  /// Destructor
  ~Driver() {}

  /*!
    @brief Specify the names of the variables to be interpolated along with the
    limiter to use for all of them
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source mesh to the target mesh.  This variable must
    exist in both meshes' state manager
    @param[in] limiter_type The type of limiter to use (NOLIMITER, BARTH_JESPERSEN)
  */
  void set_remap_var_names(std::vector<std::string> const &remap_var_names,
                           LimiterType limiter_type = NOLIMITER) {
    // All variables use the same type of limiters
    std::vector<LimiterType> limiters(remap_var_names.size(), limiter_type);

    // remap variable names same in source and target mesh
    set_remap_var_names(remap_var_names, remap_var_names, limiters);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target mesh.
    @param[in] limiter_type The limiter to use for higher order remaps (NOLIMITER, BARTH_JESPERSEN)
  */
  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      LimiterType limiter_type = NOLIMITER) {
    std::vector<LimiterType> limiters(source_remap_var_names.size(),
                                      limiter_type);
    
    set_remap_var_names(source_remap_var_names, target_remap_var_names,
                        limiters);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target mesh.
    @param[in] limiter_types Limiters to use for each remapped variable (NOLIMITER, BARTH_JESPERSEN)
  */

  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      std::vector<LimiterType> const & limiter_types) {
    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    for (int i = 0; i < nvars; ++i)
      assert(source_state_.get_entity(source_remap_var_names[i]) ==
             target_state_.get_entity(target_remap_var_names[i]));

    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    limiters_ = limiter_types;
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

  /*!
    @brief Get the dimensionality of the meshes.
    @return The dimensionality of the meshes.
  */
  unsigned int dim() const {
    return dim_;
  }

  /*!
    @brief Execute the remapping process
  */
  void run(bool distributed) {

#ifndef ENABLE_MPI
    if (distributed) {
      std::cout << "Request is for a parallel run but Portage is compiled for serial runs only\n";
      return;
    }
#endif


    int comm_rank = 0;

#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

    if (comm_rank == 0) std::printf("in Driver::run()...\n");

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;

    int nvars = source_remap_var_names_.size();

    Flat_Mesh_Wrapper<> source_mesh_flat;
    Flat_State_Wrapper<> source_state_flat;

    // Collect all cell based variables and remap them

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

    if (source_cellvar_names.size() > 0)
    {

      float tot_seconds = 0.0, tot_seconds_srch = 0.0,
          tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
      struct timeval begin_timeval, end_timeval, diff_timeval;

      int ntargetcells = target_mesh_.num_entities(CELL);

      // SEARCH
      
      Portage::vector<std::vector<int>> candidates(ntargetcells);
      Portage::vector<std::vector<Weights_t>> source_cells_and_weights(ntargetcells);

      if (distributed) {

#ifdef ENABLE_MPI
        // Our current flecsi build does not support distributed meshes,
        // so in that case don't try to build or run this code.

        // Create flat wrappers to distribute source cells 
        gettimeofday(&begin_timeval, 0);

        source_mesh_flat.initialize(source_mesh_);
        source_state_flat.initialize(source_state_, source_remap_var_names_);
        MPI_Bounding_Boxes distributor;
        distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                               target_state_);

        gettimeofday(&end_timeval, 0);
        timersub(&end_timeval, &begin_timeval, &diff_timeval);
        float tot_seconds_flat = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
        std::cout << "Redistribution Time Rank " << comm_rank << " (s): " << tot_seconds_flat << std::endl;  

        // Get an instance of the desired search algorithm type
        gettimeofday(&begin_timeval, 0);
        const Search<Dim, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper> search(source_mesh_flat, target_mesh_);
        SearchFunctor<Search<Dim, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>> searchfunctor(&search);

        Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                           target_mesh_.end(CELL, PARALLEL_OWNED),
                           candidates.begin(), searchfunctor);
#endif
      }
      else {

        // Get an instance of the desired search algorithm type
        gettimeofday(&begin_timeval, 0);
        const Search<Dim, SourceMesh_Wrapper, TargetMesh_Wrapper> search(source_mesh_, target_mesh_);
        SearchFunctor<Search<Dim, SourceMesh_Wrapper, TargetMesh_Wrapper>>
            searchfunctor(&search);

        Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                           target_mesh_.end(CELL, PARALLEL_OWNED),
                           candidates.begin(), searchfunctor);
      }

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      // INTERSECT

      gettimeofday(&begin_timeval, 0);

      if (distributed) {

#ifdef ENABLE_MPI
        // Get an instance of the desired intersect algorithm type
        const Intersect intersect;
        IntersectFunctor<Dim, Intersect,
                         Flat_Mesh_Wrapper<>,
                         TargetMesh_Wrapper>
            intersectfunctor(&intersect, source_mesh_flat, target_mesh_);

        Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                           target_mesh_.end(CELL, PARALLEL_OWNED),
                           candidates.begin(),
                           source_cells_and_weights.begin(),
                           intersectfunctor);
#endif
      }
      else {

        // Get an instance of the desired intersect algorithm type
        const Intersect intersect;
        IntersectFunctor<Dim, Intersect, SourceMesh_Wrapper, TargetMesh_Wrapper>
            intersectfunctor(&intersect, source_mesh_, target_mesh_);

        // For each cell in the target mesh get a list of candidate-weight
        // pairings (in a traditional mesh, not particle mesh, the weights
        // are moments). Note that this candidate list is different from the
        // search candidate list in that (1) it may not include some
        // candidates and (2) some candidates may occur twice to account for
        // the fact that the intersection of two cells is more than one
        // disjoint piece (if one of the cells is non-convex). Also, note
        // that for 2nd order and higher remaps, we get multiple moments
        // (0th, 1st, etc) for each target-source cell intersection

        Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                           target_mesh_.end(CELL, PARALLEL_OWNED),
                           candidates.begin(),
                           source_cells_and_weights.begin(),
                           intersectfunctor);
      }

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_xsect = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      // INTERPOLATE (one variable at a time)

      gettimeofday(&begin_timeval, 0);

      nvars = source_cellvar_names.size();
      if (comm_rank == 0) std::cout << "number of cell variables to remap is " << nvars << std::endl;

      if (distributed) {

#ifdef ENABLE_MPI
        // Get an instance of the desired interpolate algorithm type
        Interpolate<Flat_Mesh_Wrapper<>, TargetMesh_Wrapper, Flat_State_Wrapper<>, CELL, Dim>
            interpolate(source_mesh_flat, target_mesh_, source_state_flat);

        for (int i = 0; i < nvars; ++i) {
          interpolate.set_interpolation_variable(source_cellvar_names[i],
                                                   limiters_[i]);

          double *target_field_raw = nullptr;
          target_state_.get_data(CELL, target_cellvar_names[i], &target_field_raw);
          Portage::pointer<double> target_field(target_field_raw);

          Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                             target_mesh_.end(CELL, PARALLEL_OWNED),
                             source_cells_and_weights.begin(),
                             target_field, interpolate);
        }
#endif
      }
      else {

        // Get an instance of the desired interpolate algorithm type
        Interpolate<SourceMesh_Wrapper, TargetMesh_Wrapper, SourceState_Wrapper, CELL, Dim>
            interpolate(source_mesh_, target_mesh_, source_state_);
        
        for (int i = 0; i < nvars; ++i) {
          //amh: ?? add back accuracy output statement??
          if (comm_rank == 0) std::cout << "Remapping cell variable " << source_cellvar_names[i]
                                        << " to variable " << target_cellvar_names[i] << std::endl;
          interpolate.set_interpolation_variable(source_cellvar_names[i],
                                                   limiters_[i]);

          // This populates targetField with the values returned by the
          // remapper operator
          
          /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
              if (typeid(source_state_.get_type(source_var_names[i])) ==
              typeid(double)) {
          */
          double *target_field_raw = nullptr;
          target_state_.get_data(CELL, target_cellvar_names[i], &target_field_raw);
          Portage::pointer<double> target_field(target_field_raw);
          
          Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                             target_mesh_.end(CELL, PARALLEL_OWNED),
                             source_cells_and_weights.begin(),
                             target_field, interpolate);
          /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
              } else {
              std::cerr << "Cannot remap " << source_var_names[i] <<
              " because it is not a scalar double variable\n";
              continue;
              }
          */
        }
      }
      
      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_interp = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

      std::cout << "Transform Time Rank " << comm_rank << " (s): " <<
          tot_seconds << std::endl;
      std::cout << "  Search Time Rank " << comm_rank << " (s): " <<
          tot_seconds_srch << std::endl;
      std::cout << "  Intersect Time Rank " << comm_rank << " (s): " <<
          tot_seconds_xsect << std::endl;
      std::cout << "  Interpolate Time Rank " << comm_rank << " (s): " <<
          tot_seconds_interp << std::endl;
    }

    // Collect all node based variables and remap them
    MeshWrapperDual<SourceMesh_Wrapper> sourceDualWrapper(source_mesh_);
    MeshWrapperDual<TargetMesh_Wrapper> targetDualWrapper(target_mesh_);

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

      float tot_seconds = 0.0, tot_seconds_srch = 0.0,
            tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
      struct timeval begin_timeval, end_timeval, diff_timeval;

      int ntargetnodes = target_mesh_.num_entities(NODE);

      // SEARCH

      Portage::vector<std::vector<int>> candidates(ntargetnodes);
      Portage::vector<std::vector<Weights_t>> source_cells_and_weights(ntargetnodes);

      MeshWrapperDual<Flat_Mesh_Wrapper<>> sourceDualFlat(source_mesh_flat);

      if (distributed) {
#ifndef PORTAGE_SERIAL_ONLY
        // Create flat wrappers to distribute source cells
        gettimeofday(&begin_timeval, 0);

        source_mesh_flat.initialize(source_mesh_);
        source_state_flat.initialize(source_state_, source_remap_var_names_);
        MPI_Bounding_Boxes distributor;
        distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                               target_state_);
        gettimeofday(&end_timeval, 0);
        timersub(&end_timeval, &begin_timeval, &diff_timeval);
        float tot_seconds_flat = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
        std::cout << "Redistribution Time Rank " << comm_rank << " (s): " << tot_seconds_flat << std::endl;

        // Get an instance of the desired search algorithm type
        gettimeofday(&begin_timeval, 0);
        const Search<Dim, MeshWrapperDual<Flat_Mesh_Wrapper<>>, MeshWrapperDual<TargetMesh_Wrapper>> search(sourceDualFlat, targetDualWrapper);
        SearchFunctor<Search<Dim, MeshWrapperDual<Flat_Mesh_Wrapper<>>, MeshWrapperDual<TargetMesh_Wrapper>>> searchfunctor(&search);

        Portage::transform(target_mesh_.begin(NODE, PARALLEL_OWNED),
                           target_mesh_.end(NODE, PARALLEL_OWNED),
                           candidates.begin(), searchfunctor);
#endif
      }
      else {

        gettimeofday(&begin_timeval, 0);

        // Get an instance of the desired search algorithm type
        const Search<Dim, MeshWrapperDual<SourceMesh_Wrapper>,
                     MeshWrapperDual<TargetMesh_Wrapper>>
              search(sourceDualWrapper, targetDualWrapper);
        SearchFunctor<Search<Dim, MeshWrapperDual<SourceMesh_Wrapper>,
                             MeshWrapperDual<TargetMesh_Wrapper>>>
            searchfunctor(&search);
        Portage::transform(target_mesh_.begin(NODE, PARALLEL_OWNED),
                           target_mesh_.end(NODE, PARALLEL_OWNED),
                           candidates.begin(), searchfunctor);
      }

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      // INTERSECT

      gettimeofday(&begin_timeval, 0);

      if (distributed) {

#ifndef PORTAGE_SERIAL_ONLY
        // Get an instance of the desired intersect algorithm type
        const Intersect intersect;
        IntersectFunctor<Dim, Intersect,
                         MeshWrapperDual<Flat_Mesh_Wrapper<>>,
                         MeshWrapperDual<TargetMesh_Wrapper>>
            intersectfunctor(&intersect, sourceDualFlat, targetDualWrapper);

        Portage::transform(target_mesh_.begin(NODE, PARALLEL_OWNED),
                           target_mesh_.end(NODE, PARALLEL_OWNED),
                           candidates.begin(),
                           source_cells_and_weights.begin(),
                           intersectfunctor);
#endif
      }
      else {

        // Get an instance of the desired intersect algorithm type
        const Intersect intersect;
        IntersectFunctor<Dim, Intersect,
                         MeshWrapperDual<SourceMesh_Wrapper>,
                         MeshWrapperDual<TargetMesh_Wrapper>>
            intersectfunctor(&intersect, sourceDualWrapper, targetDualWrapper);

        // For each cell in the target mesh get a list of candidate-weight
        // pairings (in a traditional mesh, not particle mesh, the weights
        // are moments). Note that this candidate list is different from the
        // search candidate list in that (1) it may not include some
        // candidates and (2) some candidates may occur twice to account for
        // the fact that the intersection of two cells is more than one
        // disjoint piece (if one of the cells is non-convex). Also, note
        // that for 2nd order and higher remaps, we get multiple moments
        // (0th, 1st, etc) for each target-source cell intersection

        Portage::transform(target_mesh_.begin(NODE, PARALLEL_OWNED),
                           target_mesh_.end(NODE, PARALLEL_OWNED),
                           candidates.begin(),
                           source_cells_and_weights.begin(),
                           intersectfunctor);
      }

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_xsect = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      // INTERPOLATE (one variable at a time)

      gettimeofday(&begin_timeval, 0);

      nvars = source_nodevar_names.size();
      if (comm_rank == 0) std::cout << "number of node variables to remap is " << nvars << std::endl;

      if (distributed) {

#ifdef ENABLE_MPI
        Interpolate<Flat_Mesh_Wrapper<>, TargetMesh_Wrapper, Flat_State_Wrapper<>, NODE, Dim>
            interpolate(source_mesh_flat, target_mesh_, source_state_flat);

        for (int i = 0; i < nvars; ++i) {
          interpolate.set_interpolation_variable(source_nodevar_names[i],
                                                 limiters_[i]);

          double *target_field_raw = nullptr;
          target_state_.get_data(NODE, target_nodevar_names[i], &target_field_raw);
          Portage::pointer<double> target_field(target_field_raw);

          Portage::transform(target_mesh_.begin(NODE, PARALLEL_OWNED),
                             target_mesh_.end(NODE, PARALLEL_OWNED),
                             source_cells_and_weights.begin(),
                             target_field, interpolate);
        }
#endif
      }
      else {

        Interpolate<SourceMesh_Wrapper, TargetMesh_Wrapper,
                    SourceState_Wrapper, NODE, Dim>
              interpolate(source_mesh_, target_mesh_, source_state_);

        for (int i = 0; i < nvars; ++i) {
          if (comm_rank == 0) std::cout << "Remapping node variable " << source_nodevar_names[i]
                 << " to variable " << target_nodevar_names[i] << std::endl;

          Interpolate<SourceMesh_Wrapper, TargetMesh_Wrapper,
                      SourceState_Wrapper, NODE, Dim>
                interpolate(source_mesh_, target_mesh_, source_state_);

          interpolate.set_interpolation_variable(source_nodevar_names[i],
                                                 limiters_[i]);

          // This populates targetField with the values returned by the
          // interpolate operator

          /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
              if (typeid(source_state_.get_type(source_var_names[i])) ==
              typeid(double)) {
          */
          double *target_field_raw = nullptr;
          target_state_.get_data(NODE, target_nodevar_names[i], &target_field_raw);
          Portage::pointer<double> target_field(target_field_raw);

          Portage::transform(target_mesh_.begin(NODE, PARALLEL_OWNED),
                             target_mesh_.end(NODE, PARALLEL_OWNED),
                             source_cells_and_weights.begin(),
                             target_field, interpolate);
          /*  UNCOMMENT WHEN WE RESTORE get_type in jali_state_wrapper
              } else {
              std::cerr << "Cannot remap " << source_var_names[i] <<
              " because it is not a scalar double variable\n";
              continue;
              }
          */
        }
      }

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_interp = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

      std::cout << "Transform Time Rank " << comm_rank << " (s): " <<
          tot_seconds << std::endl;
      std::cout << "  Search Time Rank " << comm_rank << " (s): " <<
          tot_seconds_srch << std::endl;
      std::cout << "  Intersect Time Rank " << comm_rank << " (s): " <<
          tot_seconds_xsect << std::endl;
      std::cout << "  Interpolate Time Rank " << comm_rank << " (s): " <<
          tot_seconds_interp << std::endl;
    }

  }


 private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  std::vector<LimiterType> limiters_;
  unsigned int dim_;
};  // class Driver


/*!
  @struct SearchFunctor "driver.h"
  @brief This functor is used inside a Portage::transform() inside
  Driver::run() to actually do the search
  @tparam SearchType The type of search method (e.g. SearchSimple or
  SearchKDTree).
  //amh: FIXME! this is a search adapter which converts a parameter which is a reference into a return by value
*/
template <typename SearchType>
struct SearchFunctor {
  const SearchType* search_;      ///< search method (e.g. SearchSimple)

  /*!
    @brief Constructor.
    @param[in] search The search method to use (e.g. SearchSimple)
  */
  SearchFunctor(const SearchType* search) : search_(search)
  {}

  /*!
    @brief Operator for making this struct a functor

    This is called from within a Portage::transform() operation that iterates
    over the cells in a target mesh.

    @param[in] targetCellindex The cell ID in the target mesh that this functor
    is currently operating on.

    @return list of candidate cells
  */

  std::vector<int> operator() (int const targetCellIndex) {

    std::vector<int> candidates;

    // Search for candidates and return their cells indices
    (*search_)(targetCellIndex, &candidates);

    return candidates;
  }
};


/*!
  @struct IntersectFunctor "driver.h"
  @brief This functor is used inside a Portage::transform() inside
  Driver::run() to actually do the intersection of the target cell with candidate source cells
  @tparam D         Dimension of problem
  @tparam IsectType The type of intersect method (e.g. IntersectClipper).
  @tparam SourceMeshWrapper   Portage compatible wrapper class for source mesh
  @tparam TargetMeshWrapper   Portage compatible wrapper class for target mesh
*/
template <int D,
          typename IsectType,
          typename SourceMeshWrapper,
          typename TargetMeshWrapper = SourceMeshWrapper>
struct IntersectFunctor {
  const IsectType* intersect_;    ///< intersect method (e.g. IntersectR2D)
  const SourceMeshWrapper & source_mesh_;
  const TargetMeshWrapper & target_mesh_;

  /*!
    @brief Constructor.
    @param[in] intersect The intersect method to use (e.g. IntersectR2D)
  */
  IntersectFunctor(const IsectType* intersect,
                   const SourceMeshWrapper & source_mesh,
                   const TargetMeshWrapper & target_mesh)
      : intersect_(intersect),
        source_mesh_(source_mesh), target_mesh_(target_mesh)
  {}

  /*!
    @brief Operator for making this struct a functor

    This is called from within a Portage::transform() operation that iterates
    over the cells in a target mesh.

    @param[in] targetCellindex   The cell ID in the target mesh that this functor is currently operating on.
    @param[in] candidates   Candidates to intersect with
    @return paired list of contributing cells and corresponding moment vectors
  */

  std::vector<Weights_t> operator() (int const targetCellIndex,
				     std::vector<int> const & candidates) {

    // Declaration fine for serial and OpenMP, won't work with CUDA
    Portage::vector<Portage::vector<double>> moments(candidates.size());

    // first extract a simplicial (tri/tet) decomposition of the
    // targetCell and the candidate cells
    std::vector<std::array<Portage::Point<D>, D+1>> tcoords;

    // vector of tri/tet coordinates for target cell arranged linearly
    // in the form required by the intersect functor
    std::vector<Portage::Point<D>> tcoords1;

    target_mesh_.decompose_cell_into_simplices(targetCellIndex, &tcoords);
    if (D == 2) {
      int ntris = tcoords.size();
      for (int i = 0; i < ntris; i++) {
        for (int j = 0; j < 3; j++)
          tcoords1.push_back(tcoords[i][j]);
      }
    } else if (D == 3) {
      int ntets = tcoords.size();
      for (int i = 0; i < ntets; i++) {
        for (int j = 0; j < 4; j++)
          tcoords1.push_back(tcoords[i][j]);
      }
    }

    // vector of vector of tri/tet coordinates for candidate source
    // cells arranged linearly in the form required by the intersect
    // functor
    std::vector<std::vector<Portage::Point<D>>> scoords1_vec(candidates.size());
    for (int s = 0; s < candidates.size(); s++) {
      int sourceCellIndex = candidates[s];
      std::vector<std::array<Point<D>, D+1>> scoords;
      source_mesh_.decompose_cell_into_simplices(sourceCellIndex, &scoords);
      if (D == 2) {
        int ntris = scoords.size();
        scoords1_vec[s].resize(4*ntris);
        for (int i = 0, k = 0; i < ntris; i++) {
          for (int j = 0; j < 3; j++)
            scoords1_vec[s][k++] = scoords[i][j];
        }
      } else if (D == 3) {
        int ntets = scoords.size();
        scoords1_vec[s].resize(4*ntets);
        for (int i = 0, k = 0; i < ntets; i++) {
          for (int j = 0; j < 4; j++)
            scoords1_vec[s][k++] = scoords[i][j];
        }
      }
    }

    // ONCE WE UPGRADE TO BOOST 1.61 or higher (which contains
    // constant_iterator) we can merge these two together (will need
    // to edit support/portage.h)

#ifdef THRUST
    // thrust::constant_iterator repeatedly returns a pointer to the
    // same object without having to make copies of the object

    thrust::transform(scoords1_vec.begin(), scoords1_vec.end(),
                      thrust::make_constant_iterator(tcoords1),
                      moments, *intersect_);
#else
    for (int i = 0; i < candidates.size(); i++)
      moments[i] = (*intersect_)(scoords1_vec[i], tcoords1);
#endif

    // Compute new value on target cell based on source mesh
    // values and intersection moments

    std::vector<Weights_t> source_cells_and_weights(candidates.size());

    for (int i = 0; i < candidates.size(); ++i) {
      Weights_t & this_wt = source_cells_and_weights[i];
      this_wt.entityID = candidates[i];
      this_wt.weights = moments[i];
    }

    return source_cells_and_weights;
  }
};  // struct IntersectFunctor

}  // namespace Portage

#endif  // SRC_DRIVER_DRIVER_H_
