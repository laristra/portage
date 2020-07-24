/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_CORE_DRIVER_H_
#define PORTAGE_CORE_DRIVER_H_

#include <ctime>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <memory>
#include <limits>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/CoordinateSystem.h"

#include "portage/support/portage.h"

#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/driver/driver.h"
#endif

#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/interpolate/gradient.h"
#include "portage/driver/parts.h"
#include "portage/driver/fix_mismatch.h"

/*!
  @file coredriver.h
  @brief Core driver for remapping

  Remap mesh and material variables from mesh to mesh for a specific
  entity kind in serial or on each partition (without redistribution)
*/

namespace Portage {

using Wonton::Entity_kind;
using Wonton::CELL;
using Wonton::NODE;
using Wonton::UNKNOWN_KIND;
using Wonton::Entity_type;
using Wonton::PARALLEL_OWNED;
using Wonton::ALL;
using Wonton::Point;

/*!
  @class CoreDriver "driver_core.h"

  @brief CoreDriver - Core driver that remaps fields on a particular Entity_kind (ONWHAT) like CELL or NODE

  NOTE: THIS CLASS ASSUMES THAT ALL SOURCE CELLS OVERLAPPING ANY TARGET
  CELL ARE AVAILABLE ON THIS PROCESSOR. IT DOES NOT HAVE TO FETCH THE
  DATA FROM ANYWHERE

  @tparam ONWHAT     On what kind of entity are we doing the remap

  @tparam SourceMesh A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.

  @tparam SourceState A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.

  @tparam TargetMesh A lightweight wrapper to a specific output mesh
  implementation that provides certain functionality.

  @tparam TargetState A lightweight wrapper to a specific output state
  manager implementation that provides certain functionality.

  @tparam InterfaceReconstructorType  The Interface Reconstructor class we will instantiate

  @tparam Matpoly_Splitter Class to split a polyhedron into two pieces

  @tparam Matpoly_Clipper Class to clip a polyhedron with a plane and
  return the piece behind the plane

  @tprarm CoordSys  Coordinate system being used for calculations

*/
template <int D,
          Entity_kind ONWHAT,
          class SourceMesh, class SourceState,
          class TargetMesh = SourceMesh, class TargetState = SourceState,
          template <class, int, class, class> class InterfaceReconstructorType = DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void,
          class CoordSys = Wonton::DefaultCoordSys
          >
class CoreDriver {

  // useful alias
  using Gradient = Limited_Gradient<D, ONWHAT, SourceMesh, SourceState,
                                    InterfaceReconstructorType,
                                    Matpoly_Splitter, Matpoly_Clipper, CoordSys>;

 public:
  /*!
    @brief Constructor for the CORE remap driver.

    @param[in] sourceMesh A @c wrapper to the source mesh (may be
    native or redistributed source).

    @param[in] sourceState A @c wrapper for the data that lives on the
    source mesh

    @param[in] targetMesh A @c TargetMesh to the target mesh

    @param[in,out] targetState A @c TargetState for the data that will
    be mapped to the target mesh
  */
  CoreDriver(SourceMesh const& source_mesh,
             SourceState const& source_state,
             TargetMesh const& target_mesh,
             TargetState& target_state,
             Wonton::Executor_type const *executor = nullptr)
      : source_mesh_(source_mesh),
        target_mesh_(target_mesh),
        source_state_(source_state),
        target_state_(target_state),
        gradient_(source_mesh, source_state),
        executor_(executor)
  {
#ifdef WONTON_ENABLE_MPI
    mycomm_ = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm_ = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm_, &comm_rank_);
      MPI_Comm_size(mycomm_, &nprocs_);
    }
#endif
  }

  /// Copy constructor (disabled)
  CoreDriver(const CoreDriver &) = delete;

  /// Assignment operator (disabled)
  CoreDriver & operator = (const CoreDriver &) = delete;

  /// Destructor
  ~CoreDriver() = default;

  /// What entity kind is this defined on?
  Entity_kind onwhat() {return ONWHAT;}

  /*!
    Find candidates entities of a particular kind that might
    intersect each target entity of the same kind

    @tparam Search Search class templated on dimension, Entity_kind
    and both meshes

    @return Vector of intersection candidates for each target entity
  */

  template<template<int, Entity_kind, class, class> class Search>
  Wonton::vector<std::vector<int>>
  search() {
    // Get an instance of the desired search algorithm type
    const Search<D, ONWHAT, SourceMesh, TargetMesh>
        search_functor(source_mesh_, target_mesh_);

    int ntarget_ents = target_mesh_.num_entities(ONWHAT, PARALLEL_OWNED);

    // initialize search candidate vector
    Wonton::vector<std::vector<int>> candidates(ntarget_ents);

    Wonton::transform(target_mesh_.begin(ONWHAT, PARALLEL_OWNED),
                       target_mesh_.end(ONWHAT, PARALLEL_OWNED),
                       candidates.begin(), search_functor);

    return candidates;
  }


  /*! 
    Intersect source and target mesh entities of kind
    'ONWHAT' and return the intersecting entities and moments of
    intersection for each entity

    @param candidates Vector of intersection candidates for each target entity

    @return vector of intersection moments for each target entity
  */

  template<template <Entity_kind, class, class, class,
                     template <class, int, class, class> class,
                     class, class> class Intersect>
  Wonton::vector<std::vector<Wonton::Weights_t>>
  intersect_meshes(Wonton::vector<std::vector<int>> const& candidates) {

#ifdef PORTAGE_HAS_TANGRAM
    // If user did NOT set tolerances for Tangram, use Portage tolerances
    if (reconstructor_tols_.empty()) {
      reconstructor_tols_ = { {1000, num_tols_.min_absolute_distance,
                                     num_tols_.min_absolute_volume},
                              {100, num_tols_.min_absolute_distance,
                                    num_tols_.min_absolute_distance} };
    }
    // If user set tolerances for Tangram, but not for Portage,
    // use Tangram tolerances
    else if (!num_tols_.user_tolerances) {
      num_tols_.min_absolute_distance = reconstructor_tols_[0].arg_eps;
      num_tols_.min_absolute_volume = reconstructor_tols_[0].fun_eps;
    }
#endif

    int nents = target_mesh_.num_entities(ONWHAT, PARALLEL_OWNED);
    Wonton::vector<std::vector<Portage::Weights_t>> sources_and_weights(nents);
      
    Intersect<ONWHAT, SourceMesh, SourceState, TargetMesh,
              InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper>
        intersector(source_mesh_, source_state_, target_mesh_, num_tols_);

    Wonton::transform(target_mesh_.begin(ONWHAT, PARALLEL_OWNED),
                       target_mesh_.end(ONWHAT, PARALLEL_OWNED),
                       candidates.begin(),
                       sources_and_weights.begin(),
                       intersector);

    return sources_and_weights;
  }

  /**
   * @brief Deduce weights for reverse remap by transposing
   *        the weight matrix used for forward remap.
   *
   * @param forward_weights: weights list for forward remap.
   * @return weights list for reverse remap.
   */
  Wonton::vector<std::vector<Wonton::Weights_t>> deduce_reverse_weights
    (Wonton::vector<std::vector<Wonton::Weights_t>> const& forward_weights) const {

    int const num_source_entities = source_mesh_.num_entities(ONWHAT, ALL);
    int const num_target_entities = target_mesh_.num_entities(ONWHAT, PARALLEL_OWNED);
    assert(unsigned(num_target_entities) == forward_weights.size());

    Wonton::vector<std::vector<Wonton::Weights_t>> reverse_weights(num_source_entities);

    Wonton::transform(source_mesh_.begin(ONWHAT, ALL),
                      source_mesh_.end(ONWHAT, ALL),
                      reverse_weights.begin(),
                      [&](int s) {
                        entity_weights_t entries;
                        entries.reserve(10);
                        for (int t = 0; t < num_target_entities; ++t) {
                          entity_weights_t const& list = forward_weights[t];
                          for (auto const& weight : list) {
                            if (weight.entityID == s) {
                              entries.emplace_back(t, weight.weights);
                              break;
                            }
                          }
                        }
                        return entries;
                      });

    return reverse_weights;
  }

  /// Set core numerical tolerances
  void set_num_tols(const double min_absolute_distance, 
                    const double min_absolute_volume) {
    num_tols_.min_absolute_distance = min_absolute_distance;
    num_tols_.min_absolute_volume = min_absolute_volume;
    num_tols_.user_tolerances = true;
  }

  /// Set all numerical tolerances
  void set_num_tols(const NumericTolerances_t& num_tols) {
    num_tols_ = num_tols;
  }

#ifdef PORTAGE_HAS_TANGRAM
  /*!
    @brief set options for interface reconstructor driver
    @param all_convex Should be set to false if the source mesh contains
    non-convex cells.
    @param tols The vector of tolerances for each moment during reconstruction.
    By default, the values are chosen based on tolerances specified for Portage
    in NumericTolerances_t struct. If both the tolerances for Portage and for
    Tangram are explicitly set by a user, they need to make sure that selected
    values are synced. If only the tolerances for Tangram are set by a user,
    then values in Portage's NumericTolerances_t are set based on the tols
    argument.
  */
  void set_interface_reconstructor_options(bool all_convex,
                                           const std::vector<Tangram::IterativeMethodTolerances_t> &tols =
                                             std::vector<Tangram::IterativeMethodTolerances_t>()) {
    reconstructor_tols_ = tols; 
    reconstructor_all_convex_ = all_convex; 
  }

#endif


  /*! 

    Intersect target mesh cells with source material polygons
    and return the intersecting entities and moments of intersection
    for each entity

    @param[in] candidates Vector of intersection candidates for each
    target entity

    @return Material-wise vector of intersection moments for each target entity

    NOTE: WE COULD SEND IN THE MESH-MESH INTERSECTION MOMENTS AND REUSE THEM
    IF A SOURCE CELL HAS ONLY ONE MATERIAL
  */

  template<
    template <Entity_kind, class, class, class,
              template <class, int, class, class> class,
              class, class> class Intersect
    >
  std::vector<Wonton::vector<std::vector<Weights_t>>>
  intersect_materials(Wonton::vector<std::vector<int>> const& candidates) {

#ifdef PORTAGE_HAS_TANGRAM
    // set numerical tolerance if not already done.
    // it avoids the necessity of calling 'intersect_meshes',
    // but should not change anything if it was already set.
    if (reconstructor_tols_.empty()) {
      // if the user did NOT set tolerances for Tangram, use Portage tolerances
      reconstructor_tols_ = {
        {1000, num_tols_.min_absolute_distance, num_tols_.min_absolute_volume},
        { 100, num_tols_.min_absolute_distance, num_tols_.min_absolute_distance}
      };
    } else if (not num_tols_.user_tolerances) {
      // if the user has set tolerances for Tangram, but not for Portage,
      // use Tangram tolerances
      num_tols_.min_absolute_distance = reconstructor_tols_[0].arg_eps;
      num_tols_.min_absolute_volume = reconstructor_tols_[0].fun_eps;
    }

    int nmats = source_state_.num_materials();
    // Make sure we have a valid interface reconstruction method instantiated

    assert(typeid(InterfaceReconstructorType<SourceMesh, D,
                  Matpoly_Splitter, Matpoly_Clipper >) !=
           typeid(DummyInterfaceReconstructor<SourceMesh, D,
                  Matpoly_Splitter, Matpoly_Clipper>));

    // Intel 18.0.1 does not recognize std::make_unique even with -std=c++14 flag *ugh*
    // interface_reconstructor_ =
    //     std::make_unique<Tangram::Driver<InterfaceReconstructorType, D,
    //                                      SourceMesh,
    //                                      Matpoly_Splitter,
    //                                      Matpoly_Clipper>
    //                      >(source_mesh_, tols, true);
    interface_reconstructor_ =
        std::unique_ptr<Tangram::Driver<InterfaceReconstructorType, D,
                                        SourceMesh,
                                        Matpoly_Splitter,
                                        Matpoly_Clipper>
                        >(new Tangram::Driver<InterfaceReconstructorType, D,
                          SourceMesh,
                          Matpoly_Splitter,
                          Matpoly_Clipper>(source_mesh_, reconstructor_tols_,
                                           reconstructor_all_convex_));

    int ntargetcells = target_mesh_.num_entities(CELL, PARALLEL_OWNED);


    std::vector<int> cell_num_mats, cell_mat_ids;
    std::vector<double> cell_mat_volfracs;
    std::vector<Wonton::Point<D>> cell_mat_centroids;

    // Extract volume fraction and centroid data for cells in compact
    // cell-centric form (ccc)

    ccc_vfcen_data(cell_num_mats, cell_mat_ids, cell_mat_volfracs,
                   cell_mat_centroids);

    interface_reconstructor_->set_volume_fractions(cell_num_mats,
                                                   cell_mat_ids,
                                                   cell_mat_volfracs,
                                                   cell_mat_centroids);
    interface_reconstructor_->reconstruct(executor_);

    // Make an intersector which knows about the source state (to be
    // able to query the number of materials, etc) and also knows
    // about the interface reconstructor so that it can retrieve pure
    // material polygons

    Intersect<CELL, SourceMesh, SourceState, TargetMesh,
              InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper>
        intersector(source_mesh_, source_state_, target_mesh_, num_tols_,
                    interface_reconstructor_);

    // Assume (with no harm for sizing purposes) that all materials
    // in source made it into target

    std::vector<Wonton::vector<std::vector<Weights_t>>>
        source_weights_by_mat(nmats);

    for (int m = 0; m < nmats; m++) {
      std::vector<int> matcellstgt;

      intersector.set_material(m);

      // For each cell in the target mesh get a list of
      // candidate-weight pairings (in a traditional mesh, not
      // particle mesh, the weights are moments). Note that this
      // candidate list is different from the search candidate list in
      // that it may not include some of the search candidates. Also,
      // note that for 2nd order and higher remaps, we get multiple
      // moments (0th, 1st, etc) for each target-source cell
      // intersection
      //
      // NOTE: IDEALLY WE WOULD REUSE THE MESH-MESH INTERSECTIONS
      // WHEN THE SOURCE CELL CONTAINS ONLY ONE MATERIAL
      //
      // UNFORTUNATELY, THE REQUIREMENT OF THE INTERSECT FUNCTOR IS
      // THAT IT CANNOT MODIFY STATE, THIS MEANS WE CANNOT STORE THE
      // MESH-MESH INTERSECTION VALUES AND REUSE THEM AS NECESSARY
      // FOR MESH-MATERIAL INTERSECTION COMPUTATIONS. WE COULD DO
      // SOME OTHER THINGS LIKE PROCESS ONLY TARGET CELLS THAT
      // POSSIBLY INTERSECT MULTI-MATERIAL SOURCE CELLS; TARGET
      // CELLS THAT ONLY INTERSECT SINGLE MATERIAL CELLS CAN REUSE
      // MESH-MESH INTERSECTIONS (still have to determine the
      // materials coming into the target cell though - a target
      // cell intersecting pure cells of different materials will
      // become mixed).


      Wonton::vector<std::vector<Weights_t>> this_mat_sources_and_wts(ntargetcells);
      Wonton::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                         target_mesh_.end(CELL, PARALLEL_OWNED),
                         candidates.begin(),
                         this_mat_sources_and_wts.begin(),
                         intersector);

      // LOOK AT INTERSECTION WEIGHTS TO DETERMINE WHICH TARGET CELLS
      // WILL GET NEW MATERIALS

      ntargetcells = target_mesh_.num_entities(CELL, PARALLEL_OWNED);

      for (int c = 0; c < ntargetcells; c++) {
        std::vector<Weights_t> const& cell_mat_sources_and_weights =
            this_mat_sources_and_wts[c];
        int nwts = cell_mat_sources_and_weights.size();
        double cell_mat_volume = 0.0;
        for (int s = 0; s < nwts; s++) {
          std::vector<double> const& wts = cell_mat_sources_and_weights[s].weights;
          cell_mat_volume += wts[0];
        }
        // Check that the volume of material we are adding to c is not miniscule
        if (cell_mat_volume > num_tols_.min_absolute_volume)
          matcellstgt.push_back(c);
      }

      // If any processor is adding this material to the target state,
      // add it on all the processors

      int nmatcells = matcellstgt.size();
      int nmatcells_global = nmatcells;
#ifdef WONTON_ENABLE_MPI
      if (mycomm_!= MPI_COMM_NULL)
        MPI_Allreduce(&nmatcells, &nmatcells_global, 1, MPI_INT, MPI_SUM,
                      mycomm_);
#endif

      if (nmatcells_global) {
        int nmatstrg = target_state_.num_materials();
        bool found = false;
        int m2 = -1;
        for (int i = 0; i < nmatstrg; i++)
          if (target_state_.material_name(i) == source_state_.material_name(m)) {
            found = true;
            m2 = i;
            break;
          }
        if (found) {  // material already present - just update its cell list
          target_state_.mat_add_cells(m2, matcellstgt);
        } else {
          // add material along with the cell list

          // NOTE: NOT ONLY DOES THIS ROUTINE ADD A MATERIAL AND ITS
          // CELLS TO THE STATEMANAGER, IT ALSO MAKES SPACE FOR
          // FIELD VALUES FOR THIS MATERIAL IN EVERY MULTI-MATERIAL
          // VECTOR IN THE STATE MANAGER. THIS ENSURES THAT WHEN WE
          // CALL mat_get_celldata FOR A MATERIAL IN MULTI-MATERIAL
          // STATE VECTOR IT WILL ALREADY HAVE SPACE ALLOCATED FOR
          // FIELD VALUES OF THAT MATERIAL. SOME STATE WRAPPERS
          // COULD CHOOSE TO MAKE THIS A SIMPLER ROUTINE THAT ONLY
          // STORES THE NAME AND THE CELLS IN THE MATERIAL AND
          // ACTUALLY ALLOCATE SPACE FOR FIELD VALUES OF A MATERIAL
          // IN A MULTI-MATERIAL FIELD WHEN mat_get_celldata IS
          // INVOKED.

          target_state_.add_material(source_state_.material_name(m),
                                     matcellstgt);
        }
      }
      else
        continue;  // maybe the target mesh does not overlap this material

      // Add volume fractions and centroids of materials to target mesh
      //
      // Also make list of sources/weights only for target cells that are
      // getting this material - Can we avoid the copy?

      std::vector<double> mat_volfracs(nmatcells);
      std::vector<Point<D>> mat_centroids(nmatcells);

      source_weights_by_mat[m].resize(nmatcells);

      for (int ic = 0; ic < nmatcells; ic++) {
        int c = matcellstgt[ic];
        double matvol = 0.0;
        Point<D> matcen;
        std::vector<Weights_t> const &
            cell_mat_sources_and_weights = this_mat_sources_and_wts[c];
        int nwts = cell_mat_sources_and_weights.size();
        for (int s = 0; s < nwts; s++) {
          std::vector<double> const& wts = cell_mat_sources_and_weights[s].weights;
          matvol += wts[0];
          for (int d = 0; d < D; d++)
            matcen[d] += wts[d+1];
        }
        matcen /= matvol;
        mat_volfracs[ic] = matvol/target_mesh_.cell_volume(c);
        mat_centroids[ic] = matcen;

        source_weights_by_mat[m][ic] = cell_mat_sources_and_weights;
      }

      target_state_.mat_add_celldata("mat_volfracs", m, &(mat_volfracs[0]));
      target_state_.mat_add_celldata("mat_centroids", m, &(mat_centroids[0]));

    }  // for each material m

    return source_weights_by_mat;
#else
    return std::vector<Wonton::vector<std::vector<Weights_t>>>();
#endif

  }

#ifdef PORTAGE_HAS_TANGRAM
  /**
   * @brief Cache gradient stencil matrices for multi-material fields.
   *        It cannot be invoked directly in the constructor since
   *        the interface reconstructor is only initialized after
   *        the intersection step.
   */
  void cache_multimat_gradient_stencils() {
    // make sure that the interface reconstructor is initialized.
    // it is only done after the intersection step
    // this check enforces that this method is called only after that step.
    if (interface_reconstructor_) {
      if (not cached_multimat_stenc_) {
        gradient_.set_interface_reconstructor(interface_reconstructor_);
        gradient_.cache_matrices(Field_type::MULTIMATERIAL_FIELD);
        cached_multimat_stenc_ = true;
      }
    } else
      throw std::runtime_error("interface reconstructor not yet initialized");
  }

  /**
   * @brief Deduce reverse material weights by transposing
   *        the weight matrix used in forward remap for each
   *        material.
   *
   * Here we assume that the source state has already the material
   * fields and material polytopes moments (volume and centroids).
   *
   * @param forward_weights: weights list per material for forward remap.
   * @return the weights list per material for reverse remap.
   */
  std::vector<Wonton::vector<std::vector<Weights_t>>> deduce_reverse_material_weights
    (std::vector<Wonton::vector<std::vector<Wonton::Weights_t>>> const& forward_weights,
     std::vector<std::vector<int>>& index_mapping) {

    int const num_materials = forward_weights.size();
    int const num_source_cells = source_mesh_.num_entities(ONWHAT, ALL);

    std::vector<Wonton::vector<std::vector<Weights_t>>> reverse_weights(num_materials);

    for (int m = 0; m < num_materials; ++m) {
      int const num_target_material_cells = forward_weights[m].size();
      Wonton::vector<std::vector<Weights_t>> reverse_material_weights(num_source_cells);

      Wonton::transform(source_mesh_.begin(ONWHAT, ALL),
                        source_mesh_.end(ONWHAT, ALL),
                        reverse_material_weights.begin(),
                        [&](int s) {
                          entity_weights_t entries;
                          entries.reserve(10);
                          for (int t = 0; t < num_target_material_cells; ++t) {
                            entity_weights_t const& list = forward_weights[m][t];
                            for (auto const& weight : list) {
                              if (weight.entityID == s) {
                                entries.emplace_back(t, weight.weights);
                                break;
                              }
                            }
                          }
                          return entries;
                        });

      // filter weight list
      for (int s = 0; s < num_source_cells; ++s) {
        entity_weights_t const& weights = reverse_material_weights[s];
        if (not weights.empty()) {
          index_mapping[m].emplace_back(s);
        }
      }

      reverse_weights[m].resize(index_mapping[m].size());
      Wonton::transform(index_mapping[m].begin(),
                        index_mapping[m].end(),
                        reverse_weights[m].begin(),
                        [&](int s) { return reverse_material_weights[s]; });
    }

    return reverse_weights;
  }
#endif

  /**
   * @brief Compute the gradient field of the given variable on source mesh.
   *
   * @param field_name: the variable name.
   * @param limiter_type: gradient limiter to use on internal regions.
   * @param boundary_limiter_type: gradient limiter to use on boundary.
   * @param source_part: the source mesh part to consider if any.
   *
   * Remark: the gradient computation cannot be done in parallel yet
   * for multiple fields due to some side effects. Indeed the same
   * instance is used for multiple fields but limiter options
   * are specified at runtime for each field. Besides, if the stencil
   * matrices used for the least square approximation are not yet cached
   * then it will be done at this step.
   */
  Wonton::vector<Vector<D>> compute_source_gradient(
    std::string const field_name,
    Limiter_type limiter_type = NOLIMITER,
    Boundary_Limiter_type boundary_limiter_type = BND_NOLIMITER,
    int material_id = 0,
    const Part<SourceMesh, SourceState>* source_part = nullptr) {

    int nallent = 0;
#ifdef PORTAGE_HAS_TANGRAM
    // enable part-by-part only for cell-based remap
    auto const field_type = source_state_.field_type(ONWHAT, field_name);

    // multi-material remap makes only sense on cell-centered fields.
    bool const multimat =
      ONWHAT == Entity_kind::CELL and
      field_type == Field_type::MULTIMATERIAL_FIELD;

    std::vector<int> mat_cells;

    if (multimat) {
      // cache gradient stencil first
      if (not cached_multimat_stenc_) {
        std::cerr << "Warning: gradient stencil matrices for ";
        std::cerr << "multi-material fields were not cached yet." << std::endl;
        std::cerr << "Please invoke 'cache_multimat_gradient_stencils' ";
        std::cerr << "prior to 'compute_source_gradient' for optimized runs.";
        std::cerr << std::endl;
        cache_multimat_gradient_stencils();
      }

      if (interface_reconstructor_) {
        std::vector<int> mat_cells_all;
        source_state_.mat_get_cells(material_id, &mat_cells_all);
        nallent = mat_cells_all.size();

        // Filter out GHOST cells
        // SHOULD BE IN HANDLED IN THE STATE MANAGER (See ticket LNK-1589)
        mat_cells.reserve(nallent);
        for (auto const& c : mat_cells_all)
          if (source_mesh_.cell_get_type(c) == PARALLEL_OWNED)
            mat_cells.push_back(c);
      }
      else
        throw std::runtime_error("interface reconstructor not set");
    } else /* single material */ {
#endif
      nallent = source_mesh_.num_entities(ONWHAT, ALL);
#ifdef PORTAGE_HAS_TANGRAM
    }
#endif

    // use stored instance for mesh remap
    // create a new instance for part-by-part
    // nb: make_unique or make_shared would copy the object
    auto kernel = (source_part == nullptr ? &gradient_
                                          : new Gradient(source_mesh_, source_state_, source_part));
    // set gradient kernel options
    kernel->set_interpolation_variable(field_name, limiter_type, boundary_limiter_type);

    // create the field (material cell indices have owned and ghost
    // cells mixed together; so we have to have a vector of size
    // owned+ghost and fill in the right entries; the ghost entries
    // are zeroed out)
    Vector<D> zerovec;
    Wonton::vector<Vector<D>> gradient_field(nallent, zerovec);

    // populate it by invoking the kernel on each source entity.
#ifdef PORTAGE_HAS_TANGRAM
    if (multimat) {
      // no need for this to be Wonton::vector as it will be copied out
      std::vector<Vector<D>> owned_gradient_field(mat_cells.size());

      kernel->set_material(material_id);
      kernel->set_interface_reconstructor(interface_reconstructor_);
      Wonton::transform(mat_cells.begin(),
                         mat_cells.end(),
                         owned_gradient_field.begin(), *kernel);
      int i = 0;
      for (auto const& c : mat_cells) {
        int cm = source_state_.cell_index_in_material(c, material_id);
        gradient_field[cm] = owned_gradient_field[i++];
      }
    } else {
#endif
      Wonton::transform(source_mesh_.begin(ONWHAT, PARALLEL_OWNED),
                         source_mesh_.end(ONWHAT, PARALLEL_OWNED),
                         gradient_field.begin(), *kernel);
#ifdef PORTAGE_HAS_TANGRAM
    }
#endif
    if (source_part != nullptr) { delete kernel; }

    return gradient_field;
  }

  /**
   * @brief Interpolate mesh variable.
   *
   * @param[in] srcvarname          source mesh variable to remap
   * @param[in] trgvarname          target mesh variable to remap
   * @param[in] sources_and_weights weights for mesh-mesh interpolation
   * @param[in] gradients           gradients of variable on source mesh (can be nullptr for 1st order remap)
   */
  template<typename T = double,
           template<int, Entity_kind, class, class, class, class, class,
    template<class, int, class, class> class,
    class, class, class> class Interpolate
  >
  void interpolate_mesh_var(std::string srcvarname, std::string trgvarname,
                            Wonton::vector<std::vector<Weights_t>> const& sources_and_weights,
                            Wonton::vector<Vector<D>>* gradients = nullptr) {

    if (source_state_.get_entity(srcvarname) != ONWHAT) {
      std::cerr << "Variable " << srcvarname << " not defined on Entity_kind "
                << ONWHAT << ". Skipping!" << std::endl;
      return;
    }


    using Interpolator = Interpolate<D, ONWHAT,
                                     SourceMesh, TargetMesh,
                                     SourceState, TargetState,
                                     T,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper, CoordSys>;

    Interpolator interpolator(source_mesh_, target_mesh_, source_state_,
                              num_tols_);
    interpolator.set_interpolation_variable(srcvarname, gradients);

    // get a handle to a memory location where the target state
    // would like us to write this material variable into.
    T* target_mesh_field = nullptr;
    target_state_.mesh_get_data(ONWHAT, trgvarname, &target_mesh_field);

    Wonton::pointer<T> target_field(target_mesh_field);
    Wonton::transform(target_mesh_.begin(ONWHAT, PARALLEL_OWNED),
                       target_mesh_.end(ONWHAT, PARALLEL_OWNED),
                       sources_and_weights.begin(),
                       target_field, interpolator);
  }


  /**
   * @brief Interpolate mesh variable from source part to target part
   *
   * @param[in] srcvarname          source mesh variable to remap
   * @param[in] trgvarname          target mesh variable to remap
   * @param[in] sources_and_weights weights for mesh-mesh interpolation
   * @param[in] lower_bound         lower bound of variable value 
   * @param[in] upper_bound         upper bound of variable value 
   * @param[in] partition           structure containing source and target part
   * @param[in] gradients           gradients of variable on source mesh (can be nullptr for 1st order remap)

   Enable only for cells using SFINAE. Here the class, rather than the
   function is templated on ONWHAT (as opposed to the equivalent
   method in the base class); so we have to create a dummy template
   parameter ONWHAT1 and rely on that to use SFINAE with a _second_
   dummy template parameter

   **** Note ****
   If you encounter errors about not being able to find an appropriate
   overload for interpolate_mesh_var in your application code
   (particularly something like "no type named 'type' in struct
   std::enable_if<false, void>"), make sure the compiler does see the
   possiblity of calling this function with Entity_kinds that are not
   type CELL (restricting the code flow using 'if' statements will not
   be enough)
   */
  template<typename T = double,
           template<int, Entity_kind, class, class, class, class, class,
                    template<class, int, class, class> class,
                    class, class, class> class Interpolate,
           Entity_kind ONWHAT1 = ONWHAT,
           typename = typename std::enable_if<ONWHAT1 == CELL>::type>
  void
  interpolate_mesh_var(std::string srcvarname, std::string trgvarname,
                       Wonton::vector<std::vector<Weights_t>> const& sources_and_weights,
                       const PartPair<D, SourceMesh, SourceState,
                       TargetMesh, TargetState>* partition,
                       Wonton::vector<Vector<D>>* gradients = nullptr) {

    if (source_state_.get_entity(srcvarname) != ONWHAT) {
      std::cerr << "Variable " << srcvarname << " not defined on Entity_kind "
                << ONWHAT << ". Skipping!" << std::endl;
      return;
    }


    using Interpolator = Interpolate<D, ONWHAT,
                                     SourceMesh, TargetMesh,
                                     SourceState, TargetState,
                                     T,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper, CoordSys>;

    Interpolator interpolator(source_mesh_, target_mesh_, source_state_, num_tols_, partition);
    interpolator.set_interpolation_variable(srcvarname, gradients);

    // get a handle to a memory location where the target state
    // would like us to write this material variable into.
    T* target_mesh_field = nullptr;
    target_state_.mesh_get_data(ONWHAT, trgvarname, &target_mesh_field);

    // perform part-by-part interpolation
    assert(ONWHAT == Entity_kind::CELL && partition != nullptr);

    // 1. Do some basic checks on supplied source and target parts
    // to prevent bugs when interpolating values:
    // check that each entity id is within the
    // mesh entity index space.
    auto const& source_part = partition->source();
    auto const& target_part = partition->target();

#ifndef NDEBUG
    int const& max_source_id = source_mesh_.num_entities(ONWHAT, ALL);
    int const& max_target_id = target_mesh_.num_entities(ONWHAT, ALL);

    Wonton::for_each(source_part.cells().begin(),
                      source_part.cells().end(),
                      [&](int current){ assert(current <= max_source_id); });

    Wonton::for_each(target_part.cells().begin(),
                      target_part.cells().end(),
                      [&](int current){ assert(current <= max_target_id); });
#endif

    int const target_part_size = target_part.size();

    // 2. Filter intersection weights list.
    // To restrict interpolation only to source-target parts, we need
    // to filter the intersection weights list to keep only that of the
    // entities of the source part. Notice that this step can be avoided
    // when the part-by-part intersection is implemented.

    auto filter_weights = [&](int entity) {
      // For a given target entity, we aim to filter its source weights
      // list to keep only those which are in the source part list.
      // that way, we ensure that only the contribution of source part entities
      // weights are taken into account when doing the interpolation.
      // For that, we just iterate on the related weight list, and add the
      // current couple of entity/weights if it belongs to the source part.
      // nb: 'auto' may imply unexpected behavior with thrust enabled.
      entity_weights_t const& entity_weights = sources_and_weights[entity];
      entity_weights_t heap;
      heap.reserve(10); // size of a local vicinity
      for (auto&& weight : entity_weights) {
        // constant-time lookup in average case.
        if(source_part.contains(weight.entityID)) {
          heap.emplace_back(weight);
        }
      }
      heap.shrink_to_fit();
      return heap;
    };

    Wonton::vector<entity_weights_t> parts_weights(target_part_size);
    Wonton::transform(target_part.cells().begin(),
                       target_part.cells().end(),
                       parts_weights.begin(), filter_weights);

    // 3. Process interpolation.
    // Now that intersection weights is filtered, perform the interpolation.
    // Notice that we need to store the interpolated values in a temporary
    // array since they are indexed with respect to the target part
    // and not to the target mesh. Hence we need to copy them back in the
    // target state at their correct (absolute index) memory locations.
    T temporary_storage[target_part_size];
    Wonton::pointer<T> target_part_field(temporary_storage);

    Wonton::transform(target_part.cells().begin(),
                       target_part.cells().end(),
                       parts_weights.begin(), target_part_field, interpolator);

    for (int i=0; i < target_part_size; ++i) {
      auto const& j = target_part.cells()[i];
      target_mesh_field[j] = target_part_field[i];
    }
  }
  

#ifdef PORTAGE_HAS_TANGRAM
  
  /*! CoreDriver::interpolate_mat_var

    @brief interpolate a material variable

    @param[in] srcvarname  Material variable name on the source mesh

    @param[in] trgvarname  Material variable name on the target mesh

    @param[in] bnd_limiter Boundary limiter to use for variable

    @param[in] lower_bound Lower bound of variable value

    @param[in] upper_bound Upper bound of variable value

    Enable only for cells using SFINAE. Here the class, rather than
    the function is templated on ONWHAT; so we have to create a dummy
    template parameter ONWHAT1 and rely on that for SFINAE using a
    _second_ template parameter

   **** Note ****
   If you encounter errors about not being able to find an appropriate
   overload for interpolate_mesh_var in your application code
   (particularly something like "no type named 'type' in struct
   std::enable_if<false, void>"), make sure the compiler does see the
   possiblity of calling this function with Entity_kinds that are not
   type CELL (restricting the code flow using 'if' statements will not
   be enough)
  */

  template<typename T = double,
           template<int, Entity_kind, class, class, class, class, class,
                    template<class, int, class, class> class,
                    class, class, class> class Interpolate,
           Entity_kind ONWHAT1 = ONWHAT,
           typename = typename std::enable_if<ONWHAT1 == CELL>::type>
  void
  interpolate_mat_var(std::string srcvarname, std::string trgvarname,
                      std::vector<Wonton::vector<std::vector<Weights_t>>> const& sources_and_weights_by_mat,
                      std::vector<Wonton::vector<Vector<D>>>* gradients = nullptr) {
    
    using Interpolator = Interpolate<D, ONWHAT,
                                     SourceMesh, TargetMesh,
                                     SourceState, TargetState,
                                     T,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper, CoordSys>;

    Interpolator interpolator(source_mesh_, target_mesh_,
                              source_state_, num_tols_,
                              interface_reconstructor_);
      
    int const nmats = source_state_.num_materials();

    for (int m = 0; m < nmats; m++) {

      interpolator.set_material(m);    // We have to do this so we know
      //                               // which material values we have
      //                               // to grab from the source state

      auto mat_grad = (gradients != nullptr ? &((*gradients)[m]) : nullptr);
      // FEATURE ;-)  Have to set interpolation variable AFTER setting 
      // the material for multimaterial variables
      interpolator.set_interpolation_variable(srcvarname, mat_grad);

      // if the material has no cells on this partition, then don't bother
      // interpolating MM variables
      if (target_state_.mat_get_num_cells(m) == 0) continue;

      std::vector<int> matcellstgt;
      target_state_.mat_get_cells(m, &matcellstgt);

      // Get a handle to a memory location where the target state
      // would like us to write this material variable into. If it is
      // NULL, we allocate it ourself

      T *target_field_raw;
      target_state_.mat_get_celldata(trgvarname, m, &target_field_raw);
      assert (target_field_raw != nullptr);

      Wonton::pointer<T> target_field(target_field_raw);

      Wonton::transform(matcellstgt.begin(), matcellstgt.end(),
                         sources_and_weights_by_mat[m].begin(),
                         target_field, interpolator);

      // If the state wrapper knows that the target data is already
      // laid out in this way and it gave us a pointer to the array
      // where the values reside, it has to do nothing in this
      // call. If the storage format is different, however, it may
      // have to copy the values into their proper locations

      target_state_.mat_add_celldata(trgvarname, m, target_field_raw);
    }  // over all mats

  }  // CoreDriver::interpolate_mat_var

#endif  // PORTAGE_HAS_TANGRAM


  /*! 
    Check mismatch between meshes

    @param[in] sources_and_weights Intersection sources and moments
    (vols, centroids)

    @returns   Whether the meshes are mismatched
  */
  bool
  check_mismatch(Wonton::vector<std::vector<Weights_t>> const& source_weights) {

    // Instantiate mismatch fixer for later use
    if (not mismatch_fixer_) {
      // Intel 18.0.1 does not recognize std::make_unique even with -std=c++14 flag *ugh*
      // mismatch_fixer_ = std::make_unique<MismatchFixer<D, ONWHAT,
      //                                                  SourceMesh, SourceState,
      //                                                  TargetMesh,  TargetState>
      //                                    >
      //     (source_mesh_, source_state_, target_mesh_, target_state_,
      //      source_weights, executor_);

      mismatch_fixer_ = std::unique_ptr<MismatchFixer<D, ONWHAT,
                                                      SourceMesh, SourceState,
                                                      TargetMesh,  TargetState>
                                        >(new MismatchFixer<D, ONWHAT,
                                          SourceMesh, SourceState,
                                          TargetMesh,  TargetState>
                                          (source_mesh_, source_state_, target_mesh_, target_state_,
                                           executor_));
    }

    return mismatch_fixer_->check_mismatch(source_weights);
  }

  
  /*! 
    Return mismatch between meshes

    @returns   Whether the meshes are mismatched
  */
  bool
  has_mismatch() {
    assert(mismatch_fixer_ && "check_mismatch must be called first!");
    return mismatch_fixer_->has_mismatch();
  }

  /// @brief Repair the remapped field to account for boundary mismatch
  /// @param src_var_name        field variable on source mesh
  /// @param trg_var_name        field variable on target mesh
  /// @param global_lower_bound  lower limit on variable
  /// @param global_upper_bound  upper limit on variable
  /// @param partial_fixup_type  type of fixup in case of partial mismatch
  /// @param empty_fixup_type    type of fixup in empty target entities
  ///
  /// partial_fixup_type can be one of three types:
  ///
  /// CONSTANT - Fields will see no perturbations BUT REMAP WILL BE
  ///            NON-CONSERVATIVE (constant preserving, not linearity
  ///            preserving)
  /// LOCALLY_CONSERVATIVE - REMAP WILL BE LOCALLY CONSERVATIVE (target cells
  ///                        will preserve the integral quantities received from
  ///                        source mesh overlap) but perturbations will
  ///                        occur in the field (constant fields may not stay
  ///                        constant if there is mismatch)
  /// SHIFTED_CONSERVATIVE - REMAP WILL BE CONSERVATIVE and field
  ///                        perturbations will be minimum but field
  ///                        values may be shifted (Constant fields
  ///                        will be shifted to different constant; no
  ///                        guarantees on linearity preservation)
  ///
  /// empty_fixup_type can be one of two types:
  ///
  /// LEAVE_EMPTY - Leave empty cells as is
  /// EXTRAPOLATE - Fill empty cells with extrapolated values
  /// FILL        - Fill empty cells with specified values (not yet implemented)
  bool fix_mismatch(std::string const & src_var_name,
                    std::string const & trg_var_name,
                    double global_lower_bound = -std::numeric_limits<double>::max(),
                    double global_upper_bound = std::numeric_limits<double>::max(),
                    double conservation_tol = 1e2*std::numeric_limits<double>::epsilon(),
                    int maxiter = 5,
                    Partial_fixup_type partial_fixup_type =
                    Partial_fixup_type::SHIFTED_CONSERVATIVE,
                    Empty_fixup_type empty_fixup_type =
                    Empty_fixup_type::EXTRAPOLATE) {

    assert(mismatch_fixer_ && "check_mismatch must be called first!");

    if (source_state_.field_type(ONWHAT, src_var_name) ==
        Field_type::MESH_FIELD)
      return mismatch_fixer_->fix_mismatch(src_var_name, trg_var_name,
                                  global_lower_bound, global_upper_bound,
                                  conservation_tol, maxiter,
                                  partial_fixup_type, empty_fixup_type);
    return false;
  }
  
 private:
  SourceMesh const & source_mesh_;
  TargetMesh const & target_mesh_;
  SourceState const & source_state_;
  TargetState & target_state_;
  Gradient gradient_;
  NumericTolerances_t num_tols_ = DEFAULT_NUMERIC_TOLERANCES<D>;

  int comm_rank_ = 0;
  int nprocs_ = 1;

  Wonton::Executor_type const *executor_;

#ifdef WONTON_ENABLE_MPI
  MPI_Comm mycomm_ = MPI_COMM_NULL;
#endif

#ifdef PORTAGE_HAS_TANGRAM
  bool cached_multimat_stenc_ = false;

  // The following tolerances as well as the all-convex flag are
  // required for the interface reconstructor driver. The size of the
  // tols vector is currently set to two since MOF requires two
  // different set of tolerances to match the 0th-order and 1st-order
  // moments. VOF on the other does not require the second tolerance.
  // If a new IR method which requires tolerances for higher moment is
  // added to Tangram, then this vector size should be
  // generalized. The boolean all_convex flag is to specify if a mesh
  // contains only convex cells and set to true in that case.
  //
  // There is an associated method called
  // set_interface_reconstructor_options that should be invoked to set
  // user-specific values. Otherwise, the remapper will use the
  // default values.

  std::vector<Tangram::IterativeMethodTolerances_t> reconstructor_tols_;
  bool reconstructor_all_convex_ = true;  

  // Pointer to the interface reconstructor object (required by the
  // interface to be shared)
  std::shared_ptr<Tangram::Driver<InterfaceReconstructorType, D,
                                  SourceMesh,
                                  Matpoly_Splitter, Matpoly_Clipper>
                  > interface_reconstructor_;
  

  // Convert volume fraction and centroid data from compact
  // material-centric to compact cell-centric (ccc) form as needed
  // by Tangram
  void ccc_vfcen_data(std::vector<int>& cell_num_mats,
                      std::vector<int>& cell_mat_ids,
                      std::vector<double>& cell_mat_volfracs,
                      std::vector<Wonton::Point<D>>& cell_mat_centroids) {

    int nsourcecells = source_mesh_.num_entities(CELL,
                                                 ALL);

    int nmats = source_state_.num_materials();
    cell_num_mats.assign(nsourcecells, 0);

    // First build full arrays (as if every cell had every material)

    std::vector<int> cell_mat_ids_full(nsourcecells*nmats, -1);
    std::vector<double> cell_mat_volfracs_full(nsourcecells*nmats, 0.0);
    std::vector<Wonton::Point<D>> cell_mat_centroids_full(nsourcecells*nmats);

    bool have_centroids = true;
    int nvals = 0;
    for (int m = 0; m < nmats; m++) {
      std::vector<int> cellids;
      source_state_.mat_get_cells(m, &cellids);
      for (int c : cellids) {
        int nmatc = cell_num_mats[c];
        cell_mat_ids_full[c*nmats+nmatc] = m;
        cell_num_mats[c]++;
      }
      nvals += cellids.size();

      double const * matfracptr;
      source_state_.mat_get_celldata("mat_volfracs", m, &matfracptr);
      int const num_cell_ids = cellids.size();
      for (int ic = 0; ic < num_cell_ids; ic++)
        cell_mat_volfracs_full[cellids[ic]*nmats+m] = matfracptr[ic];

      Wonton::Point<D> const *matcenvec;
      
      // handle the case where we don't have centroids in the state manager at all
      if (source_state_.get_entity("mat_centroids")==Entity_kind::UNKNOWN_KIND) {
        have_centroids = false;
        continue;
      }	
      
      source_state_.mat_get_celldata("mat_centroids", m, &matcenvec);
      if (cellids.size() && !matcenvec) {
        have_centroids = false;  // VOF
      } else {
        for (int ic = 0; ic < num_cell_ids; ic++)
          cell_mat_centroids_full[cellids[ic]*nmats+m] = matcenvec[ic];
      }
    }

    // At this point nvals contains the number of non-zero volume
    // fraction entries in the full array. Use this and knowledge of
    // number of materials in each cell to compress the data into
    // linear arrays

    cell_mat_ids.resize(nvals);
    cell_mat_volfracs.resize(nvals);
    cell_mat_centroids.resize(nvals);  // dummy vals for VOF

    int idx = 0;
    for (int c = 0; c < nsourcecells; c++) {
      for (int m = 0; m < cell_num_mats[c]; m++) {
        int matid = cell_mat_ids_full[c*nmats+m];
        cell_mat_ids[idx] = matid;
        cell_mat_volfracs[idx] = cell_mat_volfracs_full[c*nmats+matid];
        if (have_centroids)
          cell_mat_centroids[idx] = cell_mat_centroids_full[c*nmats+matid];
        idx++;
      }
    }
  }

#endif

  std::unique_ptr<MismatchFixer<D, ONWHAT,
                                SourceMesh, SourceState,
                                TargetMesh, TargetState>> mismatch_fixer_;

};  // CoreDriver

}  // namespace Portage

#endif  // PORTAGE_CORE_DRIVER_H_
