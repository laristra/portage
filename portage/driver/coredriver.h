/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_CORE_DRIVER_H_
#define PORTAGE_CORE_DRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <memory>
#include <limits>

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"

#include "portage/intersect/dummy_interface_reconstructor.h"
#endif

#include "portage/support/portage.h"

#include "wonton/support/Point.h"
#include "wonton/state/state_vector_multi.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/driver/parts.h"

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

// Forward declaration of remap driver on particular entity kind

template <int D,
          Entity_kind ONWHAT,
          class SourceMesh, class SourceState,
          class TargetMesh, class TargetState,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper,
          class CoordSys
          >
class CoreDriver;


/*!
  @class CoreDriverBase "coredriver.h"

  @brief CoreDriverBase - Base class for core driver that is
  agnostic to the Entity_kind

  @tparam SourceMesh A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.

  @tparam SourceState A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.

  @tparam TargetMesh A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.

  @tparam TargetState A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.

  @tparam InterfaceReconstructorType An interface reconstruction class
  that takes the raw interface reconstruction method, the dimension of
  the problem and the source mesh class as template parameters

  @tparam Matpoly_Splitter A polygon/polyhedron splitting class (returns both pieces of the polygon)

  @tparam Matpoly_Clipper A polygon/polyhedron clipping class (returns only the piece below/behind the clipping plane)

  @tparam CoordinateSystem 

*/

template<int D,
         class SourceMesh, class SourceState,
         class TargetMesh, class TargetState,
         template <class, int, class, class> class InterfaceReconstructorType,
         class Matpoly_Splitter,
         class Matpoly_Clipper,
         class CoordSys
         >
class CoreDriverBase {

  template<Entity_kind ONWHAT>
  using CoreDriverType = CoreDriver<D, ONWHAT, SourceMesh, SourceState,
                                    TargetMesh, TargetState,
                                    InterfaceReconstructorType,
                                    Matpoly_Splitter, Matpoly_Clipper,
                                    CoordSys>;
  
 public:
  CoreDriverBase() {}
  virtual ~CoreDriverBase() {}  // Necessary
  
  // Entity kind that a derived class is defined on
  virtual Entity_kind onwhat() = 0;

  /*! @brief search for candidate source entities whose control volumes
    (cells, dual cells) overlap the control volumes of target cells
     
    @tparam Entity_kind  what kind of entity are we searching on/for

    @tparam Search       search functor

    @returns    vector of intersection candidates for each target entity
  */

  template<Entity_kind ONWHAT,
           template <int, Entity_kind, class, class> class Search>
  Portage::vector<std::vector<int>>
  search() {
    assert(ONWHAT == onwhat());
    auto derived_class_ptr = static_cast<CoreDriverType<ONWHAT> *>(this);
    return derived_class_ptr->template search<Search>();
  }
    

  /*! @brief intersect target entities with candidate source entities

    @tparam Entity_kind  what kind of entity are we searching on/for

    @tparam Intersect    intersect functor

    @param[in] intersection_candidates  vector of intersection candidates for each target entity

    @returns  vector of intersection moments for each target entity
  */

  template<
    Entity_kind ONWHAT,
    template <Entity_kind, class, class, class,
              template <class, int, class, class> class,
              class, class> class Intersect
    >
  Portage::vector<std::vector<Portage::Weights_t>>
  intersect_meshes(Portage::vector<std::vector<int>> const& intersection_candidates) {
    assert(ONWHAT == onwhat());
    auto derived_class_ptr = static_cast<CoreDriverType<ONWHAT> *>(this);
    return derived_class_ptr->template intersect_meshes<Intersect>(intersection_candidates);
  }
    



  /*! intersect target cells with source material polygons

    @tparam Intersect   intersect functor

    @param[in] intersection_candidates vector of intersection
    candidates for each target entity

    @returns material-wise vector of intersection moments for each
    target cell
  */

  template<
    template <Entity_kind, class, class, class,
              template <class, int, class, class> class,
              class, class> class Intersect
    >
  std::vector<Portage::vector<std::vector<Portage::Weights_t>>>
  intersect_materials(Portage::vector<std::vector<int>> const& intersection_candidates) {
    assert(onwhat() == CELL);
    auto derived_class_ptr = static_cast<CoreDriverType<CELL> *>(this);
    return derived_class_ptr->template intersect_materials<Intersect>(intersection_candidates);
  }


  /*!
    Interpolate a mesh variable of type T residing on entity kind ONWHAT using
    previously computed intersection weights

    @tparam T   type of variable

    @tparam ONWHAT  Entity_kind that field resides on

    @tparam Interpolate  Functor for doing the interpolate from mesh to mesh

    @param[in] srcvarname   Variable name on source mesh

    @param[in] trgvarname   Variable name on target mesh

    @param[in] lower_bound  Lower bound for variable

    @param[in] upper_bound  Upper bound for variable

    @param[in] limiter      Limiter to use for second order reconstruction

    @param[in] partial_fixup_type Method to populate fields on
    partially filled target entities (cells or dual cells)

    @param[in] empty_fixup_type Method to populate fields on empty
    target entities (cells or dual cells)

    @param[in] conservation_tol Tolerance to which source and target
    integral quantities are to be matched

    @param[in] max_fixup_iter     Max number of iterations for global repair

    See support/portage.h for options on limiter, partial_fixup_type and
    empty_fixup_type
    
  */
  
  template<typename T = double,
           Entity_kind ONWHAT,
           template<int, Entity_kind, class, class, class,
                    template <class, int, class, class> class,
                    class, class, class> class Interpolate
           >
  void interpolate_mesh_var(std::string srcvarname, std::string trgvarname,
                            Portage::vector<std::vector<Weights_t>> const& sources_and_weights,
                            T lower_bound, T upper_bound,
                            Limiter_type limiter,
                            Partial_fixup_type partial_fixup_type,
                            Empty_fixup_type empty_fixup_type,
                            double conservation_tol,
                            int max_fixup_iter,
                            PartPair<D, ONWHAT,
                              SourceMesh, SourceState,
                              TargetMesh, TargetState>* parts_pair = nullptr) {
    assert(ONWHAT == onwhat());
    auto derived_class_ptr = static_cast<CoreDriverType<ONWHAT> *>(this);
    derived_class_ptr->
        template interpolate_mesh_var<T, Interpolate>(srcvarname, trgvarname,
                                                      sources_and_weights,
                                                      lower_bound, upper_bound,
                                                      limiter,
                                                      partial_fixup_type,
                                                      empty_fixup_type,
                                                      conservation_tol,
                                                      max_fixup_iter, parts_pair);
  }


  /*!
    Interpolate a (multi-)material variable of type T residing on CELLs
    
    @param[in] srcvarname   Variable name on source mesh

    @param[in] trgvarname   Variable name on target mesh

    @param[in] lower_bound  Lower bound for variable

    @param[in] upper_bound  Upper bound for variable

    @param[in] limiter      Limiter to use for second order reconstruction

    @param[in] partial_fixup_type Method to populate fields on
    partially filled target entities (cells or dual cells)

    @param[in] empty_fixup_type Method to populate fields on empty
    target entities (cells or dual cells)

    @param[in] conservation_tol Tolerance to which source and target
    integral quantities are to be matched

    @param[in] max_fixup_iter     Max number of iterations for global repair

    See support/portage.h for options on limiter, partial_fixup_type and
    empty_fixup_type
      
  */
  
  template <typename T = double,
            template<int, Entity_kind, class, class, class,
                     template <class, int, class, class> class,
                     class, class, class> class Interpolate
            >
  void interpolate_mat_var(std::string srcvarname, std::string trgvarname,
                           std::vector<Portage::vector<std::vector<Weights_t>>> const& sources_and_weights_by_mat,
                           T lower_bound, T upper_bound,
                           Limiter_type limiter,
                           Partial_fixup_type partial_fixup_type,
                           Empty_fixup_type empty_fixup_type,
                           double conservation_tol,
                           int max_fixup_iter) {

    assert(onwhat() == CELL);
    auto derived_class_ptr = static_cast<CoreDriverType<CELL> *>(this);
     derived_class_ptr->
         template interpolate_mat_var<T, Interpolate>(srcvarname,
                                                      trgvarname,
                                                      sources_and_weights_by_mat,
                                                      lower_bound, upper_bound,
                                                      limiter,
                                                      partial_fixup_type,
                                                      empty_fixup_type,
                                                      conservation_tol,
                                                      max_fixup_iter);
  }

  /*!
    @brief Check if meshes are mismatched (don't cover identical
    portions of space)

    @tparam Entity_kind  What kind of entity are we performing intersection of

    @param[in] sources_weights  Intersection sources and moments (vols, centroids) 
    @returns   Whether the meshes are mismatched
  */

  template<Entity_kind ONWHAT>
  bool 
  check_mesh_mismatch(Portage::vector<std::vector<Weights_t>> const& source_weights) {
    assert(ONWHAT == onwhat());
    auto derived_class_ptr = static_cast<CoreDriverType<ONWHAT> *>(this);
    return derived_class_ptr->check_mesh_mismatch(source_weights);
  }


  /*!
    @brief Set numerical tolerances for small volumes, distances, etc.

    @tparam Entity_kind  what kind of entity are we setting for

    @tparam num_tols     struct of selected numerical tolerances
  */

  template<Entity_kind ONWHAT>
  void
  set_num_tols(NumericTolerances_t num_tols) {
    assert(ONWHAT == onwhat());
    auto derived_class_ptr = static_cast<CoreDriverType<ONWHAT> *>(this);
    derived_class_ptr->set_num_tols(num_tols);
  }

};




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
class CoreDriver : public CoreDriverBase<D,
                                         SourceMesh, SourceState,
                                         TargetMesh, TargetState,
                                         InterfaceReconstructorType,
                                         Matpoly_Splitter,
                                         Matpoly_Clipper,
                                         CoordSys> {
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
      : source_mesh_(source_mesh), source_state_(source_state),
        target_mesh_(target_mesh), target_state_(target_state),
        executor_(executor)
  {
#ifdef PORTAGE_ENABLE_MPI
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
  ~CoreDriver() {}

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
  Portage::vector<std::vector<int>>
  search() {
    // Get an instance of the desired search algorithm type
    const Search<D, ONWHAT, SourceMesh, TargetMesh>
        search_functor(source_mesh_, target_mesh_);

    int ntarget_ents = target_mesh_.num_entities(ONWHAT, PARALLEL_OWNED);

    // initialize search candidate vector
    Portage::vector<std::vector<int>> candidates(ntarget_ents);

    Portage::transform(target_mesh_.begin(ONWHAT, PARALLEL_OWNED),
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
  Portage::vector<std::vector<Portage::Weights_t>>
  intersect_meshes(Portage::vector<std::vector<int>> const& candidates) {

    // Use default numerical tolerances in case they were not set earlier
    if (num_tols_.tolerances_set == false) {
        NumericTolerances_t default_num_tols;
        default_num_tols.use_default();
      set_num_tols(default_num_tols);
    }

    int nents = target_mesh_.num_entities(ONWHAT, PARALLEL_OWNED);
    Portage::vector<std::vector<Portage::Weights_t>> sources_and_weights(nents);
      
    Intersect<ONWHAT, SourceMesh, SourceState, TargetMesh,
              InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper>
        intersector(source_mesh_, source_state_, target_mesh_, num_tols_);

    Portage::transform(target_mesh_.begin(ONWHAT, PARALLEL_OWNED),
                       target_mesh_.end(ONWHAT, PARALLEL_OWNED),
                       candidates.begin(),
                       sources_and_weights.begin(),
                       intersector);

    return sources_and_weights;
  }


  /// Set numerical tolerances
  void set_num_tols(NumericTolerances_t num_tols) {
    num_tols_ = num_tols;
  }


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
  std::vector<Portage::vector<std::vector<Weights_t>>>
  intersect_materials(Portage::vector<std::vector<int>> const& candidates) {
      
    int nmats = source_state_.num_materials();
    assert(nmats > 1);

#ifdef HAVE_TANGRAM

    // Make sure we have a valid interface reconstruction method instantiated

    assert(typeid(InterfaceReconstructorType<SourceMesh, D,
                  Matpoly_Splitter, Matpoly_Clipper >) !=
           typeid(DummyInterfaceReconstructor<SourceMesh, D,
                  Matpoly_Splitter, Matpoly_Clipper>));

    std::vector<Tangram::IterativeMethodTolerances_t>
        tols(2, {1000, 1e-12, 1e-12});

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
                          Matpoly_Clipper>(source_mesh_, tols, true));
    
    int nsourcecells = source_mesh_.num_entities(CELL, ALL);
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

    std::vector<Portage::vector<std::vector<Weights_t>>>
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


      std::vector<std::vector<Weights_t>> this_mat_sources_and_wts(ntargetcells);
      Portage::transform(target_mesh_.begin(CELL, PARALLEL_OWNED),
                         target_mesh_.end(CELL, PARALLEL_OWNED),
                         candidates.begin(),
                         this_mat_sources_and_wts.begin(),
                         intersector);

      // LOOK AT INTERSECTION WEIGHTS TO DETERMINE WHICH TARGET CELLS
      // WILL GET NEW MATERIALS

      int ntargetcells = target_mesh_.num_entities(CELL, ALL);

      for (int c = 0; c < ntargetcells; c++) {
        std::vector<Weights_t> const& cell_mat_sources_and_weights =
            this_mat_sources_and_wts[c];
        int nwts = cell_mat_sources_and_weights.size();
        for (int s = 0; s < nwts; s++) {
          std::vector<double> const& wts = cell_mat_sources_and_weights[s].weights;
          if (wts[0] > 0.0) {
            double vol = target_mesh_.cell_volume(c);
            // Check that the volume of material we are adding to c is not miniscule
            if (wts[0]/vol > num_tols_.driver_relative_min_mat_vol) {
              matcellstgt.push_back(c);
              break;
            }
          }
        }
      }

      // If any processor is adding this material to the target state,
      // add it on all the processors

      int nmatcells = matcellstgt.size();
      int nmatcells_global = nmatcells;
#ifdef PORTAGE_ENABLE_MPI
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
    return std::vector<Portage::vector<std::vector<Weights_t>>>();
#endif

  }



  /*! 
    Check mismatch between meshes

    @param[in] sources_and_weights Intersection sources and moments
    (vols, centroids)

    @returns   Whether the meshes are mismatched
  */

  bool
  check_mesh_mismatch(Portage::vector<std::vector<Weights_t>> const& source_weights) {

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
                                           source_weights, executor_));
    }

    return mismatch_fixer_->has_mismatch();
  }

  /**
   * @brief Interpolate mesh variable.
   *
   * @param[in] srcvarname          source mesh variable to remap
   * @param[in] trgvarname          target mesh variable to remap
   * @param[in] sources_and_weights weights for mesh-mesh interpolation
   * @param[in] lower_bound         lower bound of variable value when doing fixup
   * @param[in] upper_bound         upper bound of variable value when doing fixup
   * @param[in] limiter             limiter to use
   * @param[in] partial...          how to fixup partly filled target cells
   * @param[in] emtpy...            how to fixup empty target cells with this var
   * @param[in] cons..tol           tolerance for conservation when doing fixup
   * @param[in] max_fixup_iter      maximum number of iterations for mismatch fixup
   * @param[in] partition           source and target entities list for part-by-part
   */
  template<typename T = double,
    template<int, Entity_kind, class, class, class,
    template<class, int, class, class> class,
    class, class, class> class Interpolate
  >
  void interpolate_mesh_var(std::string srcvarname, std::string trgvarname,
                            Portage::vector<std::vector<Weights_t>> const& sources_and_weights,
                            T lower_bound, T upper_bound,
                            Limiter_type limiter = DEFAULT_LIMITER,
                            Partial_fixup_type partial_fixup_type = DEFAULT_PARTIAL_FIXUP_TYPE,
                            Empty_fixup_type empty_fixup_type = DEFAULT_EMPTY_FIXUP_TYPE,
                            double conservation_tol = DEFAULT_CONSERVATION_TOL,
                            int max_fixup_iter = DEFAULT_MAX_FIXUP_ITER,
                            const PartPair<D, ONWHAT,
                                           SourceMesh, SourceState,
                                           TargetMesh, TargetState> *const partition = nullptr) {

    if (source_state_.get_entity(srcvarname) != ONWHAT) {
      std::cerr << "Variable " << srcvarname << " not defined on Entity_kind "
                << ONWHAT << ". Skipping!" << std::endl;
      return;
    }

    // useful shortcuts
    using entity_weights_t = std::vector<Weights_t>;
    using interpolator_t =
      Interpolate<D, ONWHAT, SourceMesh, TargetMesh, SourceState,
        InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper, CoordSys>;

    interpolator_t interpolator(source_mesh_, target_mesh_, source_state_, num_tols_);
    interpolator.set_interpolation_variable(srcvarname, limiter);

    // get a handle to a memory location where the target state
    // would like us to write this material variable into.
    T* target_mesh_field = nullptr;
    target_state_.mesh_get_data(ONWHAT, trgvarname, &target_mesh_field);

    // perform part-by-part interpolation
    if (partition != nullptr) {

      // 1. Do some basic checks on supplied source and target parts
      // to prevent bugs when interpolating values:
      // - first check that source and target parts are not empty,
      // - then check that parts mismatch has already been tested.
      // - afterwards, check that each entity id is within the
      //   mesh entity index space.
      assert(partition->source_part_size() > 0);
      assert(partition->target_part_size() > 0);
      assert(partition->is_mismatch_tested());

      int const& max_source_id = source_mesh_.num_entities(ONWHAT, ALL);
      int const& max_target_id = target_mesh_.num_entities(ONWHAT, ALL);

      auto const& part_source_entities = partition->get_source_entities();
      auto const& part_target_entities = partition->get_target_entities();

      Portage::for_each(part_source_entities.begin(),
                        part_source_entities.end(),
                        [&](int current){ assert(current <= max_source_id); });

      Portage::for_each(part_target_entities.begin(),
                        part_target_entities.end(),
                        [&](int current){ assert(current <= max_target_id); });

      int const target_mesh_size = sources_and_weights.size();
      int const target_part_size = partition->target_part_size();

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
          if(partition->is_source_entity(weight.entityID)) {
            heap.emplace_back(weight);
          }
        }
        heap.shrink_to_fit();
        return heap;
      };

      Portage::vector<entity_weights_t> parts_weights(target_part_size);
      Portage::transform(part_target_entities.begin(), part_target_entities.end(),
                         parts_weights.begin(), filter_weights);

      // 3. Process interpolation.
      // Now that intersection weights is filtered, perform the interpolation.
      // Notice that we need to store the interpolated values in a temporary
      // array since they are indexed with respect to the target part
      // and not to the target mesh. Hence we need to copy them back in the
      // target state at their correct (absolute index) memory locations.
      T temporary_storage[target_part_size];
      Portage::pointer<T> target_part_field(temporary_storage);

      Portage::transform(part_target_entities.begin(), part_target_entities.end(),
                         parts_weights.begin(), target_part_field, interpolator);

      for (int i=0; i < target_part_size; ++i) {
        auto const& j = part_target_entities[i];
        target_mesh_field[j] = target_part_field[i];
      }

      // 4. Fix partially filled and empty cells values if necessary.
      // Notice that mismatch detection should have been already performed.
      if (partition->has_mismatch()) {
        #ifdef DEBUG
          std::fprintf(stderr,
            "There is a mismatch between source and target sub-meshes\n"
            "Will start fixing interpolated values\n"
          );
        #endif
        partition->fix_mismatch(srcvarname, trgvarname,
                                lower_bound, upper_bound,
                                conservation_tol, max_fixup_iter,
                                partial_fixup_type, empty_fixup_type);
      }
    } else /* mesh-mesh interpolation */ {
      Portage::pointer<T> target_field(target_mesh_field);
      Portage::transform(target_mesh_.begin(ONWHAT, PARALLEL_OWNED),
                         target_mesh_.end(ONWHAT, PARALLEL_OWNED),
                         sources_and_weights.begin(),
                         target_field, interpolator);

      assert(mismatch_fixer_ && "check_mesh_mismatch must be called first");
      if (mismatch_fixer_->has_mismatch()) {
        mismatch_fixer_->fix_mismatch(srcvarname, trgvarname,
                                      lower_bound, upper_bound,
                                      conservation_tol, max_fixup_iter,
                                      partial_fixup_type, empty_fixup_type);
      }
    }
  }


#ifdef HAVE_TANGRAM
  
  /*! CoreDriver::interpolate_mat_var

    @brief interpolate a material variable

    @param[in] srcvarname  Material variable name on the source mesh

    @param[in] trgvarname  Material variable name on the target mesh

    @param[in] limiter     Limiter to use for variable

    @param[in] partial...  How to fixup partly filled target cells (for this var)

    @param[in] emtpy...    How to fixup empty target cells with this var

    @param[in] lower_bound Lower bound of variable value when doing fixup

    @param[in] upper_bound Upper bound of variable value when doing fixup

    @param[in] cons..tol   Tolerance for conservation when doing fixup
  */

  template<typename T = double,
           template<int, Entity_kind, class, class, class,
                    template<class, int, class, class> class,
                    class, class, class> class Interpolate
           >
  void interpolate_mat_var(std::string srcvarname, std::string trgvarname,
                           std::vector<Portage::vector<std::vector<Weights_t>>> const& sources_and_weights_by_mat,
                           T lower_bound, T upper_bound,
                           Limiter_type limiter = DEFAULT_LIMITER,
                           Partial_fixup_type partial_fixup_type = DEFAULT_PARTIAL_FIXUP_TYPE,
                           Empty_fixup_type empty_fixup_type = DEFAULT_EMPTY_FIXUP_TYPE,
                           double conservation_tol = DEFAULT_CONSERVATION_TOL,
                           int max_fixup_iter = DEFAULT_MAX_FIXUP_ITER) {
    
    Interpolate<D, ONWHAT, SourceMesh, TargetMesh, SourceState,
                InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper, CoordSys>
        interpolator(source_mesh_, target_mesh_, source_state_, num_tols_, interface_reconstructor_);
      
    int nmats = source_state_.num_materials();


    for (int m = 0; m < nmats; m++) {

      interpolator.set_material(m);    // We have to do this so we know
      //                                 // which material values we have
      //                                 // to grab from the source state

      // FEATURE ;-)  Have to set interpolation variable AFTER setting 
      // the material for multimaterial variables

      interpolator.set_interpolation_variable(srcvarname, limiter);

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

      Portage::pointer<T> target_field(target_field_raw);

      Portage::transform(matcellstgt.begin(), matcellstgt.end(),
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

#endif  // HAVE_TANGRAM

  
 private:
  SourceMesh const & source_mesh_;
  TargetMesh const & target_mesh_;
  SourceState const & source_state_;
  TargetState & target_state_;

  NumericTolerances_t num_tols_;

  int comm_rank_ = 0;
  int nprocs_ = 1;

  Wonton::Executor_type const *executor_;

#ifdef PORTAGE_ENABLE_MPI
  MPI_Comm mycomm_ = MPI_COMM_NULL;
#endif
  std::unique_ptr<MismatchFixer<D, ONWHAT,
                                SourceMesh, SourceState,
                                TargetMesh, TargetState>> mismatch_fixer_;

#ifdef HAVE_TANGRAM

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

    int nvals = 0;
    for (int m = 0; m < nmats; m++) {
      std::vector<int> cellids;
      source_state_.mat_get_cells(m, &cellids);
      for (int ic = 0; ic < cellids.size(); ic++) {
        int c = cellids[ic];
        int nmatc = cell_num_mats[c];
        cell_mat_ids_full[c*nmats+nmatc] = m;
        cell_num_mats[c]++;
      }
      nvals += cellids.size();

      double const * matfracptr;
      source_state_.mat_get_celldata("mat_volfracs", m, &matfracptr);
      for (int ic = 0; ic < cellids.size(); ic++)
        cell_mat_volfracs_full[cellids[ic]*nmats+m] = matfracptr[ic];

      Wonton::Point<D> const *matcenvec;
      source_state_.mat_get_celldata("mat_centroids", m, &matcenvec);
      for (int ic = 0; ic < cellids.size(); ic++)
        cell_mat_centroids_full[cellids[ic]*nmats+m] = matcenvec[ic];
    }

    // At this point nvals contains the number of non-zero volume
    // fraction entries in the full array. Use this and knowledge of
    // number of materials in each cell to compress the data into
    // linear arrays

    cell_mat_ids.resize(nvals);
    cell_mat_volfracs.resize(nvals);
    cell_mat_centroids.resize(nvals);

    int idx = 0;
    for (int c = 0; c < nsourcecells; c++) {
      for (int m = 0; m < cell_num_mats[c]; m++) {
        int matid = cell_mat_ids_full[c*nmats+m];
        cell_mat_ids[idx] = matid;
        cell_mat_volfracs[idx] = cell_mat_volfracs_full[c*nmats+matid];
        cell_mat_centroids[idx] = cell_mat_centroids_full[c*nmats+matid];
        idx++;
      }
    }
  }

#endif

};  // CoreDriverBase


// Core Driver Factory

template <int D,
          class SourceMesh, class SourceState,
          class TargetMesh, class TargetState,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper,
          class CoordSys>
std::unique_ptr<CoreDriverBase<D, SourceMesh, SourceState,
                               TargetMesh, TargetState,
                               InterfaceReconstructorType,
                               Matpoly_Splitter, Matpoly_Clipper,
                               CoordSys>
                >
make_core_driver(Entity_kind onwhat,
                 SourceMesh const & source_mesh,
                 SourceState const & source_state,
                 TargetMesh const & target_mesh,
                 TargetState & target_state,
                 Wonton::Executor_type const *executor) {
  switch (onwhat) {
    case CELL:
      // Intel 18.0.1 does not recognize std::make_unique even with -std=c++14 flag *ugh*
      // return std::make_unique<CoreDriver<D, CELL,
      //                                        SourceMesh, SourceState,
      //                                        TargetMesh, TargetState,
      //                                        InterfaceReconstructorType,
      //                                        Matpoly_Splitter, Matpoly_Clipper,
      //                                        CoordSys>>
      //     (source_mesh, source_state, target_mesh, target_state, executor);

      return std::unique_ptr<CoreDriver<D, CELL,
                                        SourceMesh, SourceState,
                                        TargetMesh, TargetState,
                                        InterfaceReconstructorType,
                                        Matpoly_Splitter, Matpoly_Clipper,
                                        CoordSys>>(
                                            new CoreDriver<D, CELL,
                                            SourceMesh, SourceState,
                                            TargetMesh, TargetState,
                                            InterfaceReconstructorType,
                                            Matpoly_Splitter, Matpoly_Clipper,
                                            CoordSys>
                                            (source_mesh, source_state, 
                                             target_mesh, target_state,
                                             executor)
                                                        );
    case NODE:
      // Intel 18.0.1 does not recognize std::make_unique even with -std=c++14 flag *ugh*
      // return std::make_unique<CoreDriver<D, NODE,
      //                                        SourceMesh, SourceState,
      //                                        TargetMesh, TargetState,
      //                                        InterfaceReconstructorType,
      //                                        Matpoly_Splitter, Matpoly_Clipper,
      //                                        CoordSys>>
      //     (source_mesh, source_state, target_mesh, target_state, executor);
      return std::unique_ptr<CoreDriver<D, CELL,
                                        SourceMesh, SourceState,
                                        TargetMesh, TargetState,
                                        InterfaceReconstructorType,
                                        Matpoly_Splitter, Matpoly_Clipper,
                                        CoordSys>>(
                                            new CoreDriver<D, CELL,
                                            SourceMesh, SourceState,
                                            TargetMesh, TargetState,
                                            InterfaceReconstructorType,
                                            Matpoly_Splitter, Matpoly_Clipper,
                                            CoordSys>
                                            (source_mesh, source_state, 
                                             target_mesh, target_state,
                                             executor)
                                                   );
    default:
      std::cerr << "Remapping on " << Wonton::to_string(onwhat) <<
          " not implemented\n";
  }
}  // make_core_driver


}  // namespace Portage

#endif  // PORTAGE_CORE_DRIVER_H_
