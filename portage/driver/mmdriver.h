/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_MMDRIVER_H_
#define PORTAGE_MMDRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"

#include "portage/intersect/dummy_interface_reconstructor.h"
#endif

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wonton/state/flat/flat_state_wrapper.h"

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

using namespace Wonton;

/*!
  @class MMDriver "mmdriver.h"
  @brief MMDriver provides the API to mapping multi-material data from one mesh to another.

  @tparam Search  A search method that takes the dimension, source mesh class
  and target mesh class as template parameters

  @tparam Intersect  A polyhedron-polyhedron intersection class that takes
  the source and taget mesh classes as template parameters

  @tparam Interpolate An interpolation class that takes the source and
  target mesh classes, the source state class (that stores source
  field values), the kind of entity the interpolation is on and the
  dimension of the problem as template parameters

  @tparam SourceMesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.

  @tparam SourceState_Wrapper A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.

  @tparam TargetMesh_Wrapper A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.

  @tparam TargetState_Wrapper A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.

  @tparam InterfaceReconstructorType An interface reconstruction class
  that takes the raw interface reconstruction method, the dimension of
  the problem and the source mesh class as template parameters
*/
template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
          template <class, int, class, class> class, 
          class, class> class Intersect,
          template<int, Entity_kind, class, class, class,
          template<class, int, class, class> class,
          class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper = SourceMesh_Wrapper,
          class TargetState_Wrapper = SourceState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType = DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void>
class MMDriver {

  // Something like this would be very helpful to users
  // static_assert(
  //   D == Interpolate::D,
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
  MMDriver(SourceMesh_Wrapper const& sourceMesh,
           SourceState_Wrapper const& sourceState,
           TargetMesh_Wrapper const& targetMesh,
           TargetState_Wrapper& targetState)
      : source_mesh_(sourceMesh), source_state_(sourceState),
        target_mesh_(targetMesh), target_state_(targetState),
        dim_(sourceMesh.space_dimension()) {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
  }

  /// Copy constructor (disabled)
  MMDriver(const MMDriver &) = delete;

  /// Assignment operator (disabled)
  MMDriver & operator = (const MMDriver &) = delete;

  /// Destructor
  ~MMDriver() {}

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
    for (int i = 0; i < nvars; ++i) {
      Entity_kind srckind = source_state_.get_entity(source_remap_var_names[i]);
      Entity_kind trgkind = target_state_.get_entity(target_remap_var_names[i]);
      if (trgkind == Entity_kind::UNKNOWN_KIND)
        continue;  // Presumably field does not exist on target - will get added
      assert(srckind == trgkind);  // if target field exists, entity kinds
                                   // must match
    }

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
    @brief remap for a given set of MESH and MATERIAL variables on a given entity kind
    @tparam entity_kind  Kind of entity that variables live on
    @param source_meshvar_names  names of remap variables on source mesh
    @param target_meshvar_names  names of remap variables on target mesh
    @param source_matvar_names  names of remap variables on materials of source mesh
    @param target_matvar_names  names of remap variables on materials of target mesh
    @return status of remap (1 if successful, 0 if not)
  */

  template<Entity_kind onwhat>
  int remap(std::vector<std::string> const &source_meshvar_names,
            std::vector<std::string> const &target_meshvar_names,
            std::vector<std::string> const &source_matvar_names,
            std::vector<std::string> const &target_matvar_names);

#ifdef ENABLE_MPI
  /*!
    @brief remap for a given set of variables on a given entity kind in a distributed setting with redistribution of data if needed
    @tparam entity_kind  Kind of entity that variables live on
    @param source_meshvar_names  names of remap variables on source mesh
    @param target_meshvar_names  names of remap variables on target mesh
    @param source_matvar_names  names of remap variables on materials of source mesh
    @param target_matvar_names  names of remap variables on materials of target mesh
    @return status of remap (1 if successful, 0 if not)
  */

  template<Entity_kind onwhat>
  int remap_distributed(std::vector<std::string> const &source_meshvar_names,
                        std::vector<std::string> const &target_meshvar_names,
                        std::vector<std::string> const &source_matvar_names,
                        std::vector<std::string> const &target_matvar_names);
#endif

  

  /*!
    @brief Execute the remapping process
    @return status of remap (1 if successful, 0 if not)
  */
  int run(bool distributed, std::string *errmsg = nullptr) {
    std::string mesg;
#ifndef ENABLE_MPI
    if (distributed) {
      mesg = "Request is for a parallel run but Portage is compiled for serial runs only";
      if (errmsg)
        *errmsg = mesg;
      else
        std::cerr << errmsg << "\n";
      return 0;
    }
#endif

    int comm_rank = 0;
#ifdef ENABLE_MPI
    if (distributed)
      MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

    if (comm_rank == 0)
      std::cout << "in MMDriver::run()...\n";

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;


    int nvars = source_remap_var_names_.size();


    std::vector<std::string> src_meshvar_names, src_matvar_names;
    std::vector<std::string> trg_meshvar_names, trg_matvar_names;
    

    // -------- CELL VARIABLE REMAP ---------
    // Collect all cell based variables and remap them

    for (int i = 0; i < nvars; ++i) {
      std::string& srcvarname = source_remap_var_names_[i];
      Entity_kind onwhat = source_state_.get_entity(srcvarname);
      if (onwhat == CELL) {
        // Separate out mesh fields and multi-material fields - they will be
        // processed differently

        std::string& trgvarname = target_remap_var_names_[i];

        Field_type ftype = source_state_.field_type(onwhat, srcvarname);

        if (ftype == Field_type::MESH_FIELD) {
          src_meshvar_names.push_back(srcvarname);
          trg_meshvar_names.push_back(trgvarname);
        } else if (ftype == Field_type::MULTIMATERIAL_FIELD) {
          src_matvar_names.push_back(srcvarname);
          trg_matvar_names.push_back(trgvarname);
        }
      }
    }

    // ALWAYS call because we may have to remap material volume
    // fractions and centroids which are cell-based fields
    
#ifdef ENABLE_MPI
    if (distributed)
      remap_distributed<CELL>(src_meshvar_names, trg_meshvar_names,
                              src_matvar_names, trg_matvar_names);
    else
#endif
      remap<CELL>(src_meshvar_names, trg_meshvar_names,
                  src_matvar_names, trg_matvar_names);



    // -------- NODE VARIABLE REMAP ---------
    // Collect all node based variables and remap them
    // (ignore any multi-material variables on NODES - not well defined)

    src_meshvar_names.clear(); src_matvar_names.clear();
    trg_meshvar_names.clear(); trg_matvar_names.clear();
    
    for (int i = 0; i < nvars; ++i) {
      std::string& srcvarname = source_remap_var_names_[i];
      Entity_kind onwhat = source_state_.get_entity(srcvarname);
      if (onwhat == NODE) {
        std::string& trgvarname = target_remap_var_names_[i];
        
        Field_type ftype = source_state_.field_type(onwhat, srcvarname);
        
        if (ftype == Field_type::MESH_FIELD) {
          src_meshvar_names.push_back(srcvarname);
          trg_meshvar_names.push_back(trgvarname);
        } else if (ftype == Field_type::MULTIMATERIAL_FIELD)
          std::cerr << "Cannot handle multi-material fields on nodes\n";
      }
    }

    if (src_meshvar_names.size()) {
#ifdef ENABLE_MPI
      if (distributed)
        remap_distributed<NODE>(src_meshvar_names, trg_meshvar_names,
                                src_matvar_names, trg_matvar_names);
      else
#endif
        remap<NODE>(src_meshvar_names, trg_meshvar_names,
                    src_meshvar_names, trg_meshvar_names);
    }

    return 1;
  }  // run


  
 private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  std::vector<LimiterType> limiters_;
  unsigned int dim_;

  // Convert volume fraction and centroid data from compact
  // material-centric to compact cell-centric (ccc) form as needed by
  // Tangram

  void ccc_vfcen_data(std::vector<int>& cell_num_mats,
                      std::vector<int>& cell_mat_ids,
                      std::vector<double>& cell_matvolfracs,
                      std::vector<Tangram::Point<D>>& cell_mat_centroids);
  
};  // class MMDriver



// Serial remap or Distributed remap with no redistributon of data

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
          template <class, int, class, class> class,
          class, class> class Intersect,
          template<int, Entity_kind, class, class, class,
          template<class, int, class, class> class,
          class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper>
template<Entity_kind onwhat>
int MMDriver<Search, Intersect, Interpolate, D,
             SourceMesh_Wrapper, SourceState_Wrapper,
             TargetMesh_Wrapper, TargetState_Wrapper,
             InterfaceReconstructorType, Matpoly_Splitter,
             Matpoly_Clipper
             >::remap(std::vector<std::string> const &src_meshvar_names,
                      std::vector<std::string> const &trg_meshvar_names,
                      std::vector<std::string> const &src_matvar_names,
                      std::vector<std::string> const &trg_matvar_names) {

  static_assert(onwhat == NODE || onwhat == CELL,
                "Remap implemented only for CELL and NODE variables");

  int comm_rank = 0;

  int ntarget_ents_owned = target_mesh_.num_entities(onwhat, PARALLEL_OWNED);
  std::cout << "Number of target entities of kind " << onwhat <<
      " in target mesh on rank " << comm_rank << ": " <<
      ntarget_ents_owned << std::endl;

  int ntarget_ents = target_mesh_.num_entities(onwhat, ALL);

  float tot_seconds = 0.0, tot_seconds_srch = 0.0,
      tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;


  
  // SEARCH

  Portage::vector<std::vector<int>> candidates(ntarget_ents);
  Portage::vector<std::vector<Weights_t>> source_ents_and_weights(ntarget_ents);

  // Get an instance of the desired search algorithm type
  gettimeofday(&begin_timeval, 0);
  const Search<D, onwhat, SourceMesh_Wrapper, TargetMesh_Wrapper>
      search(source_mesh_, target_mesh_);

  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(), search);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  int nmats = source_state_.num_materials();

#ifdef HAVE_TANGRAM
  // Call interface reconstruction only if we got a method from the
  // calling app
  Tangram::IterativeMethodTolerances_t tol{200, 1e-12, 1e-12}; 
 
  auto interface_reconstructor =
      std::make_shared<Tangram::Driver<InterfaceReconstructorType, D,
                                       SourceMesh_Wrapper,
                                       Matpoly_Splitter,
                                       Matpoly_Clipper>
                       >(source_mesh_, tol, true);
    
  if (typeid(InterfaceReconstructorType<SourceMesh_Wrapper, D,
             Matpoly_Splitter, Matpoly_Clipper >) !=
      typeid(DummyInterfaceReconstructor<SourceMesh_Wrapper, D,
             Matpoly_Splitter, Matpoly_Clipper>)) {
    
    int nsourcecells = source_mesh_.num_entities(CELL, ALL);

    std::vector<int> cell_num_mats;
    std::vector<int> cell_mat_ids;
    std::vector<double> cell_mat_volfracs;
    std::vector<Tangram::Point<D>> cell_mat_centroids;

    // Extract volume fraction and centroid data for cells in compact
    // cell-centric form (ccc)

    ccc_vfcen_data(cell_num_mats, cell_mat_ids, cell_mat_volfracs,
                   cell_mat_centroids);
      
    interface_reconstructor->set_volume_fractions(cell_num_mats,
                                                  cell_mat_ids,
                                                  cell_mat_volfracs,
                                                  cell_mat_centroids);
    interface_reconstructor->reconstruct();
  }

  
  // Make an intersector which knows about the source state (to be able
  // to query the number of materials, etc) and also knows about the
  // interface reconstructor so that it can retrieve pure material polygons

  Intersect<onwhat, SourceMesh_Wrapper, SourceState_Wrapper,
            TargetMesh_Wrapper, InterfaceReconstructorType, 
            Matpoly_Splitter, Matpoly_Clipper>
      intersect(source_mesh_, source_state_, target_mesh_,
                interface_reconstructor);

  // Get an instance of the desired interpolate algorithm type
  Interpolate<D, onwhat, SourceMesh_Wrapper, TargetMesh_Wrapper,
              SourceState_Wrapper, InterfaceReconstructorType,
              Matpoly_Splitter, Matpoly_Clipper>
      interpolate(source_mesh_, target_mesh_, source_state_, 
                  interface_reconstructor);
#else

  Intersect<onwhat, SourceMesh_Wrapper, SourceState_Wrapper,
            TargetMesh_Wrapper, DummyInterfaceReconstructor,
            void, void>
      intersect(source_mesh_, source_state_, target_mesh_);

  // Get an instance of the desired interpolate algorithm type
  Interpolate<D, onwhat, SourceMesh_Wrapper, TargetMesh_Wrapper,
              SourceState_Wrapper, DummyInterfaceReconstructor,
              void, void>
      interpolate(source_mesh_, target_mesh_, source_state_);
#endif  // HAVE_TANGRAM

  
  //--------------------------------------------------------------------
  // REMAP MESH FIELDS FIRST (this requires just mesh-mesh intersection)
  //--------------------------------------------------------------------
  
  // INTERSECT
  
  gettimeofday(&begin_timeval, 0);
  
  // For each cell in the target mesh get a list of candidate-weight
  // pairings (in a traditional mesh, not particle mesh, the weights
  // are moments). Note that this candidate list is different from the
  // search candidate list in that it may not include some of the
  // search candidates. Also, note that for 2nd order and higher
  // remaps, we get multiple moments (0th, 1st, etc) for each
  // target-source cell intersection
  
  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(),
                     source_ents_and_weights.begin(),
                     intersect);
  
  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_xsect += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
  

  // INTERPOLATE (one variable at a time)

  gettimeofday(&begin_timeval, 0);
  
  int nvars = src_meshvar_names.size();
  if (comm_rank == 0)
    std::cout << "Number of mesh variables on entity kind " << onwhat <<
        " to remap is " << nvars << std::endl;
  
  for (int i = 0; i < nvars; ++i) {
    interpolate.set_interpolation_variable(src_meshvar_names[i], limiters_[i]);
    
    // Get a handle to a memory location where the target state
    // would like us to write this material variable into.
    
    double *target_field_raw;
    target_state_.mesh_get_data(onwhat, trg_meshvar_names[i], &target_field_raw);
    assert (target_field_raw != nullptr);
    
    
    Portage::pointer<double> target_field(target_field_raw);
    
    Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                       target_mesh_.end(onwhat, PARALLEL_OWNED),
                       source_ents_and_weights.begin(),
                       target_field, interpolate);
  }
  
  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;


  
  //--------------------------------------------------------------------
  // REMAP MULTIMATERIAL FIELDS NEXT, ONE MATERIAL AT A TIME
  //--------------------------------------------------------------------

  if (onwhat != CELL) return 1;
  
  // Material centric loop

  for (int m = 0; m < nmats; m++) {

    // INTERSECT
    
    gettimeofday(&begin_timeval, 0);

    intersect.set_material(m);

    // For each cell in the target mesh get a list of candidate-weight
    // pairings (in a traditional mesh, not particle mesh, the weights
    // are moments). Note that this candidate list is different from the
    // search candidate list in that it may not include some of the
    // search candidates. Also, note that for 2nd order and higher
    // remaps, we get multiple moments (0th, 1st, etc) for each
    // target-source cell intersection

    // NOTE: IDEALLY WE WOULD REUSE THE MESH-MESH INTERSECTIONS FROM THE
    // PREVIOUS STEP WHEN THE SOURCE MATERIAL CONTAINS ONLY ONE MATERIAL
    //
    // UNFORTUNATELY, THE REQUIREMENT OF THE INTERSECT FUNCTOR IS THAT
    // IT CANNOT MODIFY STATE, THIS MEANS WE CANNOT STORE THE MESH-MESH
    // INTERSECTION VALUES AND REUSE THEM AS NECESSARY FOR MESH-MATERIAL
    // INTERSECTION COMPUTATIONS
    
    Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                       target_mesh_.end(onwhat, PARALLEL_OWNED),
                       candidates.begin(),
                       source_ents_and_weights.begin(),
                       intersect);
    
    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    tot_seconds_xsect += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

    // LOOK AT INTERSECTION WEIGHTS TO DETERMINE WHICH TARGET CELLS
    // WILL GET NEW MATERIALS
    
    int ntargetcells = target_mesh_.num_entities(CELL, ALL);
    std::vector<int> matcellstgt;
    
    for (int c = 0; c < ntargetcells; c++) {
      std::vector<Weights_t> const& cell_sources_and_weights =
          source_ents_and_weights[c];
      for (int s = 0; s < cell_sources_and_weights.size(); s++) {
        std::vector<double> const& wts = cell_sources_and_weights[s].weights;
        if (wts[0] > 0.0) {
          double vol = target_mesh_.cell_volume(c);
          if (wts[0]/vol > 1.0e-10) {  // Check that the volume of material
                                       // we are adding to c is not miniscule
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
        // CELLS TO THE STATEMANAGER, IT ALSO MAKES SPACE FOR FIELD
        // VALUES FOR THIS MATERIAL IN EVERY MULTI-MATERIAL VECTOR IN
        // THE STATE MANAGER. THIS ENSURES THAT WHEN WE CALL
        // mat_get_celldata FOR A MATERIAL IN MULTI-MATERIAL STATE
        // VECTOR IT WILL ALREADY HAVE SPACE ALLOCATED FOR FIELD
        // VALUES OF THAT MATERIAL. SOME STATE WRAPPERS COULD CHOOSE
        // TO MAKE THIS A SIMPLER ROUTINE THAT ONLY STORES THE NAME
        // AND THE CELLS IN THE MATERIAL AND ACTUALLY ALLOCATE SPACE
        // FOR FIELD VALUES OF A MATERIAL IN A MULTI-MATERIAL FIELD
        // WHEN mat_get_celldata IS INVOKED.
        
        target_state_.add_material(source_state_.material_name(m), matcellstgt);
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
    std::vector<std::vector<Weights_t>> mat_sources_and_weights(nmatcells);
        
    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcellstgt[ic];
      double matvol = 0.0;
      Point<D> matcen;
      std::vector<Weights_t> const& cell_sources_and_weights =
          source_ents_and_weights[c];
      for (int s = 0; s < cell_sources_and_weights.size(); s++) {
        std::vector<double> const& wts = cell_sources_and_weights[s].weights;
        matvol += wts[0];
        for (int d = 0; d < D; d++)
          matcen[d] += wts[d+1];
      }
      matcen /= matvol;
      mat_volfracs[ic] = matvol/target_mesh_.cell_volume(c);
      mat_centroids[ic] = matcen;

      mat_sources_and_weights[ic] = cell_sources_and_weights;
    }

    target_state_.mat_add_celldata("mat_volfracs", m, &(mat_volfracs[0]));
    target_state_.mat_add_celldata("mat_centroids", m, &(mat_centroids[0]));


    // INTERPOLATE (one variable at a time)

    // HERE WE COULD MAKE A NEW LIST BASED ON WHICH TARGET CELLS HAVE ANY
    // INTERSECTIONS WITH SOURCE CELLS FOR THIS MATERIAL TO AVOID A NULL-OP
    // AND A WARNING MESSAGE ABOUT NO SOURCE CELLS CONTRIBUTING TO A TARGET -
    // IS IT WORTH IT?

    gettimeofday(&begin_timeval, 0);

    int nmatvars = src_matvar_names.size();
    if (comm_rank == 0)
      std::cout << "Number of multi-material variables on entity kind " <<
          onwhat << " to remap is " << nmatvars << std::endl;

    interpolate.set_material(m);    // We have to do this so we know
                                    // which material values we have
                                    // to grab from the source state

    for (int i = 0; i < nmatvars; ++i) {
      interpolate.set_interpolation_variable(src_matvar_names[i], limiters_[i]);

      // Get a handle to a memory location where the target state
      // would like us to write this material variable into. If it is
      // NULL, we allocate it ourself

      double *target_field_raw;
      target_state_.mat_get_celldata(trg_matvar_names[i], m, &target_field_raw);
      assert (target_field_raw != nullptr);


      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform(matcellstgt.begin(), matcellstgt.end(),
                         mat_sources_and_weights.begin(),
                         target_field, interpolate);

      // If the state wrapper knows that the target data is already
      // laid out in this way and it gave us a pointer to the array
      // where the values reside, it has to do nothing in this
      // call. If the storage format is different, however, it may
      // have to copy the values into their proper locations

      target_state_.mat_add_celldata(trg_matvar_names[i], m, target_field_raw);
    }  // nmatvars

    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  }  // for nmats

  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

  std::cout << "Transform Time for Entity Kind " << onwhat << " on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;

  return 1;
}



#ifdef ENABLE_MPI

// Distributed Remap with redistribution of mesh and data

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
          template <class, int, class, class> class,
          class, class> class Intersect,
          template<int, Entity_kind, class, class, class,
          template<class, int, class, class> class,
          class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper>
template<Entity_kind onwhat>
int MMDriver<Search, Intersect, Interpolate, D,
             SourceMesh_Wrapper, SourceState_Wrapper,
             TargetMesh_Wrapper, TargetState_Wrapper,
             InterfaceReconstructorType, Matpoly_Splitter, 
             Matpoly_Clipper
             >::remap_distributed(std::vector<std::string> const &src_meshvar_names,
                                  std::vector<std::string> const &trg_meshvar_names,
                                  std::vector<std::string> const &src_matvar_names,
                                  std::vector<std::string> const &trg_matvar_names) {

  static_assert(onwhat == NODE || onwhat == CELL,
                "Remap implemented only for CELL and NODE variables");

  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  int ntarget_ents_owned = target_mesh_.num_entities(onwhat, PARALLEL_OWNED);
  std::cout << "Number of target entities of kind " << onwhat <<
      " in target mesh on rank " << comm_rank << ": " <<
      ntarget_ents_owned << std::endl;

  int ntarget_ents = target_mesh_.num_entities(onwhat, ALL);

  Flat_Mesh_Wrapper<> source_mesh_flat;
  Flat_State_Wrapper<> source_state_flat;

  float tot_seconds = 0.0, tot_seconds_srch = 0.0,
      tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;


  // SEARCH

  Portage::vector<std::vector<int>> candidates(ntarget_ents);
  Portage::vector<std::vector<Weights_t>> source_ents_and_weights(ntarget_ents);

  // Create flat wrappers to distribute source cells
  gettimeofday(&begin_timeval, 0);

  source_mesh_flat.initialize(source_mesh_);
  source_state_flat.initialize(source_state_, source_remap_var_names_);
  MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat,
                         target_mesh_, target_state_);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  float tot_seconds_flat = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
  std::cout << "Redistribution Time Rank " << comm_rank << " (s): " <<
      tot_seconds_flat << std::endl;

  // Get an instance of the desired search algorithm type
  gettimeofday(&begin_timeval, 0);
  const Search<D, onwhat, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>
      search(source_mesh_flat, target_mesh_);

  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(), search);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  int nmats = source_state_.num_materials();

#ifdef HAVE_TANGRAM
  // Call interface reconstruction only if we got a method from the
  // calling app
  Tangram::IterativeMethodTolerances_t tol{200, 1e-12, 1e-12}; 
 
  auto interface_reconstructor =
      std::make_shared<Tangram::Driver<InterfaceReconstructorType, D,
                                       Flat_Mesh_Wrapper<>,
                                       Matpoly_Splitter,
                                       Matpoly_Clipper>
                       >(source_mesh_flat, tol, true);
    
  if (typeid(InterfaceReconstructorType<SourceMesh_Wrapper, D,
             Matpoly_Splitter, Matpoly_Clipper >) !=
      typeid(DummyInterfaceReconstructor<SourceMesh_Wrapper, D,
             Matpoly_Splitter, Matpoly_Clipper>)) {
    
    int nsourcecells = source_mesh_.num_entities(CELL, ALL);

    std::vector<int> cell_num_mats;
    std::vector<int> cell_mat_ids;
    std::vector<double> cell_mat_volfracs;
    std::vector<Tangram::Point<D>> cell_mat_centroids;

    // Extract volume fraction and centroid data for cells in compact
    // cell-centric form (ccc)

    ccc_vfcen_data(cell_num_mats, cell_mat_ids, cell_mat_volfracs,
                   cell_mat_centroids);
          
    interface_reconstructor->set_volume_fractions(cell_num_mats,
                                                  cell_mat_ids,
                                                  cell_mat_volfracs,
                                                  cell_mat_centroids);
    interface_reconstructor->reconstruct();
  }

  
  // Make an intersector which knows about the source state (to be able
  // to query the number of materials, etc) and also knows about the
  // interface reconstructor so that it can retrieve pure material polygons


  Intersect<onwhat, Flat_Mesh_Wrapper<>, Flat_State_Wrapper<>,
            TargetMesh_Wrapper, InterfaceReconstructorType,
            Matpoly_Splitter, Matpoly_Clipper>
      intersect(source_mesh_flat, source_state_flat, target_mesh_,
                interface_reconstructor);

  // Get an instance of the desired interpolate algorithm type
  Interpolate<D, onwhat, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper,
              Flat_State_Wrapper<>, InterfaceReconstructorType,
              Matpoly_Splitter, Matpoly_Clipper>
      interpolate(source_mesh_flat, target_mesh_, source_state_flat, 
                  interface_reconstructor);
#else

  Intersect<onwhat, Flat_Mesh_Wrapper<>, Flat_State_Wrapper<>,
            TargetMesh_Wrapper, DummyInterfaceReconstructor,
            void, void>
      intersect(source_mesh_flat, source_state_flat, target_mesh_);

  // Get an instance of the desired interpolate algorithm type
  Interpolate<D, onwhat, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper,
              Flat_State_Wrapper<>, DummyInterfaceReconstructor,
              void, void>
      interpolate(source_mesh_flat, target_mesh_, source_state_flat);
#endif  // HAVE_TANGRAM

  
  //--------------------------------------------------------------------
  // REMAP MESH FIELDS FIRST (this requires just mesh-mesh intersection)
  //--------------------------------------------------------------------
  
  // INTERSECT
  
  gettimeofday(&begin_timeval, 0);
  
  // For each cell in the target mesh get a list of candidate-weight
  // pairings (in a traditional mesh, not particle mesh, the weights
  // are moments). Note that this candidate list is different from the
  // search candidate list in that it may not include some of the
  // search candidates. Also, note that for 2nd order and higher
  // remaps, we get multiple moments (0th, 1st, etc) for each
  // target-source cell intersection
  
  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(),
                     source_ents_and_weights.begin(),
                     intersect);
  
  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_xsect += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
  

  // INTERPOLATE (one variable at a time)

  gettimeofday(&begin_timeval, 0);

  int nvars = src_meshvar_names.size();
  if (comm_rank == 0)
    std::cout << "Number of mesh variables on entity kind " << onwhat <<
        " to remap is " << nvars << std::endl;

  for (int i = 0; i < nvars; ++i) {
    interpolate.set_interpolation_variable(src_meshvar_names[i], limiters_[i]);
    
    // Get a handle to a memory location where the target state
    // would like us to write this material variable into. If it is
    // NULL, we allocate it ourself
    
    double *target_field_raw;
    target_state_.mesh_get_data(onwhat, trg_meshvar_names[i], &target_field_raw);
    assert (target_field_raw != nullptr);
    
    
    Portage::pointer<double> target_field(target_field_raw);

    Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                       target_mesh_.end(onwhat, PARALLEL_OWNED),
                       source_ents_and_weights.begin(),
                       target_field, interpolate);
  }
  
  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;


  
  //--------------------------------------------------------------------
  // REMAP MULTIMATERIAL FIELDS NEXT, ONE MATERIAL AT A TIME
  //--------------------------------------------------------------------

  if (onwhat != CELL) return 1;
  
  // Material centric loop

  for (int m = 0; m < nmats; m++) {

    // INTERSECT
    
    gettimeofday(&begin_timeval, 0);

    intersect.set_material(m);

    // For each cell in the target mesh get a list of candidate-weight
    // pairings (in a traditional mesh, not particle mesh, the weights
    // are moments). Note that this candidate list is different from the
    // search candidate list in that it may not include some of the
    // search candidates. Also, note that for 2nd order and higher
    // remaps, we get multiple moments (0th, 1st, etc) for each
    // target-source cell intersection

    // NOTE: IDEALLY WE WOULD REUSE THE MESH-MESH INTERSECTIONS FROM THE
    // PREVIOUS STEP WHEN THE SOURCE MATERIAL CONTAINS ONLY ONE MATERIAL
    //
    // UNFORTUNATELY, THE REQUIREMENT OF THE INTERSECT FUNCTOR IS THAT
    // IT CANNOT MODIFY STATE, THIS MEANS WE CANNOT STORE THE MESH-MESH
    // INTERSECTION VALUES AND REUSE THEM AS NECESSARY FOR MESH-MATERIAL
    // INTERSECTION COMPUTATIONS
    
    Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                       target_mesh_.end(onwhat, PARALLEL_OWNED),
                       candidates.begin(),
                       source_ents_and_weights.begin(),
                       intersect);
    
    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    tot_seconds_xsect += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

    // LOOK AT INTERSECTION WEIGHTS TO DETERMINE WHICH TARGET CELLS
    // WILL GET NEW MATERIALS
    
    int ntargetcells = target_mesh_.num_entities(CELL, ALL);
    std::vector<int> matcellstgt;
    
    for (int c = 0; c < ntargetcells; c++) {
      std::vector<Weights_t> const& cell_sources_and_weights =
          source_ents_and_weights[c];
      for (int s = 0; s < cell_sources_and_weights.size(); s++) {
        std::vector<double> const& wts = cell_sources_and_weights[s].weights;
        if (wts[0] > 0.0) {
          double vol = target_mesh_.cell_volume(c);
          if (wts[0]/vol > 1.0e-10) {  // Check that the volume of material
                                       // we are adding to c is not miniscule
            matcellstgt.push_back(c);
            break;
          }
        }
      }
    }
    
    // If any processor is adding this material to the target state,
    // add it on all the processors
    
    int nmatcells = matcellstgt.size();
    int nmatcells_global = 0;
    #ifdef ENABLE_MPI
    MPI_Allreduce(&nmatcells, &nmatcells_global, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    #else
    nmatcells_global=nmatcells;
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
        // CELLS TO THE STATEMANAGER, IT ALSO MAKES SPACE FOR FIELD
        // VALUES FOR THIS MATERIAL IN EVERY MULTI-MATERIAL VECTOR IN
        // THE STATE MANAGER. THIS ENSURES THAT WHEN WE CALL
        // mat_get_celldata FOR A MATERIAL IN MULTI-MATERIAL STATE
        // VECTOR IT WILL ALREADY HAVE SPACE ALLOCATED FOR FIELD
        // VALUES OF THAT MATERIAL. SOME STATE WRAPPERS COULD CHOOSE
        // TO MAKE THIS A SIMPLER ROUTINE THAT ONLY STORES THE NAME
        // AND THE CELLS IN THE MATERIAL AND ACTUALLY ALLOCATE SPACE
        // FOR FIELD VALUES OF A MATERIAL IN A MULTI-MATERIAL FIELD
        // WHEN mat_get_celldata IS INVOKED.
        
        target_state_.add_material(source_state_.material_name(m), matcellstgt);
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
    std::vector<std::vector<Weights_t>> mat_sources_and_weights(nmatcells);
        
    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcellstgt[ic];
      double matvol = 0.0;
      Point<D> matcen;
      std::vector<Weights_t> const& cell_sources_and_weights =
          source_ents_and_weights[c];
      for (int s = 0; s < cell_sources_and_weights.size(); s++) {
        std::vector<double> const& wts = cell_sources_and_weights[s].weights;
        matvol += wts[0];
        for (int d = 0; d < D; d++)
          matcen[d] += wts[d+1];
      }
      matcen /= matvol;
      mat_volfracs[ic] = matvol/target_mesh_.cell_volume(c);
      mat_centroids[ic] = matcen;

      mat_sources_and_weights[ic] = cell_sources_and_weights;
    }

    target_state_.mat_add_celldata("mat_volfracs", m, &(mat_volfracs[0]));
    target_state_.mat_add_celldata("mat_centroids", m, &(mat_centroids[0]));


    // INTERPOLATE (one variable at a time)

    // HERE WE COULD MAKE A NEW LIST BASED ON WHICH TARGET CELLS HAVE ANY
    // INTERSECTIONS WITH SOURCE CELLS FOR THIS MATERIAL TO AVOID A NULL-OP
    // AND A WARNING MESSAGE ABOUT NO SOURCE CELLS CONTRIBUTING TO A TARGET -
    // IS IT WORTH IT?

    gettimeofday(&begin_timeval, 0);

    int nmatvars = src_matvar_names.size();
    if (comm_rank == 0)
      std::cout << "Number of multi-material variables on entity kind " <<
          onwhat << " to remap is " << nmatvars << std::endl;

    interpolate.set_material(m);    // We have to do this so we know
                                    // which material values we have
                                    // to grab from the source state

    for (int i = 0; i < nmatvars; ++i) {
      interpolate.set_interpolation_variable(src_matvar_names[i], limiters_[i]);

      // Get a handle to a memory location where the target state
      // would like us to write this material variable into. If it is
      // NULL, we allocate it ourself

      double *target_field_raw;
      target_state_.mat_get_celldata(trg_matvar_names[i], m, &target_field_raw);
      assert (target_field_raw != nullptr);


      Portage::pointer<double> target_field(target_field_raw);

      Portage::transform(matcellstgt.begin(), matcellstgt.end(),
                         mat_sources_and_weights.begin(),
                         target_field, interpolate);

      // If the state wrapper knows that the target data is already
      // laid out in this way and it gave us a pointer to the array
      // where the values reside, it has to do nothing in this
      // call. If the storage format is different, however, it may
      // have to copy the values into their proper locations

      target_state_.mat_add_celldata(trg_matvar_names[i], m, target_field_raw);
    }  // nmatvars

    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  }  // for nmats

  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

  std::cout << "Transform Time for Entity Kind " << onwhat << " on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;

  return 1;
}
#endif  // ENABLE_MPI


// Convert volume fraction and centroid data from compact
// material-centric to compact cell-centric (ccc) form as needed by
// Tangram

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
          template <class, int, class, class> class, 
          class, class> class Intersect,
          template<int, Entity_kind, class, class, class,
          template<class, int, class, class> class,
          class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper>
void
MMDriver<Search, Intersect, Interpolate, D,
         SourceMesh_Wrapper, SourceState_Wrapper,
         TargetMesh_Wrapper, TargetState_Wrapper,
         InterfaceReconstructorType, Matpoly_Splitter, 
         Matpoly_Clipper
         >::ccc_vfcen_data(std::vector<int>& cell_num_mats,
                           std::vector<int>& cell_mat_ids,
                           std::vector<double>& cell_mat_volfracs,
                           std::vector<Tangram::Point<D>>& cell_mat_centroids) {
  
  int nsourcecells = source_mesh_.num_entities(CELL, ALL);
  
  int nmats = source_state_.num_materials();
  cell_num_mats.assign(nsourcecells, 0);
  
  // First build full arrays (as if every cell had every material)
  
  std::vector<int> cell_mat_ids_full(nsourcecells*nmats, -1);
  std::vector<double> cell_mat_volfracs_full(nsourcecells*nmats, 0.0);
  std::vector<Tangram::Point<D>> cell_mat_centroids_full(nsourcecells*nmats);
  
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
    
    Portage::Point<D> const *matcenvec;
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


}  // namespace Portage

#endif  // PORTAGE_DRIVER_H_
