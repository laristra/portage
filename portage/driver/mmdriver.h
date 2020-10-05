/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_DRIVER_MMDRIVER_H_
#define PORTAGE_DRIVER_MMDRIVER_H_

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <memory>
#include <limits>
#include <cmath>

#include "wonton/support/wonton.h"
#include "wonton/mesh/flat/flat_mesh_wrapper.h"
#include "wonton/state/flat/flat_state_mm_wrapper.h"
#include "wonton/support/Point.h"
#include "wonton/state/state_vector_multi.h"

#ifdef PORTAGE_HAS_TANGRAM
  #include "tangram/driver/driver.h"
#endif

#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/support/portage.h"
#include "portage/support/timer.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/driver/coredriver.h"

#ifdef WONTON_ENABLE_MPI
  #include "portage/distributed/mpi_bounding_boxes.h"
#endif

/*!
  @file mmdriver.h
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
          template<int, Entity_kind, class, class, class, class, class,
                   template<class, int, class, class> class,
                   class, class, class=Wonton::DefaultCoordSys
                   > class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper = SourceMesh_Wrapper,
          class TargetState_Wrapper = SourceState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType = DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void>
class MMDriver {

#ifdef PORTAGE_HAS_TANGRAM
  // alias for interface reconstructor parameterized on the mesh type.
  // it will be used for gradient field computation.
  template<typename SourceMesh>
  using InterfaceReconstructor = Tangram::Driver<InterfaceReconstructorType, D,
                                                 SourceMesh, Matpoly_Splitter,
                                                 Matpoly_Clipper>;
#endif

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
      : source_mesh_(sourceMesh),
        target_mesh_(targetMesh),
        source_state_(sourceState),
        target_state_(targetState),
        dim_(sourceMesh.space_dimension()) {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
  }

  /// Copy constructor (disabled)
  MMDriver(const MMDriver &) = delete;

  /// Assignment operator (disabled)
  MMDriver & operator = (const MMDriver &) = delete;

  /// Destructor
  ~MMDriver() = default;

  /// Enable move semantics
  MMDriver(MMDriver &&) noexcept = default;

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source mesh to the target mesh.  This variable must
    exist in both meshes' state manager
  */
  void set_remap_var_names(std::vector<std::string> const &remap_var_names) {
    // remap variable names same in source and target mesh
    set_remap_var_names(remap_var_names, remap_var_names);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target mesh.
  */

  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names) {

    assert(source_remap_var_names.size() == target_remap_var_names.size());

    // No appending allowed
    source_target_varname_map_.clear();

    int nvars = source_remap_var_names.size();
#ifndef NDEBUG
    for (int i = 0; i < nvars; ++i) {
      Entity_kind srckind = source_state_.get_entity(source_remap_var_names[i]);
      Entity_kind trgkind = target_state_.get_entity(target_remap_var_names[i]);
      if (trgkind == Entity_kind::UNKNOWN_KIND)
        continue;  // Presumably field does not exist on target - will get added

      assert(srckind == trgkind);  // if target field exists, entity kinds must match
    }
#endif

    for (int i = 0; i < nvars; i++) {
      source_target_varname_map_[source_remap_var_names[i]] = target_remap_var_names[i];

      // Set options so that defaults will produce something reasonable
      limiters_[source_remap_var_names[i]] = Limiter_type::BARTH_JESPERSEN;
      bnd_limiters_[source_remap_var_names[i]] = Boundary_Limiter_type::BND_NOLIMITER;
      partial_fixup_types_[target_remap_var_names[i]] =
          Partial_fixup_type::GLOBALLY_CONSERVATIVE;
      empty_fixup_types_[target_remap_var_names[i]] =
          Empty_fixup_type::EXTRAPOLATE;
    }
  }

  /*!
    @brief set limiter for all variables
    @param limiter  Limiter to use for second order reconstruction (NOLIMITER or BARTH_JESPERSEN)
  */
  void set_limiter(Limiter_type limiter) {
    for (auto const& stpair : source_target_varname_map_) {
      std::string const& source_var_name = stpair.first;
      limiters_[source_var_name] = limiter;
    }
  }

  /*!
    @brief set boundary limiter for all variables
    @param bnd_limiter  Boundary limiter to use for second order reconstruction (BND_NOLIMITER
                        BND_ZERO_GRADIENT, or BND_BARTH_JESPERSEN))
  */
  void set_bnd_limiter(Boundary_Limiter_type bnd_limiter) {
    for (auto const& stpair : source_target_varname_map_) {
      std::string const& source_var_name = stpair.first;
      bnd_limiters_[source_var_name] = bnd_limiter;
    }
  }  

  /*!
    @brief set limiter for a variable
    @param target_var_name Source mesh variable whose gradient is to be limited
    @param limiter  Limiter to use for second order reconstruction (NOLIMITER
    or BARTH_JESPERSEN)
  */
  void set_limiter(std::string const& source_var_name, Limiter_type limiter) {
    limiters_[source_var_name] = limiter;
  }

  /*!
    @brief set boundary limiter for a variable
    @param target_var_name Source mesh variable whose gradient is to be limited
    on the boundary
    @param bnd_limiter  Boundary limiter to use for second order reconstruction (BND_NOLIMITER,
                        BND_ZERO_GRADIENT, or BND_BARTH_JESPERSEN))
  */
  void set_bnd_limiter(std::string const& source_var_name, Boundary_Limiter_type bnd_limiter) {
    bnd_limiters_[source_var_name] = bnd_limiter;
  }

  /*!
    @brief set flag whether we want to check for mesh mismatch

    @param do_check_mismatch  boolean flag indicating if we want to check for
                              boundary mismatch between meshes

    This check is used to determine if the boundaries of the two
    meshes overlap exactly. If they don't conservation is
    violated. Callers can ask to compensate for the mismatch when
    interpolating a mesh variable.
  */
  void set_check_mismatch_flag(const bool do_check_mismatch) {
    do_check_mismatch_ = do_check_mismatch;
  }
  
  /*!
    @brief set repair method in partially filled cells for all variables
    @param fixup_type Can be Partial_fixup_type::CONSTANT,
                      Partial_fixup_type::LOCALLY_CONSERVATIVE,
                      Partial_fixup_type::GLOBALLY_CONSERVATIVE
  */
  void set_partial_fixup_type(Partial_fixup_type fixup_type) {
    for (auto const& stpair : source_target_varname_map_) {
      std::string const& target_var_name = stpair.second;
      partial_fixup_types_[target_var_name] = fixup_type;
    }
  }

  /*!
    @brief set repair method in partially filled cells for all variables
    @param target_var_name Target mesh variable to set fixup option for
    @param fixup_type  Can be Partial_fixup_type::CONSTANT,
                       Partial_fixup_type::LOCALLY_CONSERVATIVE,
                       Partial_fixup_type::GLOBALLY_CONSERVATIVE
  */
  void set_partial_fixup_type(std::string const& target_var_name,
                              Partial_fixup_type fixup_type) {
    partial_fixup_types_[target_var_name] = fixup_type;
  }

  /*!
    @brief set repair method in empty cells for all variables
    @param fixup_type Can be Empty_fixup_type::LEAVE_EMPTY,
                      Empty_fixup_type::EXTRAPOLATE
  */
  void set_empty_fixup_type(Empty_fixup_type fixup_type) {
    for (auto const& stpair : source_target_varname_map_) {
      std::string const& target_var_name = stpair.second;
      empty_fixup_types_[target_var_name] = fixup_type;
    }
  }

  /*!
    @brief set repair method in empty cells for all variables
    @param target_var_name Target mesh variable to set fixup option for
    @param fixup_type Can be Empty_fixup_type::LEAVE_EMPTY,
    Empty_fixup_type::EXTRAPOLATE
  */
  void set_empty_fixup_type(std::string const& target_var_name,
                            Empty_fixup_type fixup_type) {
    empty_fixup_types_[target_var_name] = fixup_type;
  }


  void set_max_fixup_iter(int maxiter) {
    max_fixup_iter_ = maxiter;
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

  /*!
    @brief set the bounds of variable to be remapped on target
    @param target_var_name Name of variable in target mesh to limit
  */
  template<typename T>
  void set_remap_var_bounds(std::string target_var_name,
                            T lower_bound, T upper_bound) {
    if (typeid(T) == typeid(double)) {
      double_lower_bounds_[target_var_name] = lower_bound;
      double_upper_bounds_[target_var_name] = upper_bound;
    } else
      throw std::runtime_error("Remap variable type not supported");
  }

  /*!
    @brief set conservation tolerance of variable to be remapped on target
    @param target_var_name Name of variable in target mesh to limit
  */
  template<typename T>
  void set_conservation_tolerance(std::string target_var_name,
                                  T conservation_tol) {
    if (typeid(T) == typeid(double))
      conservation_tol_[target_var_name] = conservation_tol;
    else
      throw std::runtime_error("Remap variable type not supported");
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
  void set_reconstructor_options(bool all_convex,
                                 const std::vector<Tangram::IterativeMethodTolerances_t> &tols = 
                                   std::vector<Tangram::IterativeMethodTolerances_t>()){
    reconstructor_tols_ = tols; 
    reconstructor_all_convex_ = all_convex; 
  }

#endif

  /*!
    @brief Get the names of the variables to be remapped from the
    source mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> source_remap_var_names() const {
    std::vector<std::string> source_var_names;
    source_var_names.reserve(source_target_varname_map_.size());
    for (auto const& stname_pair : source_target_varname_map_)
      source_var_names.push_back(stname_pair.first);
    return source_var_names;
  }

  /*!
    @brief Get the names of the variables to be remapped to the
    target mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> target_remap_var_names() const {
    std::vector<std::string> target_var_names;
    target_var_names.reserve(source_target_varname_map_.size());
    for (auto const& stname_pair : source_target_varname_map_)
      target_var_names.push_back(stname_pair.second);
    return target_var_names;
  }

  /*!
    @brief Get the dimensionality of the meshes.
    @return The dimensionality of the meshes.
  */
  unsigned int dim() const {
    return dim_;
  }


  /*!
    @brief remap for a given set of MESH and MATERIAL variables on CELLS
    @tparam SourceMesh_Wrapper2  May be the mesh wrapper sent into MMDriver or the Flat_Mesh_Wrapper created for redistribution
    @tparam SourceState_Wrapper2 May be the state wrapper sent into MMDriver or the Flat_State_Wrapper created for redistribution
    @tparam entity_kind  Kind of entity that variables live on

    @param source_meshvar_names  names of remap variables on source mesh
    @param target_meshvar_names  names of remap variables on target mesh
    @param source_matvar_names  names of remap variables on materials of source mesh
    @param target_matvar_names  names of remap variables on materials of target mesh
    @param executor             pointer to Serial Executor (generally not needed but introduced for future proofing)
    @return status of remap (1 if successful, 0 if not)
  */

  template<class SourceMesh_Wrapper2, class SourceState_Wrapper2>
  int cell_remap(SourceMesh_Wrapper2 const & source_mesh2,
                 SourceState_Wrapper2 const & source_state2,
                 std::vector<std::string> const &src_meshvar_names,
                 std::vector<std::string> const &trg_meshvar_names,
                 std::vector<std::string> const &src_matvar_names,
                 std::vector<std::string> const &trg_matvar_names,
                 Wonton::Executor_type const *executor = nullptr);


  /*!
    @brief remap for a given set of MESH variables on NODES
    @tparam SourceMesh_Wrapper2  May be the mesh wrapper sent into MMDriver or the Flat_Mesh_Wrapper created for redistribution
    @tparam SourceState_Wrapper2 May be the state wrapper sent into MMDriver or the Flat_State_Wrapper created for redistribution
    @tparam entity_kind  Kind of entity that variables live on

    @param source_meshvar_names  names of remap variables on source mesh
    @param target_meshvar_names  names of remap variables on target mesh
    @param executor             pointer to Serial Executor (generally not needed but introduced for future proofing)
    @return status of remap (1 if successful, 0 if not)
  */

  template<class SourceMesh_Wrapper2, class SourceState_Wrapper2>
  int node_remap(SourceMesh_Wrapper2 const & source_mesh2,
                 SourceState_Wrapper2 const & source_state2,
                 std::vector<std::string> const &src_meshvar_names,
                 std::vector<std::string> const &trg_meshvar_names,
                 Wonton::Executor_type const *executor = nullptr);



  /*!

    @brief Compute mismatch fixup bounds dynamically
    @tparam SourceState_Wrapper2 May be the state wrapper sent into MMDriver or the Flat_State_Wrapper created for redistribution
    @tparam entity_kind  Kind of entity that variables live on

    @param source_state2         source state wrapper
    @param source_meshvar_names  names of remap variables on source mesh
    @param target_meshvar_names  names of remap variables on target mesh
    @param sources_and_weights   mesh-mesh intersection weights
    @param executor              pointer to Executor
    
  */
  
  template<class SourceState_Wrapper2,Entity_kind onwhat>
  void compute_bounds(SourceState_Wrapper2 const& source_state2,
                   std::vector<std::string> const& src_meshvar_names,
                   std::vector<std::string> const& trg_meshvar_names,
                   Wonton::Executor_type const *executor = nullptr) {

#ifdef WONTON_ENABLE_MPI
    MPI_Comm mycomm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm = mpiexecutor->mpicomm;
    }
#endif

    int const nvars = src_meshvar_names.size();
      
    // loop over mesh variables
    for (int i = 0; i < nvars; i++) {
    
      auto& src_var = src_meshvar_names[i];
      auto& trg_var = trg_meshvar_names[i];

      // See if we have caller specified bounds, if so, keep them. 
      // Otherwise, determine sensible values
      if (not double_lower_bounds_.count(trg_var) or
          not double_upper_bounds_.count(trg_var)) {

        // Since caller has not specified bounds for variable, attempt
        // to derive them from source state. This code should go into
        // Wonton into each state manager (or better yet, deprecated :-))
        
        int nsrcents = source_mesh_.num_entities(onwhat,
                                                 Entity_type::PARALLEL_OWNED);
        
        double const *source_data;
        source_state_.mesh_get_data(onwhat, src_var, &source_data);
        double lower_bound = *std::min_element(source_data, source_data + nsrcents);
        double upper_bound = *std::max_element(source_data, source_data + nsrcents);
        
#ifdef WONTON_ENABLE_MPI
        if (mycomm != MPI_COMM_NULL) {
          double global_bounds[2] = { 0.0, 0.0 };
          MPI_Allreduce(&lower_bound, global_bounds+0, 1, MPI_DOUBLE, MPI_MIN, mycomm);
          MPI_Allreduce(&upper_bound, global_bounds+1, 1, MPI_DOUBLE, MPI_MAX, mycomm);
          lower_bound = global_bounds[0];
          upper_bound = global_bounds[1];
        }
#endif

        double relbounddiff = std::abs((upper_bound-lower_bound)/lower_bound);
        if (relbounddiff < consttol_) {
          // The field is constant over the source mesh/part. We HAVE to
          // relax the bounds to be able to conserve the integral quantity
          // AND maintain a constant.
          lower_bound -= 0.5*lower_bound;
          upper_bound += 0.5*upper_bound;
        }
        
        // set the bounds on the instance variable
        set_remap_var_bounds(trg_var, lower_bound, upper_bound);      
      }
        
      // see if caller has specified a tolerance for conservation
      if (not conservation_tol_.count(trg_var)) {
        set_conservation_tolerance(trg_var, DEFAULT_NUMERIC_TOLERANCES<D>.relative_conservation_eps);     
      }      
    }
  }  // compute_bounds



  /*!
    @brief Execute the remapping process
    @return status of remap (1 if successful, 0 if not)
  */
  int run(Wonton::Executor_type const *executor = nullptr,
          std::string *errmsg = nullptr) {
    std::string message;

#ifndef NDEBUG
    auto tic = timer::now();

    int comm_rank = 0;
#endif

#ifdef WONTON_ENABLE_MPI
    bool distributed = false;

    MPI_Comm mycomm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm = mpiexecutor->mpicomm;
#ifndef NDEBUG      
      MPI_Comm_rank(mycomm, &comm_rank);
#endif
      int nprocs = 0;
      MPI_Comm_size(mycomm, &nprocs);
      if (nprocs > 1)
        distributed = true;
    }
#endif

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
    if (comm_rank == 0)
      std::cout << "in MMDriver::run()...\n";

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;
#endif

    std::vector<std::string> src_meshvar_names, src_matvar_names;
    std::vector<std::string> trg_meshvar_names, trg_matvar_names;


    // -------- CELL VARIABLE REMAP ---------
    // Collect all cell based variables and remap them

    for (auto const& stpair : source_target_varname_map_) {
      std::string const& srcvarname = stpair.first;
      Entity_kind onwhat = source_state_.get_entity(srcvarname);
      if (onwhat == CELL) {
        // Separate out mesh fields and multi-material fields - they will be
        // processed differently

        std::string const& trgvarname = stpair.second;

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

    // Default is serial run (if MPI is not enabled or the
    // communicator is not defined or the number of processors is 1)
#ifdef WONTON_ENABLE_MPI
    // Create a new mesh wrapper that we can use for redistribution
    // of the source mesh as necessary (so that every target cell
    // sees any source cell that it overlaps with)
    
    Flat_Mesh_Wrapper<> source_mesh_flat;
    Flat_State_Wrapper<Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);

    bool redistributed_source = false;
    if (distributed) {
      MPI_Bounding_Boxes distributor(mpiexecutor);
      if (distributor.is_redistribution_needed(source_mesh_, target_mesh_)) {
#ifndef NDEBUG
        tic = timer::now();
#endif
        
        source_mesh_flat.initialize(source_mesh_);
        
        // Note the flat state should be used for everything including the
        // centroids and volume fractions for interface reconstruction
        std::vector<std::string> source_remap_var_names;
        for (auto & stpair : source_target_varname_map_)
          source_remap_var_names.push_back(stpair.first);
        source_state_flat.initialize(source_state_, source_remap_var_names);
        
        distributor.distribute(source_mesh_flat, source_state_flat,
                               target_mesh_, target_state_);
        
        redistributed_source = true;
        
#if !defined(NDEBUG)
        float tot_seconds_dist = timer::elapsed(tic);
        std::cout << "Redistribution Time Rank " << comm_rank << " (s): " <<
            tot_seconds_dist << std::endl;
#endif
      }
    }

    if (redistributed_source) {
      
      // Why is it not able to deduce the template arguments, if I don't specify
      // Flat_Mesh_Wrapper and Flat_State_Wrapper?
        
      cell_remap<Flat_Mesh_Wrapper<>, Flat_State_Wrapper<Flat_Mesh_Wrapper<>>>
          (source_mesh_flat, source_state_flat,
           src_meshvar_names, trg_meshvar_names,
           src_matvar_names,  trg_matvar_names,
           executor);
    }
    else
#endif
    {
      // Why is it not able to deduce the template arguments, if I don't specify
      // Source_Mesh_Wrapper and Source_State_Wrapper?

      cell_remap<SourceMesh_Wrapper, SourceState_Wrapper>
          (source_mesh_, source_state_,
           src_meshvar_names, trg_meshvar_names,
           src_matvar_names, trg_matvar_names,
           executor);
    }



    // -------- NODE VARIABLE REMAP ---------
    // Collect all node based variables and remap them
    // (ignore any multi-material variables on NODES - not well defined)

    src_meshvar_names.clear(); src_matvar_names.clear();
    trg_meshvar_names.clear(); trg_matvar_names.clear();

    for (auto const& stpair : source_target_varname_map_) {
      std::string const& srcvarname = stpair.first;
      Entity_kind onwhat = source_state_.get_entity(srcvarname);
      if (onwhat == NODE) {
        std::string const& trgvarname = stpair.second;

        Field_type ftype = source_state_.field_type(onwhat, srcvarname);

        if (ftype == Field_type::MESH_FIELD) {
          src_meshvar_names.push_back(srcvarname);
          trg_meshvar_names.push_back(trgvarname);
        } else if (ftype == Field_type::MULTIMATERIAL_FIELD)
          throw std::runtime_error("Cannot handle multi-material fields on nodes\n");
      }
    }

    if (not src_meshvar_names.empty()) {
#ifdef WONTON_ENABLE_MPI
      if (redistributed_source)
        node_remap<Flat_Mesh_Wrapper<>, Flat_State_Wrapper<Flat_Mesh_Wrapper<>>>
            (source_mesh_flat, source_state_flat,
             src_meshvar_names, trg_meshvar_names,
             executor);
      else
#endif
        node_remap<SourceMesh_Wrapper, SourceState_Wrapper>
            (source_mesh_, source_state_,
             src_meshvar_names, trg_meshvar_names,
             executor);
    }

    return 1;
  }  // run



 private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::unordered_map<std::string, std::string> source_target_varname_map_;
  std::unordered_map<std::string, Limiter_type> limiters_;
  std::unordered_map<std::string, Boundary_Limiter_type> bnd_limiters_;
  std::unordered_map<std::string, Partial_fixup_type> partial_fixup_types_;
  std::unordered_map<std::string, Empty_fixup_type> empty_fixup_types_;
  std::unordered_map<std::string, double> double_lower_bounds_;
  std::unordered_map<std::string, double> double_upper_bounds_;
  std::unordered_map<std::string, double> conservation_tol_;
  unsigned int dim_;
  double consttol_ =  100*std::numeric_limits<double>::epsilon();
  int max_fixup_iter_ = 5;
  NumericTolerances_t num_tols_ = DEFAULT_NUMERIC_TOLERANCES<D>;
  bool do_check_mismatch_ = true;


#ifdef PORTAGE_HAS_TANGRAM
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
  // There is an associated method called set_reconstructor_options
  // that should be invoked to set user-specific values. Otherwise,
  // the remapper will use the default values.
  std::vector<Tangram::IterativeMethodTolerances_t> reconstructor_tols_; 
  bool reconstructor_all_convex_ = true;  
#endif
  
};  // class MMDriver



// remap routine specialization for cells

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
                    template <class, int, class, class> class,
                    class, class> class Intersect,
          template<int, Entity_kind, class, class, class, class, class,
                   template<class, int, class, class> class,
                   class, class, class=Wonton::DefaultCoordSys>
          class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper>
template<class SourceMesh_Wrapper2, class SourceState_Wrapper2>
int MMDriver<Search, Intersect, Interpolate, D,
             SourceMesh_Wrapper, SourceState_Wrapper,
             TargetMesh_Wrapper, TargetState_Wrapper,
             InterfaceReconstructorType, Matpoly_Splitter,
             Matpoly_Clipper
             >::cell_remap(SourceMesh_Wrapper2 const & source_mesh2,
                           SourceState_Wrapper2 const & source_state2,
                           std::vector<std::string> const &src_meshvar_names,
                           std::vector<std::string> const &trg_meshvar_names,
                           std::vector<std::string> const &src_matvar_names,
                           std::vector<std::string> const &trg_matvar_names,
                           Wonton::Executor_type const *executor) {

#ifndef NDEBUG
  int comm_rank = 0;
#ifdef WONTON_ENABLE_MPI
  MPI_Comm mycomm = MPI_COMM_NULL;
  auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
  if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
    mycomm = mpiexecutor->mpicomm;
    MPI_Comm_rank(mycomm, &comm_rank);
  }
#endif
#endif

#ifndef NDEBUG
  float tot_seconds = 0.0;
  float tot_seconds_srch = 0.0;
  float tot_seconds_xsect = 0.0;
  float tot_seconds_interp = 0.0;
  auto tic = timer::now();
#endif

  std::vector<std::string> source_remap_var_names;
  for (auto & stpair : source_target_varname_map_)
    source_remap_var_names.push_back(stpair.first);

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

  // Instantiate core driver

  Portage::CoreDriver<D, CELL,
                      SourceMesh_Wrapper2, SourceState_Wrapper2,
                      TargetMesh_Wrapper, TargetState_Wrapper,
                      InterfaceReconstructorType,
                      Matpoly_Splitter, Matpoly_Clipper>
      coredriver_cell(source_mesh2, source_state2, target_mesh_, target_state_, executor);

  coredriver_cell.set_num_tols(num_tols_);
#ifdef PORTAGE_HAS_TANGRAM
  coredriver_cell.set_interface_reconstructor_options(reconstructor_all_convex_,
                                                      reconstructor_tols_);
#endif  
  
  // SEARCH
  auto candidates = coredriver_cell.template search<Portage::SearchKDTree>();
#ifndef NDEBUG
  tot_seconds_srch = timer::elapsed(tic, true);
#endif
#ifdef PORTAGE_HAS_TANGRAM
  int nmats = source_state2.num_materials();
#endif

  //--------------------------------------------------------------------
  // REMAP MESH FIELDS FIRST (this requires just mesh-mesh intersection)
  //--------------------------------------------------------------------

  // INTERSECT MESHES
  auto source_ents_and_weights =
      coredriver_cell.template intersect_meshes<Intersect>(candidates);
#ifndef NDEBUG
  tot_seconds_xsect += timer::elapsed(tic);
#endif
  // check for mesh mismatch
  if (do_check_mismatch_)
    coredriver_cell.check_mismatch(source_ents_and_weights);

  // compute bounds (for all variables) if required for mismatch
  if (do_check_mismatch_ && coredriver_cell.has_mismatch())
    compute_bounds<SourceState_Wrapper2, CELL>
        (source_state2, src_meshvar_names, trg_meshvar_names, executor);
  
  // INTERPOLATE (one variable at a time)
  int nvars = src_meshvar_names.size();
#if !defined(NDEBUG)
  tic = timer::now();

#if defined(VERBOSE_OUTPUT)
  if (comm_rank == 0) {
    std::cout << "Number of mesh variables on cells to remap is " <<
        nvars << std::endl;
  }
#endif
#endif

  Wonton::vector<Vector<D>> gradients;
  // to check interpolation order
  using Interpolator = Interpolate<D, CELL,
                                   SourceMesh_Wrapper2, TargetMesh_Wrapper,
                                   SourceState_Wrapper2, TargetState_Wrapper,
                                   double, InterfaceReconstructorType,
                                   Matpoly_Splitter, Matpoly_Clipper>;


  for (int i = 0; i < nvars; ++i) {
    std::string const& srcvar = src_meshvar_names[i];
    std::string const& trgvar = trg_meshvar_names[i];

    if (Interpolator::order == 2) {
      // set slope limiters
      auto limiter = (limiters_.count(srcvar) ? limiters_[srcvar] : DEFAULT_LIMITER);
      auto bndlimit = (bnd_limiters_.count(srcvar) ? bnd_limiters_[srcvar] : DEFAULT_BND_LIMITER);
      // compute gradient field
      gradients = coredriver_cell.compute_source_gradient(srcvar, limiter, bndlimit);
      // interpolate
      coredriver_cell.template interpolate_mesh_var<double, Interpolate>
        (srcvar, trgvar, source_ents_and_weights, &gradients);
    } else /* order 1 */ {
      // just interpolate
      coredriver_cell.template interpolate_mesh_var<double, Interpolate>
        (srcvar, trgvar, source_ents_and_weights);
    }

    // fix mismatch if necessary
    if (do_check_mismatch_ && coredriver_cell.has_mismatch()) {
      coredriver_cell.fix_mismatch(srcvar, trgvar,
                                   double_lower_bounds_[trgvar],
                                   double_upper_bounds_[trgvar],
                                   conservation_tol_[trgvar],
                                   max_fixup_iter_,
                                   partial_fixup_types_[trgvar],
                                   empty_fixup_types_[trgvar]
      );
    }
  }

#ifndef NDEBUG
  tot_seconds_interp += timer::elapsed(tic, true);
#endif

#ifdef PORTAGE_HAS_TANGRAM
  if (nmats > 1) {
    //--------------------------------------------------------------------
    // REMAP MULTIMATERIAL FIELDS NEXT, ONE MATERIAL AT A TIME
    //--------------------------------------------------------------------
    
    auto source_ents_and_weights_mat =
        coredriver_cell.template intersect_materials<Intersect>(candidates);

    if (Interpolator::order == 2) { coredriver_cell.cache_multimat_gradient_stencils(); }
    
    int nmatvars = src_matvar_names.size();
    std::vector<Wonton::vector<Vector<D>>> matgradients(nmats);

    for (int i = 0; i < nmatvars; ++i) {
      std::string const& srcvar = src_matvar_names[i];
      std::string const& trgvar = trg_matvar_names[i];

      if (Interpolator::order == 2) {
        // set slope limiters
        auto limiter = (limiters_.count(srcvar) ? limiters_[srcvar] : DEFAULT_LIMITER);
        auto bndlimit = (bnd_limiters_.count(srcvar) ? bnd_limiters_[srcvar] : DEFAULT_BND_LIMITER);
        // compute gradient field for each material
        for (int m = 0; m < nmats; m++) {
          matgradients[m] =
            coredriver_cell.compute_source_gradient(src_matvar_names[i],
                                                    limiter, bndlimit, m);
        }
        // interpolate
        coredriver_cell.template interpolate_mat_var<double, Interpolate>
          (srcvar, trgvar, source_ents_and_weights_mat, &matgradients);
      } else {
        // interpolate
        coredriver_cell.template interpolate_mat_var<double, Interpolate>
          (srcvar, trgvar, source_ents_and_weights_mat);
      }
    }  // nmatvars
  }
#endif


#if !defined(NDEBUG)
  tot_seconds_interp += timer::elapsed(tic);
  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

  std::cout << "Time for Cell remap on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;

#if defined(VERBOSE_OUTPUT)
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;
#endif
#endif
  return 1;
}  // remap specialization for cells




// remap routine specialization for nodes

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
                    template <class, int, class, class> class,
                    class, class> class Intersect,
          template<int, Entity_kind, class, class, class, class, class,
                   template<class, int, class, class> class,
                   class, class, class=Wonton::DefaultCoordSys>
          class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter,
          class Matpoly_Clipper>
template<class SourceMesh_Wrapper2, class SourceState_Wrapper2>
int MMDriver<Search, Intersect, Interpolate, D,
             SourceMesh_Wrapper, SourceState_Wrapper,
             TargetMesh_Wrapper, TargetState_Wrapper,
             InterfaceReconstructorType, Matpoly_Splitter,
             Matpoly_Clipper
             >::node_remap(SourceMesh_Wrapper2 const & source_mesh2,
                           SourceState_Wrapper2 const & source_state2,
                           std::vector<std::string> const &src_meshvar_names,
                           std::vector<std::string> const &trg_meshvar_names,
                           Wonton::Executor_type const *executor) {

#ifndef NDEBUG
  int comm_rank = 0;

#ifdef WONTON_ENABLE_MPI
  MPI_Comm mycomm = MPI_COMM_NULL;
  auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
  if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
    mycomm = mpiexecutor->mpicomm;
    MPI_Comm_rank(mycomm, &comm_rank);
  }
#endif
#endif

#ifndef NDEBUG
  float tot_seconds = 0.0;
  float tot_seconds_srch = 0.0;
  float tot_seconds_xsect = 0.0;
  float tot_seconds_interp = 0.0;
  auto tic = timer::now();
#endif

  std::vector<std::string> source_remap_var_names;
  for (auto & stpair : source_target_varname_map_)
    source_remap_var_names.push_back(stpair.first);

  
#ifdef PORTAGE_HAS_TANGRAM
    // If user set tolerances for Tangram, but not for Portage,
    // use Tangram tolerances
    if (!num_tols_.user_tolerances && (!reconstructor_tols_.empty()) ) {
      num_tols_.min_absolute_distance = reconstructor_tols_[0].arg_eps;
      num_tols_.min_absolute_volume = reconstructor_tols_[0].fun_eps;
    }
    // If user did NOT set tolerances for Tangram, use Portage tolerances
    if (reconstructor_tols_.empty()) {
      reconstructor_tols_ = { {1000, num_tols_.min_absolute_distance,
                                     num_tols_.min_absolute_volume},
                              {100, num_tols_.min_absolute_distance,
                                    num_tols_.min_absolute_distance} };
    }
#endif
  // Instantiate core driver

  Portage::CoreDriver<D, NODE,
                      SourceMesh_Wrapper2, SourceState_Wrapper2,
                      TargetMesh_Wrapper, TargetState_Wrapper>
      coredriver_node(source_mesh2, source_state2, target_mesh_, target_state_, executor);

  coredriver_node.set_num_tols(num_tols_);
#ifdef PORTAGE_HAS_TANGRAM
  coredriver_node.set_interface_reconstructor_options(reconstructor_all_convex_,
                                                      reconstructor_tols_);
#endif  
  
  // SEARCH

  auto candidates = coredriver_node.template search<Portage::SearchKDTree>();
#ifndef NDEBUG
  tot_seconds_srch = timer::elapsed(tic, true);
#endif
  //--------------------------------------------------------------------
  // REMAP MESH FIELDS FIRST (this requires just mesh-mesh intersection)
  //--------------------------------------------------------------------

  // INTERSECT MESHES
  auto source_ents_and_weights =
      coredriver_node.template intersect_meshes<Intersect>(candidates);
#ifndef NDEBUG
  tot_seconds_xsect += timer::elapsed(tic);
#endif
  // check for mesh mismatch
  if (do_check_mismatch_)
    coredriver_node.check_mismatch(source_ents_and_weights);

  // compute bounds if required for mismatch
  if (do_check_mismatch_ && coredriver_node.has_mismatch())
    compute_bounds<SourceState_Wrapper2, NODE>
        (source_state2, src_meshvar_names, trg_meshvar_names, executor);

  // INTERPOLATE (one variable at a time)
  int nvars = src_meshvar_names.size();
#if !defined(NDEBUG)
  tic = timer::now();

#if defined(VERBOSE_OUTPUT)
  if (comm_rank == 0) {
    std::cout << "Number of mesh variables on nodes to remap is " <<
        nvars << std::endl;
  }
#endif
#endif

  Wonton::vector<Vector<D>> gradients;
  // to check interpolation order
  using Interpolator = Interpolate<D, NODE,
                                   SourceMesh_Wrapper2, TargetMesh_Wrapper,
                                   SourceState_Wrapper2, TargetState_Wrapper,
                                   double, InterfaceReconstructorType,
                                   Matpoly_Splitter, Matpoly_Clipper>;

  for (int i = 0; i < nvars; ++i) {
    std::string const& srcvar = src_meshvar_names[i];
    std::string const& trgvar = trg_meshvar_names[i];

    if (Interpolator::order == 2) {
      // set slope limiters
      auto limiter = (limiters_.count(srcvar) ? limiters_[srcvar] : DEFAULT_LIMITER);
      auto bndlimit = (bnd_limiters_.count(srcvar) ? bnd_limiters_[srcvar] : DEFAULT_BND_LIMITER);
      // compute gradient field
      gradients = coredriver_node.compute_source_gradient(srcvar, limiter, bndlimit);
      // interpolate
      coredriver_node.template interpolate_mesh_var<double, Interpolate>
        (srcvar, trgvar, source_ents_and_weights, &gradients);
    } else /* order 1 */ {
      // just interpolate
      coredriver_node.template interpolate_mesh_var<double, Interpolate>
        (srcvar, trgvar, source_ents_and_weights);
    }

    // fix mismatch if necessary
    if (do_check_mismatch_ && coredriver_node.has_mismatch()) {
      coredriver_node.fix_mismatch(srcvar, trgvar,
                                   double_lower_bounds_[trgvar],
                                   double_upper_bounds_[trgvar],
                                   conservation_tol_[trgvar],
                                   max_fixup_iter_,
                                   partial_fixup_types_[trgvar],
                                   empty_fixup_types_[trgvar]
      );
    }
  }

#if !defined(NDEBUG)
  tot_seconds_interp += timer::elapsed(tic);
  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

  std::cout << "Time for Node remap on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;

#if defined(VERBOSE_OUTPUT)
     std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;
#endif
#endif
  return 1;
}  // remap specialization for nodes


}  // namespace Portage

#endif  // PORTAGE_DRIVER_MMDRIVER_H_
