/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_DRIVER_MMDRIVER_H_
#define PORTAGE_DRIVER_MMDRIVER_H_

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
#include <cmath>

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"
#endif

#include "portage/intersect/dummy_interface_reconstructor.h"

#include "portage/support/portage.h"

#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "wonton/mesh/flat/flat_mesh_wrapper.h"
#include "wonton/state/flat/flat_state_mm_wrapper.h"
#include "wonton/support/Point.h"
#include "wonton/state/state_vector_multi.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/driver/coredriver.h"

#ifdef PORTAGE_ENABLE_MPI
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

#if HAVE_TANGRAM
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

  /// Enable move semantics
  MMDriver(MMDriver &&) = default;

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
    for (int i = 0; i < nvars; ++i) {
      Entity_kind srckind = source_state_.get_entity(source_remap_var_names[i]);
      Entity_kind trgkind = target_state_.get_entity(target_remap_var_names[i]);
      if (trgkind == Entity_kind::UNKNOWN_KIND)
        continue;  // Presumably field does not exist on target - will get added

      assert(srckind == trgkind);  // if target field exists, entity kinds
                                   // must match
    }

    for (int i = 0; i < nvars; i++) {
      source_target_varname_map_[source_remap_var_names[i]] = target_remap_var_names[i];

      // Set options so that defaults will produce something reasonable
      limiters_[source_remap_var_names[i]] = Limiter_type::BARTH_JESPERSEN;
      bnd_limiters_[source_remap_var_names[i]] = Boundary_Limiter_type::BND_NOLIMITER;
      partial_fixup_types_[target_remap_var_names[i]] =
          Partial_fixup_type::SHIFTED_CONSERVATIVE;
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
    @brief set repair method in partially filled cells for all variables
    @param fixup_type Can be Partial_fixup_type::CONSTANT,
                      Partial_fixup_type::LOCALLY_CONSERVATIVE,
                      Partial_fixup_type::SHIFTED_CONSERVATIVE
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
                       Partial_fixup_type::SHIFTED_CONSERVATIVE
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

  void set_num_tols(NumericTolerances_t num_tols) {
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
      std::cerr << "Type not supported \n";
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
      std::cerr << "Type not supported \n";
  }

#ifdef HAVE_TANGRAM
  /*!
    @brief set options for interface reconstructor driver  
    @param tols The vector of tolerances for each moment during reconstruction
    @param all_convex Should be set to false if the source mesh contains 
    non-convex cells.  
  */
  void set_reconstructor_options(std::vector<Tangram::IterativeMethodTolerances_t> &tols, 
                                 bool all_convex){
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
                 std::vector<std::string> const &source_meshvar_names,
                 std::vector<std::string> const &target_meshvar_names,
                 std::vector<std::string> const &source_matvar_names,
                 std::vector<std::string> const &target_matvar_names,
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
                 std::vector<std::string> const &source_meshvar_names,
                 std::vector<std::string> const &target_meshvar_names,
                 Wonton::Executor_type const *executor = nullptr);



  /*!

    @brief Detect and fix if we have a mismatch between source and
    target domain boundaries
    @tparam SourceMesh_Wrapper2  May be the mesh wrapper sent into MMDriver or the Flat_Mesh_Wrapper created for redistribution
    @tparam SourceState_Wrapper2 May be the state wrapper sent into MMDriver or the Flat_State_Wrapper created for redistribution
    @tparam entity_kind  Kind of entity that variables live on

    @param source_meshvar_names  names of remap variables on source mesh
    @param target_meshvar_names  names of remap variables on target mesh
    @param sources_and_weights   mesh-mesh intersection weights
    @param executor              pointer to Executor
    
  */
  
  template<class SourceMesh_Wrapper2, class SourceState_Wrapper2,
           Entity_kind onwhat>
  int fix_mismatch(SourceMesh_Wrapper2 const& source_mesh2,
                   SourceState_Wrapper2 const& source_state2,
                   Portage::vector<std::vector<Weights_t>> const& source_ents_and_weights,
                   std::vector<std::string> const& src_meshvar_names,
                   std::vector<std::string> const& trg_meshvar_names,
                   Wonton::Executor_type const *executor = nullptr) {
    
    // Will be null if it's a parallel executor
    auto serialexecutor = dynamic_cast<Wonton::SerialExecutor_type const *>(executor);

    bool distributed = false;
    int comm_rank = 0;
    int nprocs = 1;

#ifdef PORTAGE_ENABLE_MPI
    MPI_Comm mycomm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm, &comm_rank);
      MPI_Comm_size(mycomm, &nprocs);
      if (nprocs > 1)
        distributed = true;
    }
#endif

    MismatchFixer<D, onwhat, SourceMesh_Wrapper2,
                  SourceState_Wrapper2,
                  TargetMesh_Wrapper, TargetState_Wrapper>
        mismatch_fixer(source_mesh2, source_state2,
                       target_mesh_, target_state_,
                       source_ents_and_weights, executor);

    if (mismatch_fixer.has_mismatch()) {
      int nvars = src_meshvar_names.size();
      for (int i = 0; i < nvars; i++) {
        std::string const& src_var = src_meshvar_names[i];
        std::string const& trg_var = trg_meshvar_names[i];
      
        double lower_bound, upper_bound;
        try {  // see if we have caller specified bounds
          
          lower_bound = double_lower_bounds_.at(trg_var);
          upper_bound = double_upper_bounds_.at(trg_var);
          
        } catch (const std::out_of_range& oor) {
          // Since caller has not specified bounds for variable, attempt
          // to derive them from source state. This code should go into
          // Wonton into each state manager
          
          int nsrcents = source_mesh_.num_entities(onwhat,
                                                   Entity_type::PARALLEL_OWNED);
          
          double const *source_data;
          source_state_.mesh_get_data(onwhat, src_var, &source_data);
          lower_bound = *std::min_element(source_data, source_data + nsrcents);
          upper_bound = *std::max_element(source_data, source_data + nsrcents);
          
#ifdef PORTAGE_ENABLE_MPI
          if (mycomm != MPI_COMM_NULL) {
            double global_lower_bound=0.0, global_upper_bound=0.0;
            MPI_Allreduce(&lower_bound, &global_lower_bound, 1, MPI_DOUBLE,
                          MPI_MIN, mycomm);
            lower_bound = global_lower_bound;
            
            MPI_Allreduce(&upper_bound, &global_upper_bound, 1, MPI_DOUBLE,
                          MPI_MAX, mycomm);
            upper_bound = global_upper_bound;
          }
#endif

        double relbounddiff = fabs((upper_bound-lower_bound)/lower_bound);
        if (relbounddiff < consttol_) {
          // The field is constant over the source mesh/part. We HAVE to
          // relax the bounds to be able to conserve the integral quantity
          // AND maintain a constant.
          lower_bound -= 0.5*lower_bound;
          upper_bound += 0.5*upper_bound;
        }
        }
        
        double conservation_tol = DEFAULT_CONSERVATION_TOL;
        try {  // see if caller has specified a tolerance for conservation
          conservation_tol = conservation_tol_.at(trg_var);
        } catch ( const std::out_of_range& oor) {}
        
        mismatch_fixer.fix_mismatch(src_var, trg_var, lower_bound, upper_bound,
                                    conservation_tol, max_fixup_iter_,
                                    partial_fixup_types_[trg_var],
                                    empty_fixup_types_[trg_var]);
      }
    }
  }  // fix_mismatch



  /*!
    @brief Execute the remapping process
    @return status of remap (1 if successful, 0 if not)
  */
  int run(Wonton::Executor_type const *executor = nullptr,
          std::string *errmsg = nullptr) {
    std::string message;

    struct timeval begin_timeval, end_timeval, diff_timeval;

    bool distributed = false;
    int comm_rank = 0;
    int nprocs = 1;


    // Will be null if it's a parallel executor
    auto serialexecutor = dynamic_cast<Wonton::SerialExecutor_type const *>(executor);

#ifdef PORTAGE_ENABLE_MPI
    MPI_Comm mycomm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm, &comm_rank);
      MPI_Comm_size(mycomm, &nprocs);
      if (nprocs > 1)
        distributed = true;
    }
#endif
#ifdef ENABLE_DEBUG
    if (comm_rank == 0)
      std::cout << "in MMDriver::run()...\n";

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;
#endif

    int nvars = source_target_varname_map_.size();


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
#ifdef PORTAGE_ENABLE_MPI
    Flat_Mesh_Wrapper<> source_mesh_flat;
    Flat_State_Wrapper<Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);

    if (distributed) {

      // Create a new mesh wrapper that we can use for redistribution
      // of the source mesh as necessary (so that every target cell
      // sees any source cell that it overlaps with)

      // IN FACT, WE SHOULD DO THE BOUNDING BOX OVERLAP CHECK FIRST
      // AND ONLY IF WE DETERMINE THAT THE SOURCE MESH NEEDS TO BE
      // DISTRIBUTED WE SHOULD CREATE THE FLAT MESH WRAPPER AND INVOKE
      // REDISTRIBUTION; OTHERWISE, WE JUST INVOKE REMAP WITH THE
      // ORIGINAL WRAPPER

      gettimeofday(&begin_timeval, 0);

      source_mesh_flat.initialize(source_mesh_);

      // Note the flat state should be used for everything including the
      // centroids and volume fractions for interface reconstruction
      std::vector<std::string> source_remap_var_names;
      for (auto & stpair : source_target_varname_map_)
        source_remap_var_names.push_back(stpair.first);
      source_state_flat.initialize(source_state_, source_remap_var_names);

      MPI_Bounding_Boxes distributor(mpiexecutor);
      distributor.distribute(source_mesh_flat, source_state_flat,
                             target_mesh_, target_state_);

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
#ifdef ENABLE_DEBUG
      float tot_seconds_dist = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
      std::cout << "Redistribution Time Rank " << comm_rank << " (s): " <<
          tot_seconds_dist << std::endl;
#endif
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
          std::cerr << "Cannot handle multi-material fields on nodes\n";
      }
    }

    if (src_meshvar_names.size()) {
#ifdef PORTAGE_ENABLE_MPI
      if (distributed)
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
  NumericTolerances_t num_tols_;


#ifdef HAVE_TANGRAM
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
  std::vector<Tangram::IterativeMethodTolerances_t> reconstructor_tols_ = 
  {{1000, 1e-12, 1e-12}, {1000, 1e-12, 1e-12}};
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
  
  int comm_rank = 0;
  int nprocs = 1;

  // Will be null if it's a parallel executor
  auto serialexecutor = dynamic_cast<Wonton::SerialExecutor_type const *>(executor);

#ifdef PORTAGE_ENABLE_MPI
  MPI_Comm mycomm = MPI_COMM_NULL;
  auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
  if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
    mycomm = mpiexecutor->mpicomm;
    MPI_Comm_rank(mycomm, &comm_rank);
    MPI_Comm_size(mycomm, &nprocs);
  }
#endif


  
  float tot_seconds = 0.0, tot_seconds_srch = 0.0,
      tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;

  std::vector<std::string> source_remap_var_names;
  for (auto & stpair : source_target_varname_map_)
    source_remap_var_names.push_back(stpair.first);

  
  // Use default numerical tolerances in case they were not set earlier
  if (num_tols_.tolerances_set == false) {
    NumericTolerances_t default_num_tols;
    default_num_tols.use_default();
    set_num_tols(default_num_tols);
  }


  // Instantiate core driver

  Portage::CoreDriver<D, CELL,
                      SourceMesh_Wrapper2, SourceState_Wrapper2,
                      TargetMesh_Wrapper, TargetState_Wrapper,
                      InterfaceReconstructorType,
                      Matpoly_Splitter, Matpoly_Clipper>
      coredriver_cell(source_mesh2, source_state2, target_mesh_, target_state_);

  coredriver_cell.set_num_tols(num_tols_);
#ifdef HAVE_TANGRAM
  coredriver_cell.set_interface_reconstructor_options(reconstructor_tols_,
                                                      reconstructor_all_convex_);
#endif  
  
  // SEARCH

  auto candidates = coredriver_cell.template search<Portage::SearchKDTree>();

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  int nmats = source_state2.num_materials();

  //--------------------------------------------------------------------
  // REMAP MESH FIELDS FIRST (this requires just mesh-mesh intersection)
  //--------------------------------------------------------------------

  // INTERSECT MESHES

  gettimeofday(&begin_timeval, 0);

  auto source_ents_and_weights =
      coredriver_cell.template intersect_meshes<Intersect>(candidates);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_xsect += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;


  // INTERPOLATE (one variable at a time)
  gettimeofday(&begin_timeval, 0);
  int nvars = src_meshvar_names.size();
#ifdef ENABLE_DEBUG
  if (comm_rank == 0) {
    std::cout << "Number of mesh variables on cells to remap is " <<
        nvars << std::endl;
  }
#endif

  Portage::vector<Vector<D>> gradients;

  for (int i = 0; i < nvars; ++i) {
    std::string const& srcvar = src_meshvar_names[i];
    std::string const& trgvar = trg_meshvar_names[i];

    Limiter_type limiter = DEFAULT_LIMITER;
    auto const& it1 = limiters_.find(srcvar);
    if (it1 != limiters_.end()) limiter = it1->second;

    Boundary_Limiter_type bndlimiter = DEFAULT_BND_LIMITER;
    auto const& it2 = bnd_limiters_.find(srcvar);
    if (it2 != bnd_limiters_.end()) bndlimiter = it2->second;

    auto gradients =
        coredriver_cell.compute_source_gradient(srcvar, limiter, bndlimiter);
    
    coredriver_cell.template interpolate_mesh_var<double, Interpolate>
        (srcvar, trgvar, source_ents_and_weights, &gradients);
  }

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;


  // Fix mismatch in cell variables as requested
  
  fix_mismatch<SourceMesh_Wrapper2, SourceState_Wrapper2, CELL>
      (source_mesh2, source_state2, source_ents_and_weights,
       src_meshvar_names, trg_meshvar_names, executor);
  
  if (nmats > 1) {
    //--------------------------------------------------------------------
    // REMAP MULTIMATERIAL FIELDS NEXT, ONE MATERIAL AT A TIME
    //--------------------------------------------------------------------
    
    auto source_ents_and_weights_mat =
        coredriver_cell.template intersect_materials<Intersect>(candidates);
    
    int nmatvars = src_matvar_names.size();
    for (int i = 0; i < nmatvars; ++i) {
      std::string const& srcvar = src_matvar_names[i];
      std::string const& trgvar = trg_matvar_names[i];
      
      std::vector<Portage::vector<Vector<D>>> matgradients(nmats);
      
      Limiter_type limiter = DEFAULT_LIMITER;
      auto const& it1 = limiters_.find(srcvar);
      if (it1 != limiters_.end()) limiter = it1->second;
      
      Boundary_Limiter_type bndlimiter = DEFAULT_BND_LIMITER;
      auto const& it2 = bnd_limiters_.find(srcvar);
      if (it2 != bnd_limiters_.end()) bndlimiter = it2->second;
      
      for (int m = 0; m < nmats; m++)
        matgradients[m] =
            coredriver_cell.compute_source_gradient(src_matvar_names[i],
                                                    limiter, bndlimiter,
                                                    m);
      
      coredriver_cell.template interpolate_mat_var<double, Interpolate>
          (srcvar, trgvar, source_ents_and_weights_mat, &matgradients);
    }  // nmatvars
  }
  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;
#ifdef ENABLE_DEBUG
  std::cout << "Transform Time for Cell remap on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;
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
  
  int comm_rank = 0;
  int nprocs = 1;

  // Will be null if it's a parallel executor
  auto serialexecutor = dynamic_cast<Wonton::SerialExecutor_type const *>(executor);

#ifdef PORTAGE_ENABLE_MPI
  MPI_Comm mycomm = MPI_COMM_NULL;
  auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
  if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
    mycomm = mpiexecutor->mpicomm;
    MPI_Comm_rank(mycomm, &comm_rank);
    MPI_Comm_size(mycomm, &nprocs);
  }
#endif


  
  float tot_seconds = 0.0, tot_seconds_srch = 0.0,
      tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;

  std::vector<std::string> source_remap_var_names;
  for (auto & stpair : source_target_varname_map_)
    source_remap_var_names.push_back(stpair.first);

  
  // Use default numerical tolerances in case they were not set earlier
  if (num_tols_.tolerances_set == false) {
    NumericTolerances_t default_num_tols;
    default_num_tols.use_default();
    set_num_tols(default_num_tols);
  }


  // Instantiate core driver

  Portage::CoreDriver<D, NODE,
                      SourceMesh_Wrapper2, SourceState_Wrapper2,
                      TargetMesh_Wrapper, TargetState_Wrapper>
      coredriver_node(source_mesh2, source_state2, target_mesh_, target_state_);

  coredriver_node.set_num_tols(num_tols_);
#ifdef HAVE_TANGRAM
  coredriver_node.set_interface_reconstructor_options(reconstructor_tols_,
                                                      reconstructor_all_convex_);
#endif  
  
  // SEARCH

  auto candidates = coredriver_node.template search<Portage::SearchKDTree>();

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  int nmats = source_state2.num_materials();

  //--------------------------------------------------------------------
  // REMAP MESH FIELDS FIRST (this requires just mesh-mesh intersection)
  //--------------------------------------------------------------------

  // INTERSECT MESHES

  gettimeofday(&begin_timeval, 0);

  auto source_ents_and_weights =
      coredriver_node.template intersect_meshes<Intersect>(candidates);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_xsect += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;


  // INTERPOLATE (one variable at a time)
  gettimeofday(&begin_timeval, 0);
  int nvars = src_meshvar_names.size();
#ifdef ENABLE_DEBUG
  if (comm_rank == 0) {
    std::cout << "Number of mesh variables on nodes to remap is " <<
        nvars << std::endl;
  }
#endif

  Portage::vector<Vector<D>> gradients;

  for (int i = 0; i < nvars; ++i) {
    std::string const& srcvar = src_meshvar_names[i];
    std::string const& trgvar = trg_meshvar_names[i];

    Limiter_type limiter = DEFAULT_LIMITER;
    auto const& it1 = limiters_.find(srcvar);
    if (it1 != limiters_.end()) limiter = it1->second;

    Boundary_Limiter_type bndlimiter = DEFAULT_BND_LIMITER;
    auto const& it2 = bnd_limiters_.find(srcvar);
    if (it2 != bnd_limiters_.end()) bndlimiter = it2->second;

    auto gradients =
        coredriver_node.compute_source_gradient(srcvar, limiter, bndlimiter);
    
    coredriver_node.template interpolate_mesh_var<double, Interpolate>
        (srcvar, trgvar, source_ents_and_weights, &gradients);
  }

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;


  // Fix mismatch in cell variables as requested
    
  fix_mismatch<SourceMesh_Wrapper2, SourceState_Wrapper2, NODE>
      (source_mesh2, source_state2, source_ents_and_weights,
       src_meshvar_names, trg_meshvar_names, executor);
  

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp += diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;
#ifdef ENABLE_DEBUG
  std::cout << "Transform Time for Node remap on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;
#endif
  return 1;
}  // remap specialization for cells


}  // namespace Portage

#endif  // PORTAGE_DRIVER_MMDRIVER_H_
