/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_DRIVER_PARTDRIVER_H_
#define PORTAGE_DRIVER_PARTDRIVER_H_

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
#include <map>

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"

#include "portage/intersect/dummy_interface_reconstructor.h"
#endif

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


#ifdef PORTAGE_ENABLE_MPI
#include "portage/distributed/mpi_bounding_boxes.h"
#endif

/*!
  @file mmdriver.h
  @brief Part by Part multi-material remapping driver

  Remap mesh and material variables from mesh to mesh (in serial or
  distributed settings) with the option of doing a part-by-part remap
  (remap between sets of source and target entities)
*/

namespace Portage {

using namespace Wonton;

/*!
  @class PartDriver "mmdriver.h"
  @brief PartDriver provides the API to mapping multi-material data from one mesh to another.

  @tparam Search  A search method that takes the dimension, source mesh class
  and target mesh class as template parameters

  @tparam Intersect  A polyhedron-polyhedron intersection class that takes
  the source and taget mesh classes as template parameters

  @tparam Interpolate An interpolation class that takes the source and
  target mesh classes, the source state class (that stores source
  field values), the kind of entity the interpolation is on and the
  dimension of the problem as template parameters

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
*/
template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class, class,
                    template <class, int, class, class> class,
                    class, class> class Intersect,
          template<int, Entity_kind, class, class, class,
                   template<class, int, class, class> class,
                   class, class> class Interpolate,
          int D,
          class SourceMesh,
          class SourceState,
          class TargetMesh = SourceMesh,
          class TargetState = SourceState,
          template <class, int, class, class> class InterfaceReconstructorType = DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void>
class PartDriver {

  class DriverInternalBase;
  
  // Something like this would be very helpful to users
  // static_assert(
  //   D == Interpolate::D,
  //   "The dimension of Driver and Interpolate do not match!"
  // );


 public:
  /*!
    @brief Constructor for the remap driver.
    @param[in] sourceMesh A @c  wrapper to the source mesh.
    @param[in] sourceState A @c wrapper for the data that lives on the
    source mesh
    @param[in] targetMesh A @c TargetMesh to the target mesh
    @param[in,out] targetState A @c TargetState for the data that will
    be mapped to the target mesh
  */
  PartDriver(SourceMesh const& source_mesh,
             SourceState const& source_state,
             TargetMesh const& target_mesh,
             TargetState& target_state,
             std::vector<Entity_kind> entity_kinds,
             std::vector<Field_type> field_types)
      : source_mesh_(source_mesh), source_state_(source_state),
        target_mesh_(target_mesh), target_state_(target_state),
        entity_kinds_(entity_kinds), field_types_(field_types),
        dim_(source_mesh.space_dimension()) {

    assert(source_mesh.space_dimension() == target_mesh.space_dimension());

  }

  /// Copy constructor (disabled)
  PartDriver(const PartDriver &) = delete;

  /// Assignment operator (disabled)
  PartDriver & operator = (const PartDriver &) = delete;

  /// Destructor
  ~PartDriver() {}

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
    @brief set limiter for all variables
    @param target_var_name Source mesh variable whose gradient is to be limited
    @param limiter  Limiter to use for second order reconstruction (NOLIMITER
    or BARTH_JESPERSEN)
  */
  void set_limiter(std::string const& source_var_name, Limiter_type limiter) {
    limiters_[source_var_name] = limiter;
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
    @brief Execute the remapping process
    @return status of remap (1 if successful, 0 if not)
  */
  int run(Wonton::Executor_type const *executor = nullptr,
          std::string *errmsg = nullptr) {

    // Do initialization which is an understated way of saying do all
    // the searches and intersects needed for the actual remapping
    
    init(entity_kinds_, field_types_, executor, nullptr);

    // interpolate the variables

    interpolate(source_target_varname_map_, limiters_,
                partial_fixup_types_, empty_fixup_types_,
                double_lower_bounds_, double_upper_bounds_,
                conservation_tol_, voldifftol_, consttol_, max_fixup_iter_);

    return 1;
  }  // run


  /*!
    @brief initialization routine for remap - completes search,
    intersect, interface reconstruction step
  */

  int init(std::vector<Entity_kind> const & entity_kinds,
           std::vector<Field_type> const & field_types,
           Wonton::Executor_type const *executor = nullptr,
           std::string *errmsg = nullptr) {
    std::string message;
    struct timeval begin_timeval, end_timeval, diff_timeval;

    
    // Figure out if we are running in distributed mode

    // Will be null if it's a parallel executor
    auto serialexecutor =
        dynamic_cast<Wonton::SerialExecutor_type const *>(executor);
    distributed_ = false;
    
#ifdef PORTAGE_ENABLE_MPI
    mycomm_ = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm_ = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm_, &comm_rank_);
      MPI_Comm_size(mycomm_, &nprocs_);
      if (nprocs_ > 1)
        distributed_ = true;
    }
#endif

    Portage::vector<double> test_vector;
    
    int numTargetCells = target_mesh_.num_owned_cells();

    // Default is serial run (if MPI is not enabled or the
    // communicator is not defined or the number of processors is 1)
#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {

      // Create a new mesh wrapper that we can use for redistribution
      // of the source mesh as necessary (so that every target cell
      // sees any source cell that it overlaps with)

      // IN FACT, WE SHOULD DO THE BOUNDING BOX OVERLAP CHECK FIRST
      // AND ONLY IF WE DETERMINE THAT THE SOURCE MESH NEEDS TO BE
      // DISTRIBUTED WE SHOULD CREATE THE FLAT MESH WRAPPER AND INVOKE
      // REDISTRIBUTION; OTHERWISE, WE JUST INVOKE REMAP WITH THE
      // ORIGINAL WRAPPER

      gettimeofday(&begin_timeval, 0);

      Flat_Mesh_Wrapper<> source_mesh_flat_;
      source_mesh_flat_.initialize(source_mesh_);

      std::vector<std::string> source_remap_var_names;
      for (auto & stpair : source_target_varname_map_)
        source_remap_var_names.push_back(stpair.first);

      Flat_State_Wrapper<Flat_Mesh_Wrapper<>> source_state_flat_(source_mesh_flat_);
      source_state_flat_.initialize(source_state_, source_remap_var_names);

      MPI_Bounding_Boxes distributor(mpiexecutor);
      distributor.distribute(source_mesh_flat_, source_state_flat_,
                             target_mesh_, target_state_);

      source_redistributed_ = true;

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      float tot_seconds_dist = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
      std::cout << "Redistribution Time Rank " << comm_rank_ << " (s): " <<
          tot_seconds_dist << std::endl;

      for (Entity_kind onwhat : entity_kinds_)
        internal_driver_[onwhat] =
            make_internal_driver(onwhat,
                                 source_mesh_flat_, source_state_flat_,
                                 target_mesh_, target_state_,
                                 field_types, executor);
    }
    else
#endif
    {
      for (Entity_kind onwhat : entity_kinds_)
        internal_driver_[onwhat] =
            make_internal_driver(onwhat,
                                 source_mesh_, source_state_,
                                 target_mesh_, target_state_,
                                 field_types, executor);
    }
  }  // PartDriver::init


  void interpolate(std::unordered_map<std::string, std::string> const &
                   source_target_varname_map,
                   std::unordered_map<std::string, Limiter_type> const &
                   limiters,
                   std::unordered_map<std::string, Partial_fixup_type> const &
                   partial_fixup_types,
                   std::unordered_map<std::string, Empty_fixup_type> const &
                   empty_fixup_types,
                   std::unordered_map<std::string, double> const &
                   double_lower_bounds,
                   std::unordered_map<std::string, double> const &
                   double_upper_bounds,
                   std::unordered_map<std::string, double> const &
                   conservation_tols,
                   double voldifftol, double consttol, int max_fixup_iter) {

    // INTERPOLATE MESH VARIABLE (one variable at a time)
    
    for (auto const& stpair : source_target_varname_map) {
      std::string const & srcvarname = stpair.first;
      std::string const & trgvarname = stpair.second;

      Entity_kind onwhat = source_state_.get_entity(srcvarname);
      Field_type field_type = source_state_.field_type(onwhat, srcvarname);

      if (field_type == Field_type::MESH_FIELD)
        std::cout << "Remapping mesh variable " << srcvarname <<
            " from source on " << onwhat <<
            " to " << trgvarname << " on target\n";
      else
        std::cout << "Remapping material variable " << srcvarname <<
            " from source on " << onwhat <<
            " to " << trgvarname << " on target\n";

        
      double lower_bound = std::numeric_limits<double>::max();
      double upper_bound = -lower_bound;
      try {  // see if we have caller specified bounds

        lower_bound = double_lower_bounds.at(trgvarname);
        upper_bound = double_upper_bounds.at(trgvarname);

      } catch (const std::out_of_range& oor) {
        // Since caller has not specified bounds for variable, attempt
        // to derive them from source state. This code should go into
        // Wonton into each state manager

        int nsrcents = source_mesh_.num_entities(onwhat,
                                                 Entity_type::PARALLEL_OWNED);

        double const *source_data;
        if (field_type == Field_type::MESH_FIELD) {
          source_state_.mesh_get_data(onwhat, srcvarname, &source_data);
          lower_bound = *std::min_element(source_data, source_data + nsrcents);
          upper_bound = *std::max_element(source_data, source_data + nsrcents);
        } else {
          // find the min and max over all materials
          int nmats = source_state_.num_materials();
          for (int m = 0; m < nmats; m++) {
            source_state_.mat_get_celldata(srcvarname, m, &source_data);
            lower_bound = std::min(lower_bound,
                                   *std::min_element(source_data,
                                                     source_data + nsrcents));
            upper_bound = std::max(upper_bound,
                                   *std::max_element(source_data,
                                                     source_data + nsrcents));
          }
        }
          
#ifdef PORTAGE_ENABLE_MPI
        if (mycomm_!= MPI_COMM_NULL) {
          double global_lower_bound = 0.0, global_upper_bound = 0.0;
          MPI_Allreduce(&lower_bound, &global_lower_bound, 1, MPI_DOUBLE,
                        MPI_MIN, mycomm_);
          lower_bound = global_lower_bound;

          MPI_Allreduce(&upper_bound, &global_upper_bound, 1, MPI_DOUBLE,
                        MPI_MAX, mycomm_);
          upper_bound = global_upper_bound;
        }
#endif

        double relbounddiff = fabs((upper_bound-lower_bound)/lower_bound);
        if (relbounddiff < consttol) {
          // The field is constant over the source mesh/part. We
          // HAVE to relax the bounds to be able to conserve the
          // integral quantity AND maintain a constant.
          lower_bound -= 0.5*lower_bound;
          upper_bound += 0.5*upper_bound;
        }
      }

      double conservation_tol = 100*std::numeric_limits<double>::epsilon();
      try {  // see if caller has specified a tolerance for conservation
        conservation_tol = conservation_tols.at(trgvarname);
      } catch ( const std::out_of_range& oor) {}


      interpolate_var(onwhat, srcvarname, trgvarname,
                      limiters.at(srcvarname),
                      partial_fixup_types.at(trgvarname),
                      empty_fixup_types.at(trgvarname),
                      lower_bound, upper_bound,
                      conservation_tol,
                      max_fixup_iter);
    }
  }  // interpolate

  template<typename T = double>
  void interpolate_var(Entity_kind onwhat,
                       std::string srcvarname, std::string trgvarname,
                       Limiter_type limiter,
                       Partial_fixup_type partial_fixup_type,
                       Empty_fixup_type empty_fixup_type,
                       T lower_bound, T upper_bound,
                       double conservation_tol,
                       int max_fixup_iter) {
    
    if (source_redistributed_) {
      switch (onwhat) {
        case Entity_kind::CELL: {
          auto internal_driver_onkind =
              dynamic_cast<DriverInternal<Entity_kind::CELL, Flat_Mesh_Wrapper<>, Flat_State_Wrapper<Flat_Mesh_Wrapper<>>> *>(internal_driver_[onwhat].get());

          internal_driver_onkind->interpolate_var<T>(srcvarname, trgvarname,
                                                     limiter,
                                                     partial_fixup_type,
                                                     empty_fixup_type,
                                                     lower_bound,
                                                     upper_bound,
                                                     conservation_tol,
                                                     max_fixup_iter);
          break;
        }
        case Entity_kind::NODE: {
          auto internal_driver_onkind =
              dynamic_cast<DriverInternal<Entity_kind::NODE, Flat_Mesh_Wrapper<>, Flat_State_Wrapper<Flat_Mesh_Wrapper<>>> *>(internal_driver_[onwhat].get());

          internal_driver_onkind->interpolate_var<T>(srcvarname, trgvarname,
                                                     limiter,
                                                     partial_fixup_type,
                                                     empty_fixup_type,
                                                     lower_bound,
                                                     upper_bound,
                                                     conservation_tol,
                                                     max_fixup_iter);
          break;
        }
      }
    } else {
      switch (onwhat) {
        case Entity_kind::CELL: {
          auto internal_driver_onkind =
              dynamic_cast<DriverInternal<Entity_kind::CELL, SourceMesh, SourceState> *>(internal_driver_[onwhat].get());

          internal_driver_onkind->interpolate_var<T>(srcvarname, trgvarname,
                                                     limiter,
                                                     partial_fixup_type,
                                                     empty_fixup_type,
                                                     lower_bound,
                                                     upper_bound,
                                                     conservation_tol,
                                                     max_fixup_iter);
          break;
        }
        case Entity_kind::NODE: {
          auto internal_driver_onkind =
              dynamic_cast<DriverInternal<Entity_kind::NODE, SourceMesh, SourceState> *>(internal_driver_[onwhat].get());
          
          internal_driver_onkind->interpolate_var<T>(srcvarname, trgvarname,
                                                     limiter,
                                                     partial_fixup_type,
                                                     empty_fixup_type,
                                                     lower_bound,
                                                     upper_bound,
                                                     conservation_tol,
                                                     max_fixup_iter);
          break;
        }
      }
    }
  }

 private:

  // Inputs specified by calling app
  Wonton::Executor_type const *executor_;
  SourceMesh const& source_mesh_;
  TargetMesh const& target_mesh_;
  SourceState const& source_state_;
  TargetState& target_state_;
  std::unordered_map<std::string, std::string> source_target_varname_map_;
  std::unordered_map<std::string, Limiter_type> limiters_;
  std::unordered_map<std::string, Partial_fixup_type> partial_fixup_types_;
  std::unordered_map<std::string, Empty_fixup_type> empty_fixup_types_;
  std::unordered_map<std::string, double> double_lower_bounds_;
  std::unordered_map<std::string, double> double_upper_bounds_;
  std::unordered_map<std::string, double> conservation_tol_;
  unsigned int dim_;
  double voldifftol_ = 100*std::numeric_limits<double>::epsilon();
  double consttol_ =  100*std::numeric_limits<double>::epsilon();
  int max_fixup_iter_ = 5;


  // Internal variables
  bool distributed_ = false;  // default is serial
  bool source_redistributed_ = false;  // Did we redistribute the source mesh?
  int comm_rank_ = 0;          
  int nprocs_ = 1;
  
#ifdef PORTAGE_ENABLE_MPI
  MPI_Comm mycomm_ = MPI_COMM_NULL;
#endif

  std::vector<Entity_kind> entity_kinds_;
  std::vector<Field_type> field_types_;

  // Internal driver with the original or redistributed source mesh and state

  std::map<Entity_kind, std::unique_ptr<DriverInternalBase>> internal_driver_;


  /*!
    @file driver_internal.h
    @brief Internal driver for remapping

    Remap mesh and material variables from mesh to mesh with the option
    of doing a part-by-part remap (remap between sets of source and
    target entities)
  */


  /*!
    @class DriverInternalBase "driver_internal.h"

    @brief DriverInternalBase - Base class for internal driver that is
    agnostic to the Entity_kind (see derived class DriverInternal for 
    documentation of template parameters and methods)
  */

  class DriverInternalBase {
   public:
    DriverInternalBase() {}
    virtual ~DriverInternalBase() {}   // Necessary as the instance of
                                       // the derived class will get
                                       // destroyed through a pointer
                                       // to the base class
  };


  /*!
    @class DriverInternal "driver_internal.h"

    @brief DriverInternal - Internal driver that remaps fields on a particular Entity_kind (ONWHAT) like CELL or NODE

    NOTE: THIS CLASS ASSUMES THAT ALL SOURCE CELLS OVERLAPPING ANY TARGET
    CELL ARE AVAILABLE ON THIS PROCESSOR. IT DOES NOT HAVE TO FETCH THE
    DATA FROM ANYWHERE

    @tparam ONWHAT     On what kind of entity are we doing the remap

    @tparam SourceMesh A lightweight wrapper to a specific input mesh
    implementation that provides certain functionality.

    @tparam SourceState A lightweight wrapper to a specific input state
    manager implementation that provides certain functionality.

  */
  template <Entity_kind ONWHAT,
            class SourceMesh2,
            class SourceState2>
  class DriverInternal : public DriverInternalBase {
   public:
    /*!
      @brief Constructor for the INTERNAL remap driver.
      @param[in] sourceMesh A @c  wrapper to the source mesh (may be native or redistributed source).
      @param[in] sourceState A @c wrapper for the data that lives on the
      source mesh
      @param[in] targetMesh A @c TargetMesh to the target mesh
      @param[in,out] targetState A @c TargetState for the data that will
      be mapped to the target mesh
    */
    DriverInternal(SourceMesh2 const& source_mesh,
                   SourceState2 const& source_state,
                   TargetMesh const& target_mesh,
                   TargetState& target_state,
                   std::vector<Field_type> const & field_types,
                   Wonton::Executor_type const *executor = nullptr)
        : source_mesh2_(source_mesh), source_state2_(source_state),
          target_mesh_(target_mesh), target_state_(target_state),
          field_types_(field_types), executor_(executor)
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

      initialize();   // Do all the setup, search and intersect
    }

    /// Copy constructor (disabled)
    DriverInternal(const DriverInternal &) = delete;

    /// Assignment operator (disabled)
    DriverInternal & operator = (const DriverInternal &) = delete;

    /// Destructor
    ~DriverInternal() {}


    /*! DriverInternal::interpolate_var

      @brief interpolate a variable
      @param[in] srcvarname  Variable name on the source mesh
      @param[in] trgvarname  Variable name on the target mesh
      @param[in] limiter     Limiter to use for variable
      @param[in] partial...  How to fixup partly filled target cells (for this var)
      @param[in] emtpy...    How to fixup empty target cells with this var
      @param[in] lower_bound Lower bound of variable value when doing fixup
      @param[in] upper_bound Upper bound of variable value when doing fixup
      @param[in] cons..tol   Tolerance for conservation when doing fixup
    */

    template <typename T = double>
    void interpolate_var(std::string srcvarname, std::string trgvarname,
                         Limiter_type limiter,
                         Partial_fixup_type partial_fixup_type,
                         Empty_fixup_type empty_fixup_type,
                         T lower_bound, T upper_bound,
                         double conservation_tol,
                         int max_fixup_iter) {

      // Until we fix flat_state_mm_wrapper to make get_entity a const method,
      // comment this out
      // if (source_state2_.get_entity(srcvarname) != ONWHAT) {
      //   std::cerr << "Variable " << srcvarname <<
      //       " not defined on entity kind " << ONWHAT << ". Skipping!\n";
      //   return;
      // }

      Field_type field_type = source_state2_.field_type(ONWHAT, srcvarname);

      if (field_type == Field_type::MULTIMATERIAL_FIELD)
        interpolate_mat_var<T>(srcvarname, trgvarname, limiter,
                               partial_fixup_type, empty_fixup_type,
                               lower_bound, upper_bound, conservation_tol,
                               max_fixup_iter);
      else
        interpolate_mesh_var<T>(srcvarname, trgvarname, limiter,
                                partial_fixup_type, empty_fixup_type,
                                lower_bound, upper_bound, conservation_tol,
                                max_fixup_iter);

    }  // DriverInternal::interpolate_var


   private:
    SourceMesh2 const & source_mesh2_;
    TargetMesh const & target_mesh_;
    SourceState2 const & source_state2_;
    TargetState & target_state_;

    std::vector<Field_type> const & field_types_;

    int comm_rank_ = 0;
    int nprocs_ = 1;

    Wonton::Executor_type const *executor_;

#ifdef PORTAGE_ENABLE_MPI
    MPI_Comm mycomm_ = MPI_COMM_NULL;
#endif


    // intersection candidates for target entities
    //
    // for all target entities
    //    of a particular kind
    //          ||
    //          ||     candidate list for
    //          ||     each target entity
    //          ||           ||
    //          \/           \/
    Portage::vector<std::vector<int>> candidates_;

    // Weights of intersection b/w target entities and source entities
    // Each intersection is between the control volume (cell, dual cell)
    // of a target and source entity.
    //  for all target entities
    //     of a particular kind
    //           ||
    //           ||     intersection moments list
    //           ||     for each target entity
    //           ||           ||
    //           \/           \/
    Portage::vector<std::vector<Weights_t>> source_ents_and_weights_;

    // Weights of intersection b/w target CELLS and source material polygons
    // Each intersection is between a target cell and material polygon in
    // a source cell for a particular material

    //  for each material
    //   ||
    //   ||      for all target entities
    //   ||         of a particular kind
    //   ||               ||
    //   ||               ||     intersection moments for
    //   ||               ||     each target entity
    //   ||               ||           ||
    //   \/               \/           \/
    std::vector<Portage::vector<std::vector<Weights_t>>>
    source_ents_and_weights_by_mat_;

    // cell indices of target mesh organized by material - gets
    // populated in init()

    std::vector<std::vector<int>> matcellstgt_;


    // Pointer to the interface reconstructor object (required by the
    // interface to be shared)
    std::shared_ptr<Tangram::Driver<InterfaceReconstructorType, D,
                                    SourceMesh2,
                                    Matpoly_Splitter, Matpoly_Clipper>
                    >
    interface_reconstructor_;
  

    // Pointers to interpolators (will be used as needed per variable)

    std::unique_ptr<Interpolate<D, ONWHAT,
                                SourceMesh2, TargetMesh, SourceState2,
                                InterfaceReconstructorType,
                                Matpoly_Splitter, Matpoly_Clipper>
                    >
    interpolator_;
  
    // Pointers to mismatch fixers
    std::unique_ptr<MismatchFixer<D, ONWHAT,
                                  SourceMesh2, SourceState2,
                                  TargetMesh, TargetState>
                    >
    mismatch_fixer_;

    // Flag to indicate if meshes have mismatch
    bool has_mismatch_ = false;


    // Internal driver methods

    void initialize()  {
      // Compute and store search candidates
      compute_search_candidates();

      // Compute source weights (from intersection moments)
      compute_source_weights_by_intersection();

      // Instantiate mismatch fixer for later use
      mismatch_fixer_ = std::make_unique<MismatchFixer<D, ONWHAT,
                                                       SourceMesh2, SourceState2,
                                                       TargetMesh,  TargetState>
                                         >
          (source_mesh2_, source_state2_, target_mesh_, target_state_,
           source_ents_and_weights_, executor_);

      has_mismatch_ |= mismatch_fixer_->has_mismatch();

      // Instantiate interpolator for later use    
      interpolator_ = std::make_unique<Interpolate<D, ONWHAT,
                                                   SourceMesh2, TargetMesh,
                                                   SourceState2,
                                                   InterfaceReconstructorType,
                                                   Matpoly_Splitter,
                                                   Matpoly_Clipper>
                                       >
          (source_mesh2_, target_mesh_, source_state2_, interface_reconstructor_);

    }



#ifdef HAVE_TANGRAM

    // Convert volume fraction and centroid data from compact
    // material-centric to compact cell-centric (ccc) form as needed
    // by Tangram
    void ccc_vfcen_data(std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Tangram::Point<D>>& cell_mat_centroids) {

      int nsourcecells = source_mesh2_.num_entities(Entity_kind::CELL,
                                                    Entity_type::ALL);

      int nmats = source_state2_.num_materials();
      cell_num_mats.assign(nsourcecells, 0);

      // First build full arrays (as if every cell had every material)

      std::vector<int> cell_mat_ids_full(nsourcecells*nmats, -1);
      std::vector<double> cell_mat_volfracs_full(nsourcecells*nmats, 0.0);
      std::vector<Tangram::Point<D>> cell_mat_centroids_full(nsourcecells*nmats);

      int nvals = 0;
      for (int m = 0; m < nmats; m++) {
        std::vector<int> cellids;
        source_state2_.mat_get_cells(m, &cellids);
        for (int ic = 0; ic < cellids.size(); ic++) {
          int c = cellids[ic];
          int nmatc = cell_num_mats[c];
          cell_mat_ids_full[c*nmats+nmatc] = m;
          cell_num_mats[c]++;
        }
        nvals += cellids.size();

        double const * matfracptr;
        source_state2_.mat_get_celldata("mat_volfracs", m, &matfracptr);
        for (int ic = 0; ic < cellids.size(); ic++)
          cell_mat_volfracs_full[cellids[ic]*nmats+m] = matfracptr[ic];

        Portage::Point<D> const *matcenvec;
        source_state2_.mat_get_celldata("mat_centroids", m, &matcenvec);
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

    /*!
      DriverInternal::compute_search_candidates_on_kind<Entity_kind>

      Find candidates entities of a particular kind that might
      intersect each target entity of the same kind
    */

    int compute_search_candidates() {
      // Get an instance of the desired search algorithm type
      const Search<D, ONWHAT, SourceMesh2, TargetMesh>
          search_functor(source_mesh2_, target_mesh_);

      int ntarget_ents = target_mesh_.num_entities(ONWHAT,
                                                   Entity_type::PARALLEL_OWNED);

      // initialize search candidate vector of particular entity_kind
      candidates_.resize(ntarget_ents);

      Portage::transform(target_mesh_.begin(ONWHAT, Entity_type::PARALLEL_OWNED),
                         target_mesh_.end(ONWHAT, Entity_type::PARALLEL_OWNED),
                         candidates_.begin(), search_functor);
    }


    /*! DriverInternal::compute_intersections

      Compute the moments of intersection between each target entity and
      source entity for each entity_kind on which we have to remap
    */


    int compute_source_weights_by_intersection() {

      // Compute source-target mesh intersections
      intersect_meshes();

      // Compute intersections of target mesh with source material
      // polygons No-op if it is a single material problem
      if (ONWHAT == Entity_kind::CELL)
        intersect_materials();

    }  // compute_source_weights_intersection


    /*! DriverInternal::intersect_meshes

      @brief Intersect source and target mesh entities of kind
      'ONWHAT' and return the intersecting entities and moments of
      intersection for each entity
    */

    int intersect_meshes() {

      Intersect<ONWHAT, SourceMesh2, SourceState2,
                TargetMesh, DummyInterfaceReconstructor,
                void, void>
          intersector(source_mesh2_, source_state2_, target_mesh_);

      int ntarget_ents =
          target_mesh_.num_entities(ONWHAT, Entity_type::PARALLEL_OWNED);

      source_ents_and_weights_.resize(ntarget_ents);

      Portage::transform(target_mesh_.begin(ONWHAT, Entity_type::PARALLEL_OWNED),
                         target_mesh_.end(ONWHAT, Entity_type::PARALLEL_OWNED),
                         candidates_.begin(),
                         source_ents_and_weights_.begin(),
                         intersector);
    }



    /*! DriverInternal::intersect_materials

      @brief Intersect target mesh cells with source material polygons
      and return the intersecting entities and moments of intersection
      for each entity

    */

    int intersect_materials() {

      int nmats = source_state2_.num_materials();
      if (nmats == 0) return 1;

#ifdef HAVE_TANGRAM

      // Make sure we have a valid interface reconstruction method instantiated

      assert(typeid(InterfaceReconstructorType<SourceMesh2, D,
                    Matpoly_Splitter, Matpoly_Clipper >) !=
             typeid(DummyInterfaceReconstructor<SourceMesh2, D,
                    Matpoly_Splitter, Matpoly_Clipper>));

      std::vector<Tangram::IterativeMethodTolerances_t>
          tols(2, {1000, 1e-12, 1e-12});

      interface_reconstructor_ =
          std::make_unique<Tangram::Driver<InterfaceReconstructorType, D,
                                           SourceMesh2,
                                           Matpoly_Splitter,
                                           Matpoly_Clipper>
                           >(source_mesh2_, tols, true);
    
      int nsourcecells = source_mesh2_.num_entities(Entity_kind::CELL,
                                                    Entity_type::ALL);
      int ntargetcells = target_mesh_.num_entities(Entity_kind::CELL,
                                                   Entity_type::PARALLEL_OWNED);


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

      Intersect<Entity_kind::CELL, SourceMesh2, SourceState2, TargetMesh,
                InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper>
          intersector(source_mesh2_, source_state2_, target_mesh_,
                      interface_reconstructor_);

      // Assume (with no harm for sizing purposes) that all materials
      // in source made it into target (MAYBE WE SHOULD USE A MAP)

      matcellstgt_.resize(nmats);
      source_ents_and_weights_by_mat_.resize(nmats);

      for (int m = 0; m < nmats; m++) {

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
        Portage::transform(target_mesh_.begin(Entity_kind::CELL,
                                              Entity_type::PARALLEL_OWNED),
                           target_mesh_.end(Entity_kind::CELL,
                                            Entity_type::PARALLEL_OWNED),
                           candidates_.begin(),
                           this_mat_sources_and_wts.begin(),
                           intersector);

        // LOOK AT INTERSECTION WEIGHTS TO DETERMINE WHICH TARGET CELLS
        // WILL GET NEW MATERIALS

        int ntargetcells = target_mesh_.num_entities(Entity_kind::CELL,
                                                     Entity_type::ALL);

        for (int c = 0; c < ntargetcells; c++) {
          std::vector<Weights_t> const& cell_mat_sources_and_weights =
              this_mat_sources_and_wts[c];
          int nwts = cell_mat_sources_and_weights.size();
          for (int s = 0; s < nwts; s++) {
            std::vector<double> const& wts = cell_mat_sources_and_weights[s].weights;
            if (wts[0] > 0.0) {
              double vol = target_mesh_.cell_volume(c);
              if (wts[0]/vol > 1.0e-10) {  // Check that the volume of material
                //                         // we are adding is not miniscule
                matcellstgt_[m].push_back(c);
                break;
              }
            }
          }
        }

        // If any processor is adding this material to the target state,
        // add it on all the processors

        int nmatcells = matcellstgt_[m].size();
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
            if (target_state_.material_name(i) == source_state2_.material_name(m)) {
              found = true;
              m2 = i;
              break;
            }
          if (found) {  // material already present - just update its cell list
            target_state_.mat_add_cells(m2, matcellstgt_[m]);
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

            target_state_.add_material(source_state2_.material_name(m),
                                       matcellstgt_[m]);
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

        source_ents_and_weights_by_mat_[m].resize(nmatcells);

        for (int ic = 0; ic < nmatcells; ic++) {
          int c = matcellstgt_[m][ic];
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

          source_ents_and_weights_by_mat_[m][ic] = cell_mat_sources_and_weights;
        }

        target_state_.mat_add_celldata("mat_volfracs", m, &(mat_volfracs[0]));
        target_state_.mat_add_celldata("mat_centroids", m, &(mat_centroids[0]));

      }  // for each material m

#endif
    }


    /*! DriverInternal::interpolate_mesh_var

      @brief interpolate a mmesh variable
      @param[in] srcvarname  Mesh variable name on the source mesh
      @param[in] trgvarname  Mesh variable name on the target mesh
      @param[in] limiter     Limiter to use for variable
      @param[in] partial...  How to fixup partly filled target cells (for this var)
      @param[in] emtpy...    How to fixup empty target cells with this var
      @param[in] lower_bound Lower bound of variable value when doing fixup
      @param[in] upper_bound Upper bound of variable value when doing fixup
      @param[in] cons..tol   Tolerance for conservation when doing fixup
    */

    template <typename T = double>
    void interpolate_mesh_var(std::string srcvarname, std::string trgvarname,
                              Limiter_type limiter,
                              Partial_fixup_type partial_fixup_type,
                              Empty_fixup_type empty_fixup_type,
                              T lower_bound, T upper_bound,
                              double conservation_tol,
                              int max_fixup_iter) {

      // Until flat_state_mm_wrapper fixes get_entity to be const,
      // comment this out
      // if (source_state2_.get_entity(srcvarname) != ONWHAT) {
      //   std::cerr << "Variable " << srcvarname << " not defined on Entity_kind "
      //             << ONWHAT << ". Skipping!\n";
      //   return;
      // }


      // Get a handle to a memory location where the target state
      // would like us to write this material variable into. If it is
      // NULL, we allocate it ourself

      T *target_field_raw;
      target_state_.mesh_get_data(ONWHAT, trgvarname, &target_field_raw);
      assert(target_field_raw != nullptr);
      Portage::pointer<T> target_field(target_field_raw);


      interpolator_->set_interpolation_variable(srcvarname, limiter);

      Portage::transform(target_mesh_.begin(ONWHAT, Entity_type::PARALLEL_OWNED),
                         target_mesh_.end(ONWHAT, Entity_type::PARALLEL_OWNED),
                         source_ents_and_weights_.begin(),
                         target_field, *interpolator_);

      if (has_mismatch_)
        mismatch_fixer_->fix_mismatch(srcvarname, trgvarname,
                                      lower_bound, upper_bound,
                                      conservation_tol,
                                      max_fixup_iter,
                                      partial_fixup_type,
                                      empty_fixup_type);
    }


    /*! DriverInternal::interpolate_mat_var

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

    template<typename T = double>
    void interpolate_mat_var(std::string srcvarname, std::string trgvarname,
                             Limiter_type limiter,
                             Partial_fixup_type partial_fixup_type,
                             Empty_fixup_type empty_fixup_type,
                             T lower_bound, T upper_bound,
                             double conservation_tol,
                             int max_fixup_iter) {

      int nmats = source_state2_.num_materials();


      for (int m = 0; m < nmats; m++) {

        interpolator_->set_material(m);    // We have to do this so we know
        //                                 // which material values we have
        //                                 // to grab from the source state

        // FEATURE ;-)  Have to set interpolation variable AFTER setting 
        // the material for multimaterial variables

        interpolator_->set_interpolation_variable(srcvarname, limiter);

        // if the material has no cells on this partition, then don't bother
        // interpolating MM variables
        if (target_state_.mat_get_num_cells(m) == 0) continue;

        // Get a handle to a memory location where the target state
        // would like us to write this material variable into. If it is
        // NULL, we allocate it ourself

        T *target_field_raw;
        target_state_.mat_get_celldata(trgvarname, m, &target_field_raw);
        assert (target_field_raw != nullptr);

        Portage::pointer<T> target_field(target_field_raw);

        Portage::transform(matcellstgt_[m].begin(), matcellstgt_[m].end(),
                           source_ents_and_weights_by_mat_[m].begin(),
                           target_field, *interpolator_);

        // If the state wrapper knows that the target data is already
        // laid out in this way and it gave us a pointer to the array
        // where the values reside, it has to do nothing in this
        // call. If the storage format is different, however, it may
        // have to copy the values into their proper locations

        target_state_.mat_add_celldata(trgvarname, m, target_field_raw);
      }  // over all mats

    }  // DriverInternal::interpolate_mat_var
  };  // DriverInternal

  
  template <class SourceMesh2, class SourceState2>
  std::unique_ptr<DriverInternalBase>    // return type
  make_internal_driver(Entity_kind onwhat,
                       SourceMesh2 const & source_mesh2,
                       SourceState2 const & source_state2,
                       TargetMesh const & target_mesh,
                       TargetState & target_state,
                       std::vector<Field_type> const & field_types,
                       Wonton::Executor_type const *executor) {

    switch (onwhat) {
      case Entity_kind::CELL:
        return
            std::make_unique<DriverInternal<Entity_kind::CELL,
                                            SourceMesh2, SourceState2>>
            (source_mesh2, source_state2, target_mesh, target_state,
             field_types, executor);
      case Entity_kind::NODE:
        return
            std::make_unique<DriverInternal<Entity_kind::NODE,
                                            SourceMesh2, SourceState2>>
            (source_mesh2, source_state2, target_mesh, target_state,
             field_types, executor);
      default:
        std::cerr << "Remapping on entity kind " << onwhat << " not implemented\n";
    }
  }

};  // PartDriver

}  // namespace Portage

#endif  // PORTAGE_DRIVER_PARTDRIVER_H_
