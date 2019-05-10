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
#include "portage/driver/driver_internal.h"


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

      for (auto onwhat : entity_kinds_) {
        switch (onwhat) {
          case Entity_kind::CELL:
            internal_driver_[onwhat] =
                make_internal_driver<D, Entity_kind::CELL,
                                     Flat_Mesh_Wrapper<>,
                                     Flat_State_Wrapper<Flat_Mesh_Wrapper<>>,
                                     TargetMesh, TargetState,
                                     Search, Intersect, Interpolate,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper>
                (onwhat, source_mesh_flat_, source_state_flat_,
                 target_mesh_, target_state_, field_types, executor);
            break;
          case Entity_kind::NODE:
            internal_driver_[onwhat] =
                make_internal_driver<D, Entity_kind::NODE,
                                     Flat_Mesh_Wrapper<>,
                                     Flat_State_Wrapper<Flat_Mesh_Wrapper<>>,
                                     TargetMesh, TargetState,
                                     Search, Intersect, Interpolate,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper>
                (onwhat, source_mesh_flat_, source_state_flat_,
                 target_mesh_, target_state_, field_types, executor);
            break;
          default:
            std::cerr << "Cannot handle remap on entity kind " << onwhat << "\n";
        }
      }
    }
    else
#endif
    {
      for (Wonton::Entity_kind onwhat : entity_kinds_) {
        switch (onwhat) {
          case Entity_kind::CELL:
            internal_driver_[onwhat] =
                make_internal_driver<D, Entity_kind::CELL,
                                     SourceMesh, SourceState,
                                     TargetMesh, TargetState,
                                     Search, Intersect, Interpolate,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper>
                (onwhat, source_mesh_, source_state_,
                 target_mesh_, target_state_, field_types, executor);
            break;
          case Entity_kind::NODE:
            internal_driver_[onwhat] =
                make_internal_driver<D, Entity_kind::NODE,
                                     SourceMesh, SourceState,
                                     TargetMesh, TargetState,
                                     Search, Intersect, Interpolate,
                                     InterfaceReconstructorType,
                                     Matpoly_Splitter, Matpoly_Clipper>
                (onwhat, source_mesh_, source_state_,
                 target_mesh_, target_state_, field_types, executor);
            break;
        }
      }
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


      internal_driver_[onwhat]->interpolate_var(srcvarname, trgvarname,
                                                limiters.at(srcvarname),
                                                partial_fixup_types.at(trgvarname),
                                                empty_fixup_types.at(trgvarname),
                                                lower_bound, upper_bound,
                                                conservation_tol,
                                                max_fixup_iter);
    }
  }  // interpolate_var

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
  
};  // PartDriver

}  // namespace Portage

#endif  // PORTAGE_DRIVER_PARTDRIVER_H_
