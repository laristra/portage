



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
#include "portage/support/basis.h"
#include "portage/support/weight.h"
#include "portage/search/search_simple_points.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

/*!
  @file driver.h
  @brief Example driver for mapping between two meshes.

  This should serve as a good example for how to write your own driver routine
  and datastructures.
*/

namespace Portage {
namespace Meshfree {

/*!
  @class Driver "driver.h"
  @brief Driver provides the API to mapping from one mesh to another.
  @tparam SourceSwarm A lightweight wrapper to a specific input Swarm
  implementation that provides certain functionality.
  @tparam SourceState A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.
  @tparam TargetMesh A lightweight wrapper to a specific target Swarm
  implementation that provides certain functionality.
  @tparam TargetState A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.
*/
template <template <int, class, class> class Search,
          template <size_t, class, class> class Accumulate,
          template<size_t, class> class Estimate,
          int Dim,
          class SourceSwarm,
          class SourceState,
          class TargetSwarm = SourceSwarm,
          class TargetState = SourceState>
class SwarmDriver {
  
  // Something like this would be very helpful to users
  // static_assert(
  //   Dim == Interpolate::Dim,
  //   "The dimension of Driver and Interpolate do not match!"
  // );
  

 public:
  /*!
    @brief Constructor for running the Swarm driver where the kernel_type and
    support geometry are the same for every particle 

    @param[in] sourceSwarm A @c SourceSwarm to the source swarm.
    @param[in] sourceState A @c SourceState for the data that lives on
    the source swarm.
    @param[in] targetSwarm A @c TargetSwarm to the target swarm.
    @param[in,out] targetState A @c TargetState for the data that will
    be mapped to the target swarm.
    @param[in] smoothing_lengths Vector of smoothing lengths for each
    target particle - the three levels of vectors allow for a
    polygonal/polyhedral shaped support around each target particle and
    a vector of smoothing lengths for each facet of the
    polygonal/polyhedral support
    @param[in] kernel_type The type of weight function kernel (B4, SQUARE, etc)
    for any target particle
    @param[in] geom_type The geometry of the support (ELLIPTIC, TENSOR,
    FACETED) for any target particle
  */
  SwarmDriver(SourceSwarm const& sourceSwarm,
              SourceState const& sourceState,
              TargetSwarm const& targetSwarm,
              TargetState& targetState,
              std::vector<std::vector<std::vector<double>>> const& smoothing_lengths,
              Weight::Kernel const& kernel_type=Weight::B4,
              Weight::Geometry const& support_geom_type=Weight::ELLIPTIC,
	      WeightCenter center=Gather)
      : source_swarm_(sourceSwarm), source_state_(sourceState),
        target_swarm_(targetSwarm), target_state_(targetState),
    smoothing_lengths_(smoothing_lengths) {
           
    assert(Dim == sourceSwarm.space_dimension());
    assert(Dim == targetSwarm.space_dimension());

    weight_center_ = center;

    if (weight_center_ == Gather) {
      assert(smoothing_lengths_.size() == target_swarm_.num_particles());
      kernel_types_ = std::vector<Weight::Kernel>(target_swarm_.num_particles(),
						  kernel_type);
      geom_types_ = std::vector<Weight::Geometry>(target_swarm_.num_particles(),
						  support_geom_type);
    } else if (weight_center_ == Scatter) {
      assert(smoothing_lengths_.size() == source_swarm_.num_particles());
      kernel_types_ = std::vector<Weight::Kernel>(source_swarm_.num_particles(),
						  kernel_type);
      geom_types_ = std::vector<Weight::Geometry>(source_swarm_.num_particles(),
						  support_geom_type);
    }
  }

  /*!
    @brief Constructor for running the Swarm driver.
    @param[in] sourceSwarm A @c SourceSwarm to the source swarm.
    @param[in] sourceState A @c SourceState for the data that lives on
    the source swarm.
    @param[in] targetSwarm A @c TargetSwarm to the target swarm.
    @param[in,out] targetState A @c TargetState for the data that will
    be mapped to the target swarm.
    @param[in] smoothing_lengths Vector of smoothing lengths for each
    target particle - the three levels of vectors allow for a
    polygonal/polyhedral shaped support around each target particle and
    a vector of smoothing lengths for each facet of the
    polygonal/polyhedral support
    @param[in] kernel_types The type of weight function kernel (B4, SQUARE,
    etc) for each target particle
    @param[in] geom_types The geometry of the support (ELLIPTIC, TENSOR,
    FACETED) for each target particle
  */
  SwarmDriver(SourceSwarm const& sourceSwarm,
              SourceState const& sourceState,
              TargetSwarm const& targetSwarm,
              TargetState& targetState,
              std::vector<std::vector<std::vector<double>>> const& smoothing_lengths,
              std::vector<Weight::Kernel> const& kernel_types,
              std::vector<Weight::Geometry> const& geom_types,
	      WeightCenter center=Gather)
      : source_swarm_(sourceSwarm), source_state_(sourceState),
        target_swarm_(targetSwarm), target_state_(targetState),
    kernel_types_(kernel_types),
    geom_types_(geom_types),
    smoothing_lengths_(smoothing_lengths) {

    assert(sourceSwarm.space_dimension() == targetSwarm.space_dimension());
    
    weight_center_ = center;

    if (weight_center_ == Gather) {
      assert(smoothing_lengths_.size() == target_swarm_.num_particles());
      assert(kernel_types_.size() == target_swarm_.num_particles());
      assert(geom_types_.size() == target_swarm_.num_particles());
    } else if (weight_center_ == Scatter) {
      assert(smoothing_lengths_.size() == source_swarm_.num_particles());
      assert(kernel_types_.size() == source_swarm_.num_particles());
      assert(geom_types_.size() == source_swarm_.num_particles());
    }
  }

  /// Copy constructor (disabled)
  SwarmDriver(const SwarmDriver &) = delete;

  /// Assignment operator (disabled)
  SwarmDriver & operator = (const SwarmDriver &) = delete;

  /// Destructor
  ~SwarmDriver() {}

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source swarm to the target swarm.  This variable must
    exist in both swarms' state manager
  */
  void set_remap_var_names(std::vector<std::string> const &remap_var_names) {
    // remap variable names same in source and target swarm
    set_remap_var_names(remap_var_names, remap_var_names);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source swarm.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target swarm.
    @param[in] estimator_type Estimator to be used to remap variables
    (KernelDensity, LocalRegression)
    @param[in] basis   Order of the basis used to remap variables 
    (UNITARY, LINEAR, QUADRATIC)
  */
  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      EstimateType const& estimator_type = LocalRegression,
      Basis::Type const& basis_type = Basis::Unitary) {

    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    estimator_type_ = estimator_type;
    basis_type_ = basis_type;
  }


  /*!
    @brief Get the names of the variables to be remapped from the
    source swarm.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> source_remap_var_names() const {
    return source_remap_var_names_;
  }

  /*!
    @brief Get the names of the variables to be remapped to the
    target swarm.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> target_remap_var_names() const {
    return target_remap_var_names_;
  }

  /*!
    @brief Execute the remapping process
  */
  void run(bool distributed) {

    int comm_rank = 0;
#ifndef ENABLE_MPI
    if (distributed) {
      std::cout << "Request is for a parallel run but Portage is compiled for serial runs only\n";
      return;
    }
#endif

    // EVEN OTHERWISE WE DON'T HAVE THE ABILITY TO DO DISTRIBUTED SWARMS
    if (distributed) {
      std::cout << "CANNOT DO DISTRIBUTED REMAPPING WITH SWARMS\n";
      return;
    }

    int numSourcePts = source_swarm_.num_owned_particles();
    int numTargetPts = target_swarm_.num_owned_particles();
    std::cout << "Number of target cells in target swarm on rank "
              << comm_rank << ": "
              << numTargetPts << std::endl;

    int nvars = source_remap_var_names_.size();

    // Collect all variables and remap them

    std::vector<std::string> source_var_names;
    std::vector<std::string> target_var_names;
    for (int i = 0; i < nvars; ++i) {
      source_var_names.emplace_back(source_remap_var_names_[i]);
      target_var_names.emplace_back(target_remap_var_names_[i]);
    }

    if (source_var_names.size() > 0)
    {
      
      float tot_seconds = 0.0, tot_seconds_srch = 0.0,
          tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
      struct timeval begin_timeval, end_timeval, diff_timeval;
      
      // SEARCH
      Portage::vector<std::vector<unsigned int>> candidates(numTargetPts);
      
      // Get an instance of the desired search algorithm type which is expected
      // to be a functor with an operator() of the right form

      gettimeofday(&begin_timeval, 0);

      // search boxes around source points are of zero dimension but
      // those around target points are determined by particle
      // smoothing lengths

      auto sourceExtents =
          std::make_shared<std::vector<Point<Dim>>>(numSourcePts);
      auto targetExtents =
          std::make_shared<std::vector<Point<Dim>>>(numTargetPts);
      for (int i = 0; i < numTargetPts; i++)
        (*targetExtents)[i] = Point<Dim>(smoothing_lengths_[i][0]);

      const Search<Dim, SourceSwarm, TargetSwarm>
          searchfunctor(source_swarm_, target_swarm_,
                        sourceExtents, targetExtents,
			weight_center_);
      
      Portage::transform(target_swarm_.begin(PARTICLE, PARALLEL_OWNED),
                         target_swarm_.end(PARTICLE, PARALLEL_OWNED),
                         candidates.begin(), searchfunctor);

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
      


      // ACCUMULATE (build moment matrix, calculate shape functions)
      // EQUIVALENT TO INTERSECT IN MESH-MESH REMAP

      gettimeofday(&begin_timeval, 0);
    
      // Get an instance of the desired accumulate algorithm type which is 
      // expected to be a functor with an operator() of the right form

      const Accumulate<Dim, SourceSwarm, TargetSwarm>
          accumulateFunctor(source_swarm_, target_swarm_,
                            estimator_type_, weight_center_,
                            kernel_types_, geom_types_, smoothing_lengths_,
                            basis_type_);
    
      Portage::vector<std::vector<Weights_t>> source_pts_and_mults(numTargetPts);
      
      // For each particle in the target swarm get the shape functions
      // (multipliers for source particle values)
      
      Portage::transform(target_swarm_.begin(PARTICLE, PARALLEL_OWNED),
                         target_swarm_.end(PARTICLE, PARALLEL_OWNED),
                         candidates.begin(),
                         source_pts_and_mults.begin(),
                         accumulateFunctor);
      
      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds_xsect = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
      
      
      // ESTIMATE (one variable at a time)
      
      gettimeofday(&begin_timeval, 0);
      
      nvars = source_var_names.size();
      if (comm_rank == 0)
        std::cout << "number of variables to remap is " << nvars <<
            std::endl;
      
      // Get an instance of the desired interpolate algorithm type
      Estimate<Dim, SourceState> estimateFunctor(source_state_);
      
      for (int i = 0; i < nvars; ++i) {
        //amh: ?? add back accuracy output statement??
        if (comm_rank == 0) std::cout << "Remapping cell variable " <<
                                source_var_names[i] <<
                                " to variable " << target_var_names[i] <<
                                std::endl;

        estimateFunctor.set_variable(source_var_names[i]);
        
        // This populates targetField with the values returned by the
        // remapper operator
        
        // ***************** NOTE NOTE NOTE NOTE ********************
        // THE CURRENT SWARM_STATE DOES NOT HAVE AN OPERATOR TO RETURN
        // A RAW POINTER TO THE DATA. INSTEAD IT RETURNS THIS REQUIRES
        // A SPECIFIC TYPE OF THE SWARM_STATE CLASS. THIS IS A BAD
        // IDEA BECAUSE IT FORCES ALL OTHER SWARM_STATE CLASSES TO
        // HAVE THIS SAME DEFINITION. IT ALSO LEADS TO THE UGLY USAGE
        // WE SEE BELOW BECAUSE WE NEED A PORTAGE POINTER TO BE ABLE
        // TO USE THRUST

        typename SwarmState<Dim>::DblVecPtr target_field_shared_ptr;
        
        target_state_.get_field(target_var_names[i], target_field_shared_ptr);
        Portage::pointer<double> target_field(&((*target_field_shared_ptr)[0]));
        
        Portage::transform(target_swarm_.begin(PARTICLE, PARALLEL_OWNED),
                           target_swarm_.end(PARTICLE, PARALLEL_OWNED),
                           source_pts_and_mults.begin(),
                           target_field, estimateFunctor);
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
  SourceSwarm const& source_swarm_;
  TargetSwarm const& target_swarm_;
  SourceState const& source_state_;
  TargetState& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  EstimateType estimator_type_;
  WeightCenter weight_center_ = Gather;  // smoothing len. centered on trgt. pts
  std::vector<std::vector<std::vector<double>>> smoothing_lengths_;
  std::vector<Weight::Kernel> kernel_types_;
  std::vector<Weight::Geometry> geom_types_;
  Basis::Type basis_type_;
};  // class Driver_Swarm


}  // namespace Meshfree
}  // namespace Portage

#endif  // SRC_DRIVER_DRIVER_H_
