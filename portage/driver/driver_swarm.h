/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_DRIVER_SWARM_H_
#define SRC_DRIVER_SWARM_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>

#include "portage/support/portage.h"

#include "portage/support/basis.h"
#include "portage/support/weight.h"
#include "portage/support/operator.h"
#include "portage/search/search_simple_points.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

#ifdef ENABLE_MPI
#include "portage/distributed/mpi_particle_distribute.h"
#endif 

/*!
  @file driver_swarm.h
  @brief Example driver for mapping between two swarms.

  This should serve as a good example for how to write your own driver routine
  and datastructures.
*/

namespace Portage {
namespace Meshfree {

/*!
  @class Driver "driver_swarm.h"
  @brief Driver provides the API to remap variablesfrom one swarm to another.
  @tparam SourceSwarm A lightweight wrapper to a specific input Swarm
  implementation that provides certain functionality.
  @tparam SourceState A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.
  @tparam TargetSwarm A lightweight wrapper to a specific target Swarm
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

    @c smoothing_lengths must have the size of @c SourceSwarm if center is @c Scatter and 
    of @c TargetSwarm if center is @c Gather.
  */
  SwarmDriver(SourceSwarm& sourceSwarm,
              SourceState& sourceState,
              TargetSwarm const& targetSwarm,
              TargetState& targetState,
              vector<std::vector<std::vector<double>>> const& smoothing_lengths,
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
      assert(smoothing_lengths_.size() == target_swarm_.num_particles(PARALLEL_OWNED));
      kernel_types_ = vector<Weight::Kernel>(target_swarm_.num_particles(PARALLEL_OWNED),
                                                  kernel_type);
      geom_types_ = vector<Weight::Geometry>(target_swarm_.num_particles(PARALLEL_OWNED),
                                                  support_geom_type);
    } else if (weight_center_ == Scatter) {
      assert(smoothing_lengths_.size() == source_swarm_.num_particles(PARALLEL_OWNED));
      kernel_types_ = vector<Weight::Kernel>(source_swarm_.num_particles(PARALLEL_OWNED),
                                                  kernel_type);
      geom_types_ = vector<Weight::Geometry>(source_swarm_.num_particles(PARALLEL_OWNED),
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
  SwarmDriver(SourceSwarm& sourceSwarm,
              SourceState& sourceState,
              TargetSwarm const& targetSwarm,
              TargetState& targetState,
              vector<std::vector<std::vector<double>>> const& smoothing_lengths,
              vector<Weight::Kernel> const& kernel_types,
              vector<Weight::Geometry> const& geom_types,
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
      std::vector<std::string> const &source_remap_var_names,
      std::vector<std::string> const &target_remap_var_names,
      EstimateType const estimator_type = LocalRegression,
      Basis::Type const basis_type = Basis::Unitary,
      Operator::Type const operator_spec = Operator::LastOperator,
      Portage::vector<Operator::Domain> const &operator_domains = vector<Operator::Domain>(0),
      Portage::vector<std::vector<Point<Dim>>> const &operator_data=
        vector<std::vector<Point<Dim>>>(0,std::vector<Point<Dim>>(0))) 
  {
    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    estimator_type_ = estimator_type;
    basis_type_ = basis_type;
    operator_spec_ = operator_spec;
    operator_domains_ = operator_domains;
    operator_data_ = operator_data;
    if (operator_spec_ != Operator::LastOperator) { 
      assert(operator_domains_.size() == target_swarm_.num_owned_particles());
      assert(operator_data_.size() == target_swarm_.num_owned_particles());
    }
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

  void remap(std::vector<std::string> const &src_varnames,
             std::vector<std::string> const &trg_varnames,
             bool distributed,
             bool report_time);

  /*!
    @brief Execute the remapping process
  */
  void run(bool distributed, bool report_time=true) {

#ifndef ENABLE_MPI
    if (distributed) {
      std::cout << "Request is for a parallel run but Portage is compiled " <<
          "for serial runs only\n";
      return;
    }
#endif

    int comm_rank = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

    if (comm_rank == 0)
      std::cout << "in SwarmDriver::run()...\n";

    int numTargetPts = target_swarm_.num_owned_particles();
    std::cout << "Number of target particles in target swarm on rank "
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

    remap(source_var_names, target_var_names, distributed, report_time);
    
  }  // run

 private:
  SourceSwarm& source_swarm_;
  TargetSwarm const& target_swarm_;
  SourceState& source_state_;
  TargetState& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  WeightCenter weight_center_ = Gather;
  vector<std::vector<std::vector<double>>> smoothing_lengths_;
  vector<Weight::Kernel> kernel_types_;
  vector<Weight::Geometry> geom_types_;
  EstimateType estimator_type_;
  Basis::Type basis_type_;
  Operator::Type operator_spec_;
  Portage::vector<Operator::Domain> operator_domains_;
  Portage::vector<std::vector<Point<Dim>>> operator_data_;
};  // class Driver_Swarm

template <template <int, class, class> class Search,
          template <size_t, class, class> class Accumulate,
          template<size_t, class> class Estimate,
          int Dim,
          class SourceSwarm,
          class SourceState,
          class TargetSwarm,
          class TargetState>
void SwarmDriver<Search,
                 Accumulate, 
                 Estimate, 
                 Dim, 
                 SourceSwarm, 
                 SourceState,
                 TargetSwarm, 
                 TargetState>::
remap(std::vector<std::string> const &src_varnames,
      std::vector<std::string> const &trg_varnames,
      bool distributed,
      bool report_time)
{

  int comm_rank = 0;

#ifdef ENABLE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

  int numTargetPts = target_swarm_.num_particles(PARALLEL_OWNED);

  int nvars = source_remap_var_names_.size();
   
  float tot_seconds = 0.0, tot_seconds_dist = 0.0, 
        tot_seconds_srch = 0.0, tot_seconds_xsect = 0.0, 
        tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;

  //DISTRIBUTE
  // This step would change the input source swarm and its state
  // if after distribution it receives particles from other 
  // ranks. 
  // For the scatter scheme, the smoothing_lengths will also 
  // be changed. 
#ifdef ENABLE_MPI
  if (distributed) {
  gettimeofday(&begin_timeval, 0);
  MPI_Particle_Distribute<Dim> distributor;
  
  //For scatter scheme, the smoothing_lengths_, kernel_types_
  //and geom_types_  are also communicated and changed for the
  //source swarm. 
  distributor.distribute(source_swarm_, source_state_,
                         target_swarm_, target_state_,
                         smoothing_lengths_, kernel_types_,
                         geom_types_, weight_center_);
  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_dist = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
}
#endif 

  // SEARCH
  Portage::vector<std::vector<unsigned int>> candidates(numTargetPts);
    
  // Get an instance of the desired search algorithm type which is expected
  // to be a functor with an operator() of the right form

  gettimeofday(&begin_timeval, 0);

  // search boxes around source points are of zero dimension but
  // those around target points are determined by particle
  // smoothing lengths

  // code below does not work with facted weights
  std::shared_ptr<vector<Point<Dim>>> sourceExtents;
  std::shared_ptr<vector<Point<Dim>>> targetExtents;
  targetExtents = std::make_shared<vector<Point<Dim>>>(numTargetPts);
  if (weight_center_ == Portage::Meshfree::Gather) {
    for (int i = 0; i < numTargetPts; i++) {
      if (geom_types_[i] == Weight::FACETED) {
        throw std::runtime_error("FACETED geometry is not available here");
      }
      {Point<Dim> pt=(*targetExtents)[i];
       std::vector<std::vector<double>> vv=smoothing_lengths_[i];
       pt=Point<Dim>(vv[0]); (*targetExtents)[i]=pt;}
    }
  }
  
  int numSourcePts = source_swarm_.num_particles();
  sourceExtents = std::make_shared<vector<Point<Dim>>>(numSourcePts);
  if (weight_center_ == Portage::Meshfree::Scatter) {
    for (int i = 0; i < numSourcePts; i++) {
      if (geom_types_[i] == Weight::FACETED) {
        throw std::runtime_error("FACETED geometry is not available here");
      }
      {Point<Dim> pt=(*sourceExtents)[i];
       std::vector<std::vector<double>> vv=smoothing_lengths_[i];
       pt=Point<Dim>(vv[0]); (*sourceExtents)[i]=pt;}
    }
  }

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
                        basis_type_, 
                        operator_spec_, operator_domains_, operator_data_);

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
  
  nvars = src_varnames.size();
  if (comm_rank == 0)
    std::cout << "number of variables to remap is " << nvars <<
        std::endl;
  
  // Get an instance of the desired interpolate algorithm type
  Estimate<Dim, SourceState> estimateFunctor(source_state_);
  
  for (int i = 0; i < nvars; ++i) {
    //amh: ?? add back accuracy output statement??
    if (comm_rank == 0) std::cout << "Remapping swarm variable " <<
                            src_varnames[i] <<
                            " to variable " << trg_varnames[i] <<
                            std::endl;

    estimateFunctor.set_variable(src_varnames[i]);
    
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
    
    target_state_.get_field(trg_varnames[i], target_field_shared_ptr);
    Portage::pointer<double> target_field(&((*target_field_shared_ptr)[0]));
    
    Portage::transform(target_swarm_.begin(PARTICLE, PARALLEL_OWNED),
                       target_swarm_.end(PARTICLE, PARALLEL_OWNED),
                       source_pts_and_mults.begin(),
                       target_field, estimateFunctor);
    
    
    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    tot_seconds_interp = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
    
    tot_seconds = tot_seconds_dist + tot_seconds_srch +
                  tot_seconds_xsect + tot_seconds_interp;
    
    if (report_time) {
      std::cout << "Swarm Transform Time Rank " << comm_rank << " (s): " <<
        tot_seconds << std::endl;
      std::cout << "  Swarm Distribution Time Rank " << comm_rank << " (s): " <<
        tot_seconds_dist << std::endl;
      std::cout << "  Swarm Search Time Rank " << comm_rank << " (s): " <<
        tot_seconds_srch << std::endl;
      std::cout << "  Swarm Accumulate Time Rank " << comm_rank << " (s): " <<
        tot_seconds_xsect << std::endl;
      std::cout << "  Swarm Estimate Time Rank " << comm_rank << " (s): " <<
        tot_seconds_interp << std::endl;
    }
  }

}//remap


}  // namespace Meshfree
}  // namespace Portage

#endif  // SRC_DRIVER_SWARM_H_
