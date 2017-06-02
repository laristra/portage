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
          template<size_t, class, class, class> class Estimate,
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
              Weight::Kernel const& kernel_type=B4,
              Weight::Geometry const& support_geom_type=ELLIPTIC)
      : source_swarm_(sourceSwarm), source_state_(sourceState),
        target_swarm_(targetSwarm), target_state_(targetState),
    dim_(sourceSwarm.space_dimension(),
    smoothing_lengths_(smoothing_lengths) {

    assert(sourceSwarm.space_dimension() == targetSwarm.space_dimension());
    int num_target_particles = targetSwarm.num_owned_particles();
    kernel_types_ = std::vector<Weight::Kernel>(num_target_particles,
                                                kernel_type);
    geom_types_ = std::vector<Weight::Geometry>(num_target_particles,
                                                support_geom_type);
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
              std::vector<Weight::Geometry> const& geom_types)
      : source_swarm_(sourceSwarm), source_state_(sourceState),
        target_swarm_(targetSwarm), target_state_(targetState),
    dim_(sourceSwarm.space_dimension(),
    kernel_types_(kernel_types),
    geom_types_(support_geom_types),
    smoothing_lengths_(smoothing_lengths) {

    assert(sourceSwarm.space_dimension() == targetSwarm.space_dimension());
    assert(targetSwarm.size() == smoothing_lengths.size());
    assert(targetSwarm.size() == kernel_types.size());
    assert(targetSwarm.size() == support_geom_types.size());
  }

  /// Copy constructor (disabled)
  SwarmDriver(const SwarmDriver &) = delete;

  /// Assignment operator (disabled)
  SwarmDriver & operator = (const Driver &) = delete;

  /// Destructor
  ~SwarmDriver() {}

  /*!
    @brief Specify the names of the variables to be interpolated along with the
    limiter to use for all of them
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source swarm to the target swarm.  This variable must
    exist in both swarmes' state manager
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
      MeshFree::EstimateType const& estimator_type = MeshFree::LocalRegression,
      MeshFree::Basis::Type const& basis = MeshFree::Basis::Unitary
                           ) {
    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    for (int i = 0; i < nvars; ++i)
      assert(source_state_.get_entity(source_remap_var_names[i]) ==
             target_state_.get_entity(target_remap_var_names[i]));

    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    estimator_types_.resize(nvars, estimator_type);
    bases_.resize(nvars, basis);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source swarm.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target swarm.
    @param[in] estimator_types A list indicating what estimators should be 
    used for each remapped variable (KernelDensity, LocalRegression)
    @param[in] bases A list indicating what basis order to use for each variable
    (UNITARY, LINEAR, QUADRATIC)
  */
  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      std::vector<MeshFree::EstimateType> const& estimator_types,
      std::vector<MeshFree::Basis> const& bases) {

    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    for (int i = 0; i < nvars; ++i)
      assert(source_state_.get_entity(source_remap_var_names[i]) ==
             target_state_.get_entity(target_remap_var_names[i]));

    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    estimator_types_ = estimator_types;
    bases_ = bases;
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
    @brief Get the dimensionality of the swarmes.
    @return The dimensionality of the swarmes.
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

    // EVEN OTHERWISE WE DON'T HAVE THE ABILITY TO DO DISTRIBUTED SWARMS
    if (distributed) {
      std::cout << "CANNOT DO DISTRIBUTED REMAPPING WITH SWARMS\n";
      return;
    }

    int numTargetPts = target_swarm_.num_owned_particles();
    std::cout << "Number of target cells in target swarm on rank "
              << comm_rank << ": "
              << numTargetPts << std::endl;

    int nvars = source_remap_var_names_.size();

    // Collect all variables and remap them

    std::vector<std::string> source_var_names;
    std::vector<std::string> target_var_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == PARTICLE) {
        source_var_names.emplace_back(source_remap_var_names_[i]);
        target_var_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_var_names.size() > 0)
    {
      
      float tot_seconds = 0.0, tot_seconds_srch = 0.0,
          tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
      struct timeval begin_timeval, end_timeval, diff_timeval;
      
      // SEARCH
      Portage::vector<std::vector<int>> candidates(numTargetPts);
      
      // Get an instance of the desired search algorithm type which is expected
      // to be a functor with an operator() of the right form

      gettimeofday(&begin_timeval, 0);
      const Search<Dim, SourceSwarm, TargetSwarm>
          searchfunctor(source_swarm_, target_swarm_);
      
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
                            estimate_type_, weight_center_,
                            kernel_types_, geom_types_, smoothing_lengths_,
                            basis_type_);
    
      Portage::vector<std::vector<Weights_t>> source_pts_and_mults(nTargetPts);
      
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

        estimateFunctor.set_variable(source_var_names[i], limiters_[i]);
        
        // This populates targetField with the values returned by the
        // remapper operator
        
        double *target_field_raw = nullptr;
        target_state_.get_data(PARTICLE, target_nodevar_names[i],
                               &target_field_raw);
        Portage::pointer<double> target_field(target_field_raw);
        
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
  std::vector<MeshFree::EstimateType> estimator_types_;
  std::vector<MeshFree::Weight::Kernel> kernel_types_;
  std::vector<MeshFree::Weight::Geometry> support_geom_types_;
  std::vector<MeshFree::Basis> bases_;
  std::vector<std::vector<std::vector<double>>> smoothing_lengths_;
  unsigned int dim_;
};  // class Driver_Swarm


}  // namespace Portage

#endif  // SRC_DRIVER_DRIVER_H_
