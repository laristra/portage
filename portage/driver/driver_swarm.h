/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 *  https://github.com/laristra/portage/blob/master/LICENSE
 */
#ifndef SRC_DRIVER_SWARM_H_
#define SRC_DRIVER_SWARM_H_

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <cmath>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "portage/support/portage.h"
#include "portage/support/timer.h"
#include "portage/support/basis.h"
#include "portage/support/weight.h"
#include "portage/support/operator.h"
#include "portage/search/search_simple_points.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

#ifdef WONTON_ENABLE_MPI
  #include "portage/distributed/mpi_particle_distribute.h"
#endif

namespace Portage { namespace Meshfree {
  // avoid very long type names.
  using SmoothingLengths = Wonton::vector<std::vector<std::vector<double>>>;

/**
 * @brief Provides an interface to remap variables from one swarm to another.
 *
 * @tparam Search: the particle search type to use.
 * @tparam Accumulate: the particle accumulate type to use.
 * @tparam Estimate: the particle estimate type to use.
 * @tparam dim: spatial dimension of the problem.
 * @tparam SourceSwarm: the wrapper type of the input source swarm.
 * @tparam SourceState: the wrapper type of the input source state.
 * @tparam TargetSwarm: the wrapper type of the input target swarm.
 * @tparam TargetState: the wrapper type of the input target state.
 */
template <template <int, class, class> class Search,
          template <int, class, class> class Accumulate,
          template<int, class> class Estimate,
          int dim,
          class SourceSwarm, class SourceState,
          class TargetSwarm = SourceSwarm,
          class TargetState = SourceState>
class SwarmDriver {
public:
  /**
   * @brief Constructor using a unique kernel_type and support geometry for all particles.
   *
   * @param source_swarm: a reference to the source swarm
   * @param source_state: a reference to the source state
   * @param target_swarm: a reference to the target swarm
   * @param target_state: a reference to the target state
   * @param smoothing_lengths: a vector of smoothing lengths for each
   *                           target particle - the three levels of vectors
   *                           allow for a polygonal/polyhedral shaped support
   *                           around each target particle and a vector of
   *                           smoothing lengths for each facet of the
   *                           polygonal/polyhedral support.
   * @param kernel_types: types of weight kernels (B4, SQUARE, etc) to use.
   * @param geom_types: geometry of the supports (ELLIPTIC, TENSOR, FACETED)
   *                    of the weight functions to use.
   * @param center: specify whether gather-form or scatter-form weights are used.
   *
   * NOTE: 'smoothing_lengths' must have the size of 'source_swarm' if 'scatter'
   *       and that of 'target_swarm' if 'gather'.
   */
  SwarmDriver(SourceSwarm& source_swarm,
              SourceState& source_state,
              TargetSwarm const& target_swarm,
              TargetState& target_state,
              SmoothingLengths const& smoothing_lengths,
              Weight::Kernel const kernel_type = Weight::B4,
              Weight::Geometry const support_geom_type = Weight::ELLIPTIC,
              WeightCenter const center=Gather)
      : source_swarm_(source_swarm),
        target_swarm_(target_swarm),
        source_state_(source_state),
        target_state_(target_state),
        smoothing_lengths_(smoothing_lengths) {

    assert(dim == source_swarm.space_dimension());
    assert(dim == target_swarm.space_dimension());
    weight_center_ = center;
    check_sizes_and_set_types(center, kernel_type, support_geom_type);
    set_extents_from_smoothing_lengths();
  }

  /**
   * @brief Constructor with search extents computed inline.
   *
   *  This is used for any weights except faceted.
   *
   * @param source_swarm: a reference to the source swarm
   * @param source_state: a reference to the source state
   * @param target_swarm: a reference to the target swarm
   * @param target_state: a reference to the target state
   * @param smoothing_lengths: a vector of smoothing lengths for each
   *                           target particle - the three levels of vectors
   *                           allow for a polygonal/polyhedral shaped support
   *                           around each target particle and a vector of
   *                           smoothing lengths for each facet of the
   *                           polygonal/polyhedral support.
   * @param kernel_types: types of weight kernels (B4, SQUARE, etc) to use.
   * @param geom_types: geometry of the supports (ELLIPTIC, TENSOR, FACETED)
   *                    of the weight functions to use.
   * @param center: specify whether gather-form or scatter-form weights are used.
   *
   * NOTE: 'smoothing_lengths' must have the size of 'source_swarm' if 'scatter'
   *       and that of 'target_swarm' if 'gather'.
   */
  SwarmDriver(SourceSwarm& source_swarm,
              SourceState& source_state,
              TargetSwarm const& target_swarm,
              TargetState& target_state,
              SmoothingLengths const& smoothing_lengths,
              Wonton::vector<Weight::Kernel> const& kernel_types,
              Wonton::vector<Weight::Geometry> const& geom_types,
              WeightCenter const center = Gather)
      : source_swarm_(source_swarm),
        target_swarm_(target_swarm),
        source_state_(source_state),
        target_state_(target_state),
        smoothing_lengths_(smoothing_lengths),
        kernel_types_(kernel_types),
        geom_types_(geom_types) {

    assert(dim == source_swarm.space_dimension());
    assert(dim == target_swarm.space_dimension());
    weight_center_ = center;
    check_sizes(center);
    set_extents_from_smoothing_lengths();
  }

  /**
   * @brief Constructor with search extents specified.
   *
   * This is used for faceted weights, which will be uniformly used throughout
   * the swarm specified by center. Kernel and geometry type are set internally.
   *
   * @param source_swarm: a reference to the source swarm
   * @param source_state: a reference to the source state
   * @param target_swarm: a reference to the target swarm
   * @param target_state: a reference to the target state
   * @param smoothing_lengths: a vector of smoothing lengths for each
   *                           target particle - the three levels of vectors
   *                           allow for a polygonal/polyhedral shaped support
   *                           around each target particle and a vector of
   *                           smoothing lengths for each facet of the
   *                           polygonal/polyhedral support.
   * @param source_extents: extents for source swarm
   * @param target_extents: extents for target swarm
   * @param center: specify whether gather-form or scatter-form weights are used.
   *
   * NOTE: 'smoothing_lengths' must have the size of 'source_swarm' if 'scatter'
   *       and that of 'target_swarm' if 'gather'.
   */
  SwarmDriver(SourceSwarm& source_swarm,
              SourceState& source_state,
              TargetSwarm const& target_swarm,
              TargetState& target_state,
              SmoothingLengths const& smoothing_lengths,
              Wonton::vector<Point<dim>> const& source_extents,
              Wonton::vector<Point<dim>> const& target_extents,
              WeightCenter const center = Gather)
      : source_swarm_(source_swarm),
        target_swarm_(target_swarm),
        source_state_(source_state),
        target_state_(target_state),
        smoothing_lengths_(smoothing_lengths) {
    
   assert(dim == source_swarm.space_dimension());
   assert(dim == target_swarm.space_dimension());
   weight_center_ = center;

   int const swarm_size = get_swarm_size();

#ifdef DEBUG
   int const nb_source = source_swarm_.num_particles(Wonton::PARALLEL_OWNED);
   int const nb_target = target_swarm_.num_particles(Wonton::PARALLEL_OWNED);

   switch (weight_center_) {
     case Gather:  assert(target_extents.size() == unsigned(nb_target)); break;
     case Scatter: assert(source_extents.size() == unsigned(nb_source)); break;
     default: throw std::runtime_error("invalid weight center type");
   }

   assert(smoothing_lengths.size() == unsigned(swarm_size));
#endif

    kernel_types_.resize(swarm_size, Weight::POLYRAMP);
    geom_types_.resize(swarm_size, Weight::FACETED);
    source_extents_ = source_extents;
    target_extents_ = target_extents;
  }

  /**
   * @brief Disabled copy constructor
   */
  SwarmDriver(const SwarmDriver&) = delete;

  /**
   * @brief Disabled assignment operator
   */
  SwarmDriver& operator = (const SwarmDriver &) = delete;

  /**
   * @brief Destructor
   */
  ~SwarmDriver() = default;

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source swarm to the target swarm.  This variable must
    exist in both swarms' state manager
  */
  void set_remap_var_names(std::vector<std::string> const& remap_var_names) {
    // remap variable names same in source and target swarm
    set_remap_var_names(remap_var_names, remap_var_names);
  }

  /**
   * @brief Specify the names of the variables to be interpolated.
   *
   * @param source_vars: a list of source swarm variables names.
   * @param target_vars: a list of target swarm variables names.
   * @param estimator_type: estimator to be used (KernelDensity, LocalRegression).
   * @param basis_type: order of the basis used (UNITARY, LINEAR, QUADRATIC).
   * @param operator_spec: operator specification.
   * @param operator_domains: operator domains.
   * @param operator_data: operator data.
   * @param part_field: name of the field for part-by-particle.
   * @param part_tolerance: tolerance threshold for part-by-particle.
   * @param part_smoothing: smoothing lengths matrix for part-by-particle.
   */
  void set_remap_var_names(std::vector<std::string> const& source_vars,
                           std::vector<std::string> const& target_vars,
                           EstimateType const estimator_type = LocalRegression,
                           basis::Type const basis_type = basis::Unitary,
                           oper::Type const operator_spec = oper::LastOperator,
                           Wonton::vector<oper::Domain> const& operator_domains = {},
                           Wonton::vector<std::vector<Point<dim>>> const& operator_data = {},
                           std::string part_field = "NONE",
                           double part_tolerance = 0.0,
                           SmoothingLengths const& part_smoothing = {}) {

    assert(source_vars.size() == target_vars.size());
    // variables names
    source_vars_ = source_vars;
    target_vars_ = target_vars;

    // estimator and operator
    estimator_type_   = estimator_type;
    basis_type_       = basis_type;
    operator_spec_    = operator_spec;
    operator_domains_ = operator_domains;
    operator_data_    = operator_data;
#ifdef DEBUG
    if (operator_spec_ != oper::LastOperator) {
      unsigned const num_target_particles = target_swarm_.num_owned_particles();
      assert(operator_domains_.size() == num_target_particles);
      assert(operator_data_.size() == num_target_particles);
    }
#endif
    // part-by-particle
    part_field_     = part_field;
    part_tolerance_ = part_tolerance;
    part_smoothing_ = part_smoothing;
  }

  /**
   * @brief Get all source swarm variables names.
   *
   * @return a list of source variables names.
   */
  std::vector<std::string> source_vars() const { return source_vars_; }

  /**
   * @brief Get all target swarm variables names.
   *
   * @return a list of target variables names.
   */
  std::vector<std::string> target_vars() const { return target_vars_; }

  /**
   * @brief Perform the remap.
   *
   * @param executor: the executor type: serial or parallel.
   * @param report_time: specify if saving timings or not.
   */
  void run(Wonton::Executor_type const* executor = nullptr, bool report_time = true) {

    int rank   = 0;

#ifdef WONTON_ENABLE_MPI
    bool distributed = false;
    MPI_Comm comm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const*>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      comm = mpiexecutor->mpicomm;
      MPI_Comm_rank(comm, &rank);
      int nprocs = 0;
      MPI_Comm_size(comm, &nprocs);
      distributed = nprocs > 1;
    }
#endif

    if (rank == 0)
      std::cout << "in SwarmDriver::run() ... " << std::endl;

    // useful aliases
    using Searcher = Search<dim, SourceSwarm, TargetSwarm>;
    using Accumulator = Accumulate<dim, SourceSwarm, TargetSwarm>;
    using Estimator = Estimate<dim, SourceState>;

    int nb_source = source_swarm_.num_particles(Wonton::PARALLEL_OWNED);
    int nb_target = target_swarm_.num_particles(Wonton::PARALLEL_OWNED);
    int nb_fields = source_vars_.size();

    float tot_seconds = 0.0, tot_seconds_dist = 0.0,
      tot_seconds_srch = 0.0, tot_seconds_xsect = 0.0,
      tot_seconds_interp = 0.0;
    //struct timeval begin_timeval, end_timeval, diff_timeval;
    auto tic = timer::now();

    //DISTRIBUTE
    // This step would change the input source swarm and its state
    // if after distribution it receives particles from other
    // ranks.
    // For the scatter scheme, the smoothing_lengths will also
    // be changed.
#ifdef WONTON_ENABLE_MPI
    if (distributed) {
      MPI_Particle_Distribute<dim> distributor(mpiexecutor);
      //For scatter scheme, the smoothing_lengths_, kernel_types_
      //and geom_types_  are also communicated and changed for the
      //source swarm.
      distributor.distribute(source_swarm_, source_state_,
                             target_swarm_, target_state_,
                             smoothing_lengths_, source_extents_, target_extents_,
                             kernel_types_, geom_types_, weight_center_);

      tot_seconds_dist = timer::elapsed(tic, true);
    }
#endif

    // SEARCH
    Wonton::vector<std::vector<int>> candidates(nb_target);

    // Get an instance of the desired search algorithm type which is expected
    // to be a functor with an operator() of the right form
    Searcher search(source_swarm_, target_swarm_,
                    source_extents_, target_extents_, weight_center_);

    Wonton::transform(target_swarm_.begin(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                       target_swarm_.end(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                       candidates.begin(), search);

    tot_seconds_srch = timer::elapsed(tic);

    // In case of faceted, scatter, and parts, obtain part assignments on target particles.
    // Eliminate neighbors that aren't in the same part.
    // It is assumed that the faceted weight function will be non-zero only on the
    // source cell it came from, which is achieved by using a smoothing factor of 1/2.
    // It is also assumed no target point will appear in more than one source cell.
    if (geom_types_[0] == Weight::FACETED and weight_center_ == Scatter and part_field_!="NONE") {


      // make storage for target part assignments
      nb_source = source_swarm_.num_particles(Wonton::PARALLEL_OWNED);
      nb_target = target_swarm_.num_particles(Wonton::PARALLEL_OWNED);

      // get source part assignments
      auto source_field_part = source_state_.get_field_dbl(part_field_);
      std::vector<double> target_field_part(nb_target);

      // create accumulator to evaluate weight function on source cells
      Wonton::vector<Weight::Kernel> step_kern(nb_source, Weight::STEP);
      Accumulator accumulator(source_swarm_, target_swarm_,
                              estimator_type_, weight_center_,
                              step_kern, geom_types_, part_smoothing_, basis_type_,
                              operator_spec_, operator_domains_, operator_data_);

      // loop over target points and set part assignments

      for (int i = 0; i < nb_target; i++) {
        std::vector<int> possibles = candidates[i];
        for (auto&& neighbor : possibles) {
          double weight = accumulator.weight(i, neighbor);
          if (weight > 0.) {
            double part_val = source_field_part[neighbor];
            target_field_part[i] = part_val;
#undef DEBUG_HERE
#ifdef DEBUG_HERE
            if (dim == 2) {
              auto p1 = source_swarm_.get_particle_coordinates(nbr);
              auto p2 = target_swarm_.get_particle_coordinates(i);
              double h = part_smoothing_[neighbor][0][2];

              std::cout << "i=" << i << " " << " neighbor=" << neighbor << " ";
              std::cout << std::abs(p2[0]-p1[0])/h <<", "<< std::abs(p2[1]-p1[1])/h;
              std::cout << " w= "<< weight <<" part= "<< part_val << std::endl;
            }
#endif
#undef DEBUG_HERE
          }
        }
      }

      // loop over target points and dismiss neighbors that aren't in the same part
      for (int i = 0; i < nb_target; i++) {
        std::vector<int> new_candidates;
        std::vector<int> possibles = candidates[i];
        for (auto&& neighbor : possibles) {
          if (std::abs(target_field_part[i] - source_field_part[neighbor]) < part_tolerance_) {
            new_candidates.emplace_back(neighbor);
          }
        }
        candidates[i] = new_candidates;
      }
    }

    // ACCUMULATE (build moment matrix, calculate shape functions)
    // EQUIVALENT TO INTERSECT IN MESH-MESH REMAP
    tic = timer::now();

    // Get an instance of the desired accumulate algorithm type which is
    // expected to be a functor with an operator() of the right form
    Accumulator accumulate(source_swarm_, target_swarm_,
                           estimator_type_, weight_center_, kernel_types_,
                           geom_types_, smoothing_lengths_, basis_type_,
                           operator_spec_, operator_domains_, operator_data_);

    Wonton::vector<std::vector<Weights_t>> source_points_and_multipliers(nb_target);

    // For each particle in the target swarm get the shape functions
    // (multipliers for source particle values)
    Wonton::transform(target_swarm_.begin(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                       target_swarm_.end(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                       candidates.begin(), source_points_and_multipliers.begin(),
                       accumulate);

    tot_seconds_xsect = timer::elapsed(tic, true);

    // ESTIMATE (one variable at a time)
    nb_fields = source_vars_.size();
    if (rank == 0)
      std::cout << "number of variables to remap is " << nb_fields << std::endl;

    // Get an instance of the desired interpolate algorithm type
    Estimator estimator(source_state_);

    for (int i = 0; i < nb_fields; ++i) {
      //amh: ?? add back accuracy output statement??
      if (rank == 0)
        std::cout << "Remap "<< source_vars_[i] <<" to "<< target_vars_[i] << std::endl;

      estimator.set_variable(source_vars_[i]);

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

      // TODO: perform a deep-copy back to target state
      auto& target_data = target_state_.get_field(target_vars_[i]);
      Wonton::pointer<double> target_field(target_data.data());

      Wonton::transform(target_swarm_.begin(Entity_kind::PARTICLE, Entity_type::PARALLEL_OWNED),
                         target_swarm_.end(Entity_kind::PARTICLE, Entity_type::PARALLEL_OWNED),
                         source_points_and_multipliers.begin(),
                         target_field, estimator);

      tot_seconds_interp = timer::elapsed(tic, true);
      tot_seconds = tot_seconds_dist + tot_seconds_srch +
                    tot_seconds_xsect + tot_seconds_interp;

      if (report_time) {
        std::cout << "Swarm Transform Time Rank " << rank << " (s): " << tot_seconds << std::endl;
        std::cout << "  Swarm Distribution Time Rank " << rank << " (s): " << tot_seconds_dist << std::endl;
        std::cout << "  Swarm Search Time Rank " << rank << " (s): " << tot_seconds_srch << std::endl;
        std::cout << "  Swarm Accumulate Time Rank " << rank << " (s): " << tot_seconds_xsect << std::endl;
        std::cout << "  Swarm Estimate Time Rank " << rank << " (s): " << tot_seconds_interp << std::endl;

        // put out neighbor statistics
        int nnbrmax = 0;
        int nnbrmin = nb_target;
        int nnbravg = 0;
        int nnbrsum = 0;
        double nnbrsdev = 0;

        for (int j=0; j < nb_target; j++) {
          int n = 0;
          std::vector<Weights_t> wts = source_points_and_multipliers[j];
          for (auto&& wt : wts) {
            if (std::abs(wt.weights[0]) > 0.0)
              n++;
          }
          if (n > nnbrmax) nnbrmax = n;
          if (n < nnbrmin) nnbrmin = n;
          nnbrsum += n;
        }
        nnbravg = nnbrsum / nb_target;

        for (int j=0; j < nb_target; j++) {
          int n = 0;
          std::vector<Weights_t> wts = source_points_and_multipliers[j];
          for (auto&& wt : wts) {
            if (std::abs(wt.weights[0]) > 0.0)
              n++;
          }
          nnbrsdev += std::pow(n - nnbravg, 2);
        }
        nnbrsdev = std::sqrt(nnbrsdev / nb_target);

        std::cout << "Max number of neighbors: " << nnbrmax << std::endl;
        std::cout << "Min number of neighbors: " << nnbrmin << std::endl;
        std::cout << "Avg number of neighbors: " << nnbravg << std::endl;
        std::cout << "Std Dev for number of neighbors: " << nnbrsdev << std::endl;
      }
    }

  }  // run

protected:
  /**
   * @brief Get swarm size according to weight center type.
   *
   * @return the number of particles of the swarm.
   */
  int get_swarm_size() const {
    return (weight_center_ == Scatter
      ? source_swarm_.num_particles(Wonton::PARALLEL_OWNED)
      : target_swarm_.num_particles(Wonton::PARALLEL_OWNED));
  }

  /**
   * @brief Check sizes according to weight center type.
   *
   * @param weight_center: weight center type.
   */
  void check_sizes(WeightCenter const weight_center) {
#ifdef DEBUG
    unsigned const swarm_size = get_swarm_size();
    assert(smoothing_lengths_.size() == swarm_size);
    assert(kernel_types_.size() == swarm_size);
    assert(geom_types_.size() == swarm_size);
#endif
  }

  /**
   * @brief Check sizes according to weight center, set kernel and geometry.
   *
   * @param weight_center: weight center type.
   * @param kernel_type: kernel type
   * @param support_geom_type: support geometry type.
   */
  void check_sizes_and_set_types(WeightCenter const weight_center,
                                 Weight::Kernel const kernel_type,
                                 Weight::Geometry const support_geom_type) {
    int const swarm_size = get_swarm_size();
    assert(smoothing_lengths_.size() == unsigned(swarm_size));
    kernel_types_.resize(swarm_size, kernel_type);
    geom_types_.resize(swarm_size, support_geom_type);
  }

  /**
   * @brief Setup search boxes around source and target points
   *        for non-faceted weights.
   */
  void set_extents_from_smoothing_lengths() {
    if (weight_center_ == Gather) {
      int const nb_target = target_swarm_.num_particles(Wonton::PARALLEL_OWNED);
      assert(smoothing_lengths_.size() == unsigned(nb_target));
      target_extents_.resize(nb_target);

      for (int i = 0; i < nb_target; i++) {
        if (geom_types_[i] == Weight::FACETED) {
          throw std::runtime_error("FACETED geometry is not available here");
        }
        std::vector<std::vector<double>> temp = smoothing_lengths_[i];
        target_extents_[i] = Point<dim>(temp[0]);
      }
    } else if (weight_center_ == Scatter) {
      int const nb_source = source_swarm_.num_particles(Wonton::PARALLEL_OWNED);
      assert(smoothing_lengths_.size() == unsigned(nb_source));
      source_extents_.resize(nb_source);

      for (int i = 0; i < nb_source; i++) {
        if (geom_types_[i] == Weight::FACETED) {
          throw std::runtime_error("FACETED geometry is not available here");
        }
        std::vector<std::vector<double>> temp = smoothing_lengths_[i];
        source_extents_[i] = Point<dim>(temp[0]);
      }
    }
  }

private:
  SourceSwarm& source_swarm_;
  TargetSwarm const& target_swarm_;
  SourceState& source_state_;
  TargetState& target_state_;
  std::vector<std::string> source_vars_ {};
  std::vector<std::string> target_vars_ {};
  WeightCenter weight_center_ = Gather;
  SmoothingLengths smoothing_lengths_ {};
  Wonton::vector<Weight::Kernel> kernel_types_ {};
  Wonton::vector<Weight::Geometry> geom_types_ {};
  Wonton::vector<Point<dim>> source_extents_ {};
  Wonton::vector<Point<dim>> target_extents_ {};
  EstimateType estimator_type_ {};
  basis::Type basis_type_ {};
  oper::Type operator_spec_ {};
  Wonton::vector<oper::Domain> operator_domains_ {};
  Wonton::vector<std::vector<Point<dim>>> operator_data_ {};
  std::string part_field_ = "";
  double part_tolerance_ = 0.0;
  Wonton::vector<std::vector<std::vector<double>>> part_smoothing_ {};
};  // class SwarmDriver

}}  // namespace Portage::Meshfree

#endif  // SRC_DRIVER_SWARM_H_
