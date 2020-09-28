/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_DRIVER_MESH_SWARM_MESH_H_
#define SRC_DRIVER_MESH_SWARM_MESH_H_

//#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <string>
#include <limits>

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/swarm/swarm.h"
#include "wonton/swarm/swarm_state.h"

// portage includes
#include "portage/support/timer.h"
#include "portage/support/portage.h"
#include "portage/support/basis.h"
#include "portage/support/weight.h"
#include "portage/support/operator.h"
#include "portage/support/faceted_setup.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/driver/driver_swarm.h"

#ifdef WONTON_ENABLE_MPI
  #include "portage/distributed/mpi_bounding_boxes.h"
#endif

/*!
  @file driver_mesh_swarm_mesh.h
  @brief Example driver for mapping between two meshes with a particle swarm
  intermediary. (Mesh-Swarm-Mesh = MSM).

  This should serve as a good example for how to write your own driver routine
  and datastructures.
*/

namespace Portage {

/*!
  @class MSM_Driver "driver_mesh_swarm_mesh.h"
  @brief MSM_Driver provides the API to mapping from one mesh to another using particles.
  @tparam Search compatible particle search class
  @tparam Accumulate Meshfree accumulate class (performs particle sums)
  @tparam Estimate Meshfree estimator class (performs calculus operations
  @tparam dim Spatial dimension of the source and target meshes
  @tparam SourceMesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.
  @tparam SourceState_Wrapper A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.
  @tparam TargetMesh_Wrapper A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.
  @tparam TargetState_Wrapper A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.
*/
template <template <int, class, class> class Search,
          template<int, class, class> class Accumulate,
          template<int, class> class Estimate,
          int dim,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper = SourceMesh_Wrapper,
          class TargetState_Wrapper = SourceState_Wrapper>
class MSM_Driver {
public:
  /*!
    @brief Constructor for running the interpolation driver.
    @param[in] source_mesh A @c SourceMesh_Wrapper to the source mesh.
    @param[in] source_state A @c SourceState_Wrapperfor the data that lives on the
    source mesh.
    @param[in] target_mesh A @c TargetMesh_Wrapper to the target mesh.
    @param[in,out] target_state A @c TargetState_Wrapper for the data that will
    be mapped to the target mesh.
    @param[in] smoothing_factor multiplies cell sizes to get smoothing lengths
    @param[in] boundary_factor with faceted weights only, multiplies center-to-face distance on boundary to get smoothing length
    @param[in] geometry weight function geometric configuration
    @param[in] kernel weight function kernel
    @param[in] center weight function centering (scatter for sources, gather for targets)
    @param[in] part_field name of field to use for part assignments, with faceted weights only
    @param[in] part_tolerance tolerance for determining if part assignment matches a neighbor's assignment
  */
  MSM_Driver(SourceMesh_Wrapper const& source_mesh,
             SourceState_Wrapper const& source_state,
             TargetMesh_Wrapper const& target_mesh,
             TargetState_Wrapper& target_state,
             double smoothing_factor = 1.25,
             double boundary_factor  = 0.5,
             Meshfree::Weight::Geometry geometry = Meshfree::Weight::TENSOR,
             Meshfree::Weight::Kernel   kernel   = Meshfree::Weight::B4,
             Meshfree::WeightCenter     center   = Meshfree::Gather,
             std::string part_field = "NONE",
             double part_tolerance = std::numeric_limits<double>::infinity())

      : source_mesh_(source_mesh),
        target_mesh_(target_mesh),
        source_state_(source_state),
        target_state_(target_state),
        smoothing_factor_(smoothing_factor),
        boundary_factor_(boundary_factor),
        kernel_(kernel),
        geometry_(geometry),
        center_(center),
        part_field_(part_field), 
        part_tolerance_(part_tolerance),
        dim_(source_mesh.space_dimension()) {

    assert(source_mesh.space_dimension() == target_mesh.space_dimension());
    assert(source_mesh.space_dimension() == dim);
    if (geometry == Meshfree::Weight::FACETED) {
      assert(kernel == Meshfree::Weight::POLYRAMP);
    }
  }

  /// Copy constructor (disabled)
  MSM_Driver(const MSM_Driver &) = delete;

  /// Assignment operator (disabled)
  MSM_Driver& operator = (const MSM_Driver &) = delete;

  /// Destructor
  ~MSM_Driver() = default;

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source mesh to the target mesh.  This variable must
    exist in both meshes' state manager
  */
  void set_remap_var_names(std::vector<std::string> const& remap_var_names) {
    // remap variable names same in source and target mesh
    set_remap_var_names(remap_var_names, remap_var_names);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_vars A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] target_vars  A list of the variables names of the
    variables to interpolate to the target mesh.
    @param[in] estimator_type What type of estimator to apply in particle stage
    @param[in] basis_type what regression basis to use
    @param[in] operator_spec what type of operator to use in regression, if any
    @param[in] operator_domains if using an integral operator, what domains of integration to use
    @param[in] operator_data node data for integral domains, if needed
  */

  void set_remap_var_names(std::vector<std::string> const& source_vars,
                           std::vector<std::string> const& target_vars,
                           Meshfree::EstimateType const& estimator_type = Meshfree::LocalRegression,
                           Meshfree::basis::Type const& basis_type = Meshfree::basis::Unitary,
                           Meshfree::oper::Type operator_spec = Meshfree::oper::LastOperator,
                           Wonton::vector<Meshfree::oper::Domain> const& operator_domains = {},
                           Wonton::vector<std::vector<Point<dim>>> const& operator_data = {}) {
#ifndef NDEBUG
    assert(source_vars.size() == target_vars.size());

    int nvars = source_vars.size();
    for (int i = 0; i < nvars; ++i) {
      auto const& source_entity = source_state_.get_entity(source_vars[i]);
      auto const& target_entity = target_state_.get_entity(target_vars[i]);
      assert(source_entity == target_entity);
    }
#endif

    source_vars_      = source_vars;
    target_vars_      = target_vars;
    estimate_         = estimator_type;
    basis_            = basis_type;
    operator_spec_    = operator_spec;
    operator_domains_ = operator_domains;
    operator_data_    = operator_data;
  }

  /*!
    @brief Get the names of the variables to be remapped from the
    source mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> source_vars() const { return source_vars_; }

  /*!
    @brief Get the names of the variables to be remapped to the
    target mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> target_vars() const { return target_vars_; }

  /*!
    @brief Get the dimensionality of the meshes.
    @return The dimensionality of the meshes.
  */
  int dimension() const { return dim_; }

  /*!
    @brief Execute the remapping process
  */
  void run(Wonton::Executor_type const *executor = nullptr) {

    using namespace Meshfree;
    // useful aliases
    using SwarmRemap = SwarmDriver<Search, Accumulate, Estimate, dim,
                                   Wonton::Swarm<dim>, Wonton::SwarmState<dim>>;

    int rank = 0;
    int nprocs = 1;
    
#ifdef WONTON_ENABLE_MPI
    MPI_Comm mycomm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm, &rank);
      MPI_Comm_size(mycomm, &nprocs);
      if (nprocs > 1) {
#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
        std::cerr << "Cannot run Mesh-Swarm-Mesh driver in distributed mode yet";
        std::cerr << std::endl;
#endif
        return;
      }
    }
#endif
#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
    if (rank == 0)
      std::cout << "in MSM_Driver::run() ... " << std::endl;

    int ncells = target_mesh_.num_owned_cells();
    std::cout << "Number of owned target cells on rank "<< rank <<": "<< ncells;
    std::cout << std::endl;

    auto tic = timer::now();
#endif
    int nvars = source_vars_.size();

    // Wonton::CELL VARIABLE SECTION ---------------------------------------------------

    // get cell variable names
    std::vector<std::string> source_cellvar_names;
    std::vector<std::string> target_cellvar_names;

    for (int i = 0; i < nvars; ++i) {
      auto onwhat = source_state_.get_entity(source_vars_[i]);
      if (onwhat == Wonton::CELL) {
        source_cellvar_names.emplace_back(source_vars_[i]);
        target_cellvar_names.emplace_back(target_vars_[i]);
      }
    }

    // Collect all cell based variables and remap them
    if (not source_cellvar_names.empty()) {
      // convert mesh and state wrappers to swarm ones
      Wonton::Swarm<dim> source_swarm(source_mesh_, Wonton::CELL);
      Wonton::Swarm<dim> target_swarm(target_mesh_, Wonton::CELL);
      Wonton::SwarmState<dim> source_swarm_state(source_state_, Wonton::CELL);
      Wonton::SwarmState<dim> target_swarm_state(target_state_, Wonton::CELL);

      // set up smoothing lengths and extents
      Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths;
      std::vector<std::vector<double>> default_lengths(1, std::vector<double>(dim));
      Wonton::vector<Wonton::Point<dim>> weight_extents, other_extents;
      Wonton::vector<std::vector<std::vector<double>>> part_smoothing; // only for faceted,scatter,parts

      if (geometry_ == Weight::FACETED) {
        using Weight::faceted_setup_cell;

        if (part_field_ == "NONE") {
          switch (center_) {
            case Scatter: faceted_setup_cell<dim>(source_mesh_,
                                                  smoothing_lengths,
                                                  weight_extents,
                                                  smoothing_factor_,
                                                  boundary_factor_); break;
            case Gather:  faceted_setup_cell<dim>(target_mesh_,
                                                  smoothing_lengths,
                                                  weight_extents,
                                                  smoothing_factor_,
                                                  boundary_factor_); break;
            default: break;
          }
        } else {
          Wonton::vector<Wonton::Point<dim>> dummy_extents;
          switch (center_) {
            case Scatter:
              // Set smoothing_factor to 1/4 to make weight support exactly equal to cell volume.
              // Store these smoothing lengths off on the side for later use in the swarm driver.
              faceted_setup_cell<dim>(source_mesh_, part_smoothing,
                                      dummy_extents, 0.25, 0.25);
              // Get usual smoothing lengths and extents
              faceted_setup_cell<dim>(source_mesh_,
                                      smoothing_lengths, weight_extents,
                                      smoothing_factor_, boundary_factor_);
              break;
            case Gather: faceted_setup_cell<dim>(target_mesh_, target_state_,
                                                 part_field_, part_tolerance_,
                                                 smoothing_lengths, weight_extents,
                                                 smoothing_factor_, boundary_factor_);
            break;
            default: break;
          }
        }
      } else /* part_field_ != NONE */ {
        int ncells = (center_ == Scatter ? source_mesh_.num_owned_cells()
                                         : target_mesh_.num_owned_cells());

        smoothing_lengths.resize(ncells, default_lengths);

        for (int i = 0; i < ncells; i++) {
          double radius = 0.0;
          switch (center_) {
            case Scatter: Wonton::cell_radius<dim>(source_mesh_, i, &radius); break;
            case Gather:  Wonton::cell_radius<dim>(target_mesh_, i, &radius); break;
            default: break;
          }

          std::vector<std::vector<double>> h = smoothing_lengths[i];
          h[0] = std::vector<double>(dim, 2. * radius * smoothing_factor_);
          smoothing_lengths[i] = h;
        }
      }

      // create swarm remap driver
      SwarmRemap* swarm_remap_ptr = nullptr;

      if (geometry_ == Weight::FACETED) {
        if (center_ == Scatter) {
          swarm_remap_ptr = new SwarmRemap(source_swarm, source_swarm_state,
                                           target_swarm, target_swarm_state,
                                           smoothing_lengths,
                                           weight_extents, other_extents,
                                           center_);
        } else if (center_ == Gather) {
          swarm_remap_ptr = new SwarmRemap(source_swarm, source_swarm_state,
                                           target_swarm, target_swarm_state,
                                           smoothing_lengths,
                                           other_extents, weight_extents,
                                           center_);
        }
      } else {
        swarm_remap_ptr = new SwarmRemap(source_swarm, source_swarm_state,
                                         target_swarm, target_swarm_state,
                                         smoothing_lengths,
                                         kernel_, geometry_, center_);
      }

      assert(swarm_remap_ptr != nullptr);
      SwarmRemap& swarm_remap(*swarm_remap_ptr);
      swarm_remap.set_remap_var_names(source_cellvar_names, target_cellvar_names,
                                      estimate_, basis_, operator_spec_,
                                      operator_domains_, operator_data_,
                                      part_field_, part_tolerance_, part_smoothing);
      // do the remap
      swarm_remap.run(executor, true);

      // transfer data back to target mesh
      for (auto&& name : target_cellvar_names) {
        auto& swarm_field = target_swarm_state.get_field(name);
        double* mesh_field;
        target_state_.mesh_get_data(Wonton::CELL, name, &mesh_field);
        for (int i = 0; i < target_swarm_state.get_size(); i++)
          mesh_field[i] = swarm_field[i];
      }

      delete swarm_remap_ptr;
    }

    // Wonton::NODE VARIABLE SECTION ------------------------------------------

    // get node variable names
    std::vector<std::string> source_nodevar_names;
    std::vector<std::string> target_nodevar_names;
    for (int i = 0; i < nvars; ++i) {
      auto onwhat = source_state_.get_entity(source_vars_[i]);
      if (onwhat == Wonton::NODE) {
        source_nodevar_names.emplace_back(source_vars_[i]);
        target_nodevar_names.emplace_back(target_vars_[i]);
      }
    }

    if (not source_nodevar_names.empty()) {
      // convert mesh and state wrappers to swarm ones
      Wonton::Swarm<dim> source_swarm(source_mesh_, Wonton::NODE);
      Wonton::Swarm<dim> target_swarm(target_mesh_, Wonton::NODE);
      Wonton::SwarmState<dim> source_swarm_state(source_state_, Wonton::NODE);
      Wonton::SwarmState<dim> target_swarm_state(target_state_, Wonton::NODE);

      // create smoothing lengths
      Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths;
      std::vector<std::vector<double>> default_lengths(1, std::vector<double>(dim));

      if (geometry_ == Weight::FACETED) {
        throw std::runtime_error("Cannot do FACETED weights for nodal variables yet.");
      }

      int nnodes = (center_ == Scatter ? source_mesh_.num_owned_nodes()
                                       : target_mesh_.num_owned_nodes());

      smoothing_lengths.resize(nnodes, default_lengths);

      for (int i = 0; i < nnodes; i++) {
        double radius = 0.0;
        switch (center_) {
          case Scatter: Wonton::node_radius<dim>(source_mesh_, i, &radius); break;
          case Gather:  Wonton::node_radius<dim>(target_mesh_, i, &radius); break;
          default: break;
        }

        std::vector<std::vector<double>> h = smoothing_lengths[i];
        h[0] = std::vector<double>(dim, radius * smoothing_factor_);
        smoothing_lengths[i] = h;
      }

      // create swarm remap driver
      SwarmRemap swarm_remap(source_swarm, source_swarm_state,
                             target_swarm, target_swarm_state,
                             smoothing_lengths, kernel_,
                             geometry_, center_);

      swarm_remap.set_remap_var_names(source_nodevar_names, target_nodevar_names,
                                      LocalRegression, basis_);

      // do the remap
      swarm_remap.run(executor, true);

      // transfer data back to target mesh
      for (auto&& name : target_nodevar_names) {
        auto& swarm_field = target_swarm_state.get_field(name);
        double* mesh_field;
        target_state_.mesh_get_data(Wonton::NODE, name, &mesh_field);
        for (int i = 0; i < target_swarm_state.get_size(); i++)
          mesh_field[i] = swarm_field[i];
      }
    }

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
    float elapsed = timer::elapsed(tic);
    std::cout << "Mesh-Swarm-Mesh Time for Rank " << rank << " (s): " << elapsed << std::endl;
#endif
  }

private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::vector<std::string> source_vars_;
  std::vector<std::string> target_vars_;
  double smoothing_factor_ = 0.0;
  double boundary_factor_ = 0.0;
  Meshfree::Weight::Kernel kernel_ {};
  Meshfree::Weight::Geometry geometry_ {};
  Meshfree::WeightCenter center_ {};
  std::string part_field_;
  double part_tolerance_ = 0.0;
  Meshfree::EstimateType estimate_ {};
  Meshfree::basis::Type basis_ {};
  Meshfree::oper::Type operator_spec_ {};
  Wonton::vector<Meshfree::oper::Domain> operator_domains_ {};
  Wonton::vector<std::vector<Point<dim>>> operator_data_ {};
  int dim_ = 2;
};  // class MSM_Driver

}  // namespace Portage

#endif  // SRC_DRIVER_MESH_SWARM_MESH_H_
