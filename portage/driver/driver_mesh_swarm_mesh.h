/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_DRIVER_MESH_SWARM_MESH_H_
#define SRC_DRIVER_MESH_SWARM_MESH_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <string>
#include <limits>

#include "portage/support/portage.h"

#include "portage/support/basis.h"
#include "portage/support/weight.h"
#include "portage/support/operator.h"
#include "portage/support/faceted_setup.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
//#include "portage/search/search_simple_points.h"
//#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/driver/driver_swarm.h"

#ifdef PORTAGE_ENABLE_MPI
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
          template <size_t, class, class> class Accumulate,
          template<size_t, class> class Estimate,
          int Dim,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper = SourceMesh_Wrapper,
          class TargetState_Wrapper = SourceState_Wrapper>
class MSM_Driver {

 public:
  /*!
    @brief Constructor for running the interpolation driver.
    @param[in] sourceMesh A @c SourceMesh_Wrapper to the source mesh.
    @param[in] sourceState A @c SourceState_Wrapperfor the data that lives on the
    source mesh.
    @param[in] targetMesh A @c TargetMesh_Wrapper to the target mesh.
    @param[in,out] targetState A @c TargetState_Wrapper for the data that will
    be mapped to the target mesh.
    @param[in] smoothing_factor multiplies cell sizes to get smoothing lengths
    @param[in] boundary_factor with faceted weights only, multiplies center-to-face distance on boundary to get smoothing length
    @param[in] geometry weight function geometric configuration
    @param[in] kernel weight function kernel
    @param[in] center weight function centering (scatter for sources, gather for targets)
    @param[in] part_field name of field to use for part assignments, with faceted weights only
    @param[in] part_tolerance tolerance for determining if part assignment matches a neighbor's assignment
  */
  MSM_Driver
        (SourceMesh_Wrapper const& sourceMesh,
         SourceState_Wrapper const& sourceState,
         TargetMesh_Wrapper const& targetMesh,
	 TargetState_Wrapper& targetState, 
	 double smoothing_factor             = 1.25,
	 double boundary_factor              = 0.5,
	 Meshfree::Weight::Geometry geometry = Meshfree::Weight::TENSOR,
         Meshfree::Weight::Kernel   kernel   = Meshfree::Weight::B4,
         Meshfree::WeightCenter     center   = Meshfree::Gather, 
         std::string part_field = "NONE", 
         double part_tolerance = std::numeric_limits<double>::infinity())
      : source_mesh_(sourceMesh), source_state_(sourceState),
        target_mesh_(targetMesh), target_state_(targetState),
        smoothing_factor_(smoothing_factor),
        boundary_factor_(boundary_factor),
        geometry_(geometry),
        kernel_(kernel),
        center_(center),
        dim_(sourceMesh.space_dimension()),
        part_field_(part_field), 
        part_tolerance_(part_tolerance)
  {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
    assert(sourceMesh.space_dimension() == Dim);
    if (geometry == Meshfree::Weight::FACETED) {
      assert(kernel == Meshfree::Weight::POLYRAMP);
    }
  }

  /// Copy constructor (disabled)
  MSM_Driver(const MSM_Driver &) = delete;

  /// Assignment operator (disabled)
  MSM_Driver & operator = (const MSM_Driver &) = delete;

  /// Destructor
  ~MSM_Driver() {}

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
    @param[in] estimator_type What type of estimator to apply in particle stage
    @param[in] basis_type what regression basis to use
    @param[in] operator_spec what type of operator to use in regression, if any
    @param[in] operator_domains if using an integral operator, what domains of integration to use
    @param[in] operator_data node data for integral domains, if needed
  */

  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      Meshfree::EstimateType const& estimator_type = Meshfree::LocalRegression,
      Meshfree::Basis::Type const& basis_type = Meshfree::Basis::Unitary,
      Meshfree::Operator::Type operator_spec = Meshfree::Operator::LastOperator,
      Portage::vector<Meshfree::Operator::Domain> operator_domains = 
        vector<Meshfree::Operator::Domain>(0),
      Portage::vector<std::vector<Point<Dim>>> const& operator_data=
        vector<std::vector<Point<Dim>>>(0,std::vector<Point<Dim>>(0))) {
    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    for (int i = 0; i < nvars; ++i)
      assert(source_state_.get_entity(source_remap_var_names[i]) ==
             target_state_.get_entity(target_remap_var_names[i]));

    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    estimate_ = estimator_type;
    basis_ = basis_type;
    operator_spec_ = operator_spec;
    operator_domains_ = operator_domains;
    operator_data_ = operator_data;
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
    @brief Execute the remapping process
  */
  void run(Wonton::Executor_type const *executor = nullptr) {

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
      if (nprocs > 1) {
        distributed = true;
        std::cerr << "Cannot run Mesh-Swarm-Mesh driver in distributed mode yet\n";
        return;
      }
    }
#endif
    
    if (comm_rank == 0) std::printf("in MSM_Driver::run()...\n");

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;

    int nvars = source_remap_var_names_.size();

    float tot_seconds = 0.0, tot_seconds_srch = 0.0,
        tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
    struct timeval begin_timeval, end_timeval, diff_timeval;

    gettimeofday(&begin_timeval, 0);

    // CELL VARIABLE SECTION ---------------------------------------------------

    // get cell variable names
    std::vector<std::string> source_cellvar_names;
    std::vector<std::string> target_cellvar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == Entity_kind::CELL) {
        source_cellvar_names.emplace_back(source_remap_var_names_[i]);
        target_cellvar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    // Collect all cell based variables and remap them
    if (source_cellvar_names.size() > 0)
    {
      // convert mesh wrappers to swarms
      std::shared_ptr<Meshfree::Swarm<Dim>> source_swarm_ptr = 
        Meshfree::SwarmFactory<Dim>(source_mesh_, Entity_kind::CELL);
      std::shared_ptr<Meshfree::Swarm<Dim>> target_swarm_ptr = 
        Meshfree::SwarmFactory<Dim>(target_mesh_, Entity_kind::CELL);
      Meshfree::Swarm<Dim> &source_swarm(*source_swarm_ptr);
      Meshfree::Swarm<Dim> &target_swarm(*target_swarm_ptr);

      // convert state wrappers to swarm variety
      std::shared_ptr<Meshfree::SwarmState<Dim>> source_swarm_state_ptr = 
        Meshfree::SwarmStateFactory<Dim>(source_state_, Entity_kind::CELL);
      std::shared_ptr<Meshfree::SwarmState<Dim>> target_swarm_state_ptr = 
        Meshfree::SwarmStateFactory<Dim>(target_state_, Entity_kind::CELL);
      Meshfree::SwarmState<Dim> &source_swarm_state(*source_swarm_state_ptr);
      Meshfree::SwarmState<Dim> &target_swarm_state(*target_swarm_state_ptr);

      // set up smoothing lengths and extents
      Portage::vector<std::vector<std::vector<double>>> smoothing_lengths;
      Portage::vector<Wonton::Point<Dim>> weight_extents, other_extents;
      Portage::vector<std::vector<std::vector<double>>> part_smoothing; // only for faceted,scatter,parts
      if (geometry_ == Meshfree::Weight::FACETED) {
        if (part_field_ == "NONE") {
          if (center_ == Meshfree::Scatter) {
            Meshfree::Weight::faceted_setup_cell<Dim,SourceMesh_Wrapper>
              (source_mesh_, smoothing_lengths, weight_extents, smoothing_factor_, boundary_factor_);
          } else if (center_ == Meshfree::Gather) {
            Meshfree::Weight::faceted_setup_cell<Dim,TargetMesh_Wrapper>
              (target_mesh_, smoothing_lengths, weight_extents, smoothing_factor_, boundary_factor_);
          }
        } else {
          if (center_ == Meshfree::Scatter) {
            // Set smoothing_factor to 1/4 to make weight support exactly equal to cell volume.
            // Store these smoothing lengths off on the side for later use in the swarm driver.
            Portage::vector<Wonton::Point<Dim>> dummy_extents;
            Meshfree::Weight::faceted_setup_cell<Dim,SourceMesh_Wrapper>
              (source_mesh_, part_smoothing, dummy_extents, 0.25, 0.25);

            // Get usual smoothing lengths and extents 
            Meshfree::Weight::faceted_setup_cell<Dim,SourceMesh_Wrapper>
              (source_mesh_, smoothing_lengths, weight_extents, smoothing_factor_, boundary_factor_);
          } else if (center_ == Meshfree::Gather) {
            Meshfree::Weight::faceted_setup_cell<Dim,TargetMesh_Wrapper,TargetState_Wrapper>
              (target_mesh_, target_state_, part_field_, part_tolerance_,
               smoothing_lengths, weight_extents, smoothing_factor_, boundary_factor_);
          }
        }
      } else {
        int ncells;
        if      (center_ == Meshfree::Scatter) ncells = source_mesh_.num_owned_cells();
        else if (center_ == Meshfree::Gather)  ncells = target_mesh_.num_owned_cells();
        smoothing_lengths = vector<std::vector<std::vector<double>>> 
          (ncells, std::vector<std::vector<double>>(1, std::vector<double>(Dim)));
        for (int i=0; i<ncells; i++) {
          double radius;
          if      (center_ == Meshfree::Scatter) 
            Wonton::cell_radius<Dim>(source_mesh_, i, &radius);
          else if (center_ == Meshfree::Gather)  
            Wonton::cell_radius<Dim>(target_mesh_, i, &radius);
          std::vector<std::vector<double>> h=smoothing_lengths[i];
          h[0] = std::vector<double>(Dim, 2.*radius*smoothing_factor_);
          smoothing_lengths[i]=h;
        }
      }

      // create swarm remap driver
      using SwarmDriverType=
        Meshfree::SwarmDriver<
        Search,
	Accumulate,
	Estimate,
        Dim,
        Meshfree::Swarm<Dim>,
        Meshfree::SwarmState<Dim>,
        Meshfree::Swarm<Dim>,
        Meshfree::SwarmState<Dim>
        >;

      SwarmDriverType *swarm_driver_ptr;

      if (geometry_ == Meshfree::Weight::FACETED) {
        if (center_ == Meshfree::Scatter) {
          swarm_driver_ptr = new SwarmDriverType
            ( source_swarm, source_swarm_state,
              target_swarm, target_swarm_state,
              smoothing_lengths, 
              weight_extents, other_extents, 
              center_);
        } else if  (center_ == Meshfree::Gather) {
          swarm_driver_ptr = new SwarmDriverType
            ( source_swarm, source_swarm_state,
              target_swarm, target_swarm_state,
              smoothing_lengths, 
              other_extents, weight_extents,  
              center_);
        }
      } else {
        swarm_driver_ptr = new SwarmDriverType
          ( source_swarm, source_swarm_state,
            target_swarm, target_swarm_state,
            smoothing_lengths,
            kernel_,
            geometry_,
            center_);
      }
      SwarmDriverType &swarm_driver(*swarm_driver_ptr);

      swarm_driver.set_remap_var_names(source_cellvar_names,
                                       target_cellvar_names,
                                       estimate_,
                                       basis_,
				       operator_spec_,
				       operator_domains_,
                                       operator_data_,
                                       part_field_,
                                       part_tolerance_, 
                                       part_smoothing);

      // do the remap
      swarm_driver.run(executor, true);

      // transfer data back to target mesh
      for (auto name=target_cellvar_names.begin();
                name!=target_cellvar_names.end(); name++)
      {
        typename Meshfree::SwarmState<Dim>::DblVecPtr sfield;
        target_swarm_state.get_field(*name, sfield);
        double *mfield;
        target_state_.mesh_get_data(Entity_kind::CELL, *name, &mfield);
        for (int i=0; i<target_swarm_state.get_size(); i++) {
          mfield[i] = (*sfield)[i];
        }
      }

      delete swarm_driver_ptr;
    }

    // NODE VARIABLE SECTION ---------------------------------------------------

    // get node variable names
    std::vector<std::string> source_nodevar_names;
    std::vector<std::string> target_nodevar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == Entity_kind::NODE) {
        source_nodevar_names.emplace_back(source_remap_var_names_[i]);
        target_nodevar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_nodevar_names.size() > 0) {
      // convert mesh wrappers to swarm variety
      std::shared_ptr<Meshfree::Swarm<Dim>> source_swarm_ptr = 
        Meshfree::SwarmFactory<Dim>(source_mesh_, Entity_kind::NODE);
      std::shared_ptr<Meshfree::Swarm<Dim>> target_swarm_ptr = 
        Meshfree::SwarmFactory<Dim>(target_mesh_, Entity_kind::NODE);
      Meshfree::Swarm<Dim> &source_swarm(*source_swarm_ptr);
      Meshfree::Swarm<Dim> &target_swarm(*target_swarm_ptr);

      // convert state wrappers to swarm variety
      std::shared_ptr<Meshfree::SwarmState<Dim>> source_swarm_state_ptr = 
        Meshfree::SwarmStateFactory<Dim>(source_state_, Entity_kind::NODE);
      std::shared_ptr<Meshfree::SwarmState<Dim>> target_swarm_state_ptr = 
        Meshfree::SwarmStateFactory<Dim>(target_state_, Entity_kind::NODE);
      Meshfree::SwarmState<Dim> &source_swarm_state(*source_swarm_state_ptr);
      Meshfree::SwarmState<Dim> &target_swarm_state(*target_swarm_state_ptr);

      // create smoothing lengths
      if (geometry_ == Meshfree::Weight::FACETED) {
        std::cerr << "Cannot do FACETED weights for nodal variables yet.\n";
        return;
      }
      int nnodes;
      if      (center_ == Meshfree::Scatter) nnodes = source_mesh_.num_owned_nodes();
      else if (center_ == Meshfree::Gather)  nnodes = target_mesh_.num_owned_nodes();
      vector<std::vector<std::vector<double>>> smoothing_lengths
        (nnodes, std::vector<std::vector<double>>(1, std::vector<double>(Dim)));
      for (int i=0; i<nnodes; i++) {
        double radius;
        if      (center_ == Meshfree::Scatter) Wonton::node_radius<Dim>(source_mesh_, i, &radius);
	else if (center_ == Meshfree::Gather)  Wonton::node_radius<Dim>(target_mesh_, i, &radius);
	std::vector<std::vector<double>> h=smoothing_lengths[i];
        h[0] = std::vector<double>(Dim, radius*smoothing_factor_);
	smoothing_lengths[i]=h;
      }

      // create swarm remap driver
      Meshfree::SwarmDriver<
        Search,
        Accumulate,
        Estimate,
        Dim,
        Meshfree::Swarm<Dim>,
        Meshfree::SwarmState<Dim>,
        Meshfree::Swarm<Dim>,
        Meshfree::SwarmState<Dim>
      > swarm_driver(source_swarm, source_swarm_state,
                     target_swarm, target_swarm_state,
                     smoothing_lengths,
                     kernel_,
                     geometry_,
                     center_);

      swarm_driver.set_remap_var_names(source_nodevar_names,
                                       target_nodevar_names,
                                       Meshfree::LocalRegression,
                                       basis_);

      // do the remap
      swarm_driver.run(executor, true);

      // transfer data back to target mesh
      for (auto name=target_nodevar_names.begin();
                name!=target_nodevar_names.end(); name++)
      {
        typename Meshfree::SwarmState<Dim>::DblVecPtr sfield;
        target_swarm_state.get_field(*name, sfield);
        double *mfield;
        target_state_.mesh_get_data(Entity_kind::NODE, *name, &mfield);
        for (int i=0; i<target_swarm_state.get_size(); i++) {
          mfield[i] = (*sfield)[i];
        }
      }
    }

    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    tot_seconds_interp = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

    std::cout << "Mesh-Swarm-Mesh Time for Rank " << comm_rank << " (s): " <<
        tot_seconds << std::endl;
  }


 private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  double smoothing_factor_;
  double boundary_factor_;
  Meshfree::Weight::Kernel kernel_;
  Meshfree::Weight::Geometry geometry_;
  Meshfree::WeightCenter center_;
  std::string part_field_;
  double part_tolerance_;
  Meshfree::EstimateType estimate_;
  Meshfree::Basis::Type basis_;
  Meshfree::Operator::Type operator_spec_;
  Portage::vector<Meshfree::Operator::Domain> operator_domains_;
  Portage::vector<std::vector<Point<Dim>>> operator_data_;
  unsigned int dim_;
};  // class MSM_Driver



}  // namespace Portage

#endif  // SRC_DRIVER_MESH_SWARM_MESH_H_
