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

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"
#include "portage/support/basis.h"
#include "portage/support/weight.h"
#include "portage/search/search_simple_points.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

#ifdef ENABLE_MPI
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
  @brief MSM_Driver provides the API to mapping from one mesh to another.
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
template <int Dim,
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
  */
  MSM_Driver(SourceMesh_Wrapper const& sourceMesh,
         SourceState_Wrapper const& sourceState,
         TargetMesh_Wrapper const& targetMesh,
         TargetState_Wrapper& targetState)
      : source_mesh_(sourceMesh), source_state_(sourceState),
        target_mesh_(targetMesh), target_state_(targetState),
        dim_(sourceMesh.space_dimension()) {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
    assert(sourceMesh.space_dimension() == Dim);
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
  */

  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names) {
    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    for (int i = 0; i < nvars; ++i)
      assert(source_state_.get_entity(source_remap_var_names[i]) ==
             target_state_.get_entity(target_remap_var_names[i]));

    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
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
  void run(bool distributed) {

#ifndef ENABLE_MPI
    if (distributed) {
      std::cout << "Request is for a parallel run but Portage is compiled for serial runs only\n";
      return;
    }
#endif


    int comm_rank = 0;

#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

    if (comm_rank == 0) std::printf("in MSM_Driver::run()...\n");

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;

    int nvars = source_remap_var_names_.size();

    // get cell variable names
    std::vector<std::string> source_cellvar_names;
    std::vector<std::string> target_cellvar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == CELL) {
        source_cellvar_names.emplace_back(source_remap_var_names_[i]);
        target_cellvar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    // convert incoming and outgoing mesh wrappers to flat variety
    Flat_Mesh_Wrapper<> source_mesh_flat;
    Flat_State_Wrapper<> source_state_flat;

    Flat_Mesh_Wrapper<double> source_mesh_flat;
    source_mesh_flat.initialize(source_mesh_);

    Flat_Mesh_Wrapper<double> target_mesh_flat;
    target_mesh_flat.initialize(target_mesh_);

    // convert flat wrappers to swarm variety
    Meshfree::Swarm<Dim> source_swarm(source_mesh_flat, CELL);
    Meshfree::Swarm<Dim> target_swarm(target_mesh_flat, CELL);

    // Collect all cell based variables and remap them
    if (source_cellvar_names.size() > 0)
    {
      // convert incoming and outgoing state wrappers to flat variety
      Flat_State_Wrapper<double> source_state_flat;
      source_state_flat.initialize(source_state_, source_cellvar_names);

      Flat_State_Wrapper<double> target_state_flat;
      target_state_flat.initialize(target_state_, target_cellvar_names);

      // convert flat wrappers to swarm variety
      Meshfree::SwarmState<Dim> source_swarm_state(source_mesh_flat, CELL,
                                                   source_state_flat);
      Meshfree::SwarmState<Dim> target_swarm_state(target_mesh_flat, CELL,
                                                   target_state_flat);

      // create smoothing lengths
      using std::vector;
      vector<vector<vector<double>>> smoothing_lengths
        (source_mesh_flat.num_owned_cells(),
         vector<vector<double>>(1, vector<double>(Dim)));
      for (int i=0; i<target_mesh_flat.num_owned_cells(); i++) {
        double radius;
        cell_radius<Dim>(source_mesh_flat, i, &radius);
        smoothing_lengths[i][0] = vector<double>(Dim, radius*1.5);
      }

      // create swarm remap driver
      Meshfree::SwarmDriver<
        SearchPointsByCells,
        Meshfree::Accumulate,
        Meshfree::Estimate,
        Dim,
        Meshfree::Swarm<Dim>,
        Meshfree::SwarmState<Dim>>
      > swarm_driver(source_swarm, source_swarm_state,
                     target_swarm, target_swarm_state,
                     smoothing_lengths,
                     Meshfree::Weight::B4,
                     Meshfree::Weight::ELLIPTIC,
                     Meshfree::Scatter);

      swarm_driver.set_remap_var_names(source_remap_var_names_,
                                       target_remap_var_names_,
                                       Meshfree::LocalRegression,
                                       Meshfree::Basis::Linear);

      // do the remap
      swarm_driver.run(false);

      // transfer data back to target mesh
      for (auto name=target_remap_var_names_.begin();
                name!=target_remap_var_names_.end(); name++)
      {
        typename Meshfree::SwarmState::DblVecPtr sfield;
        target_swarm_state.get_field(*name, sfield);
        double *mfield;
        target_state_flat.get_data(CELL, *name, &mfield);
        for (int i=0; i<target_swarm_state.get_size(); i++) {
          mfield[i] = (*sfield)[i];
        }
      }
    }

    // get node variable names
    std::vector<std::string> source_nodevar_names;
    std::vector<std::string> target_nodevar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat =
          source_state_.get_entity(source_remap_var_names_[i]);

      if (onwhat == NODE) {
        source_nodevar_names.emplace_back(source_remap_var_names_[i]);
        target_nodevar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_nodevar_names.size() > 0) {
      // convert incoming and outgoing state wrappers to flat variety
      Flat_State_Wrapper<double> source_state_flat;
      source_state_flat.initialize(source_state_, source_nodevar_names);

      Flat_State_Wrapper<double> target_state_flat;
      target_state_flat.initialize(target_state_, target_nodevar_names);

      // convert flat wrappers to swarm variety
      Meshfree::SwarmState<Dim> source_swarm_state(source_mesh_flat, NODE,
                                                   source_state_flat);
      Meshfree::SwarmState<Dim> target_swarm_state(target_mesh_flat, NODE,
                                                   target_state_flat);

      // create smoothing lengths
      using std::vector;
      vector<vector<vector<double>>> smoothing_lengths
        (source_mesh_flat.num_owned_nodes(),
         vector<vector<double>>(1, vector<double>(Dim)));
      for (int i=0; i<target_mesh_flat.num_owned_nodes(); i++) {
        double radius;
        node_radius<Dim>(source_mesh_flat, i, &radius);
        smoothing_lengths[i][0] = vector<double>(Dim, radius*1.5);
      }

      // create swarm remap driver
      Meshfree::SwarmDriver<
        SearchPointsByCells,
        Meshfree::Accumulate,
        Meshfree::Estimate,
        Dim,
        Meshfree::Swarm<Dim>,
        Meshfree::SwarmState<Dim>>
      > swarm_driver(source_swarm, source_swarm_state,
                     target_swarm, target_swarm_state,
                     smoothing_lengths,
                     Meshfree::Weight::B4,
                     Meshfree::Weight::ELLIPTIC,
                     Meshfree::Scatter);

      swarm_driver.set_remap_var_names(source_remap_var_names_,
                                       target_remap_var_names_,
                                       Meshfree::LocalRegression,
                                       Meshfree::Basis::Linear);

      // do the remap
      swarm_driver.run(false);

      // transfer data back to target mesh
      for (auto name=target_remap_var_names_.begin();
                name!=target_remap_var_names_.end(); name++)
      {
        typename Meshfree::SwarmState::DblVecPtr sfield;
        target_swarm_state.get_field(*name, sfield);
        double *mfield;
        target_state_flat.get_data(CELL, *name, &mfield);
        for (int i=0; i<target_swarm_state.get_size(); i++) {
          mfield[i] = (*sfield)[i];
        }
      }
    }

  }


 private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  unsigned int dim_;
};  // class MSM_Driver



}  // namespace Portage

#endif  // SRC_DRIVER_MESH_SWARM_MESH_H_
