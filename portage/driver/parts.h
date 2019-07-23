//
// Created by Hoby Rakotarivelo on 2019-07-22.
//

#ifndef PORTAGE_DRIVER_PARTS_H
#define PORTAGE_DRIVER_PARTS_H

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <unordered_set>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <limits>

#include "portage/support/portage.h"

/*!
  @file detect_mismatch.h
  @brief Detect if the boundaries of the source and target mesh are not exactly coincident
*/

namespace Portage {

  using Wonton::Entity_kind;
  using Wonton::Entity_type;
  using Wonton::Weights_t;
  using entity_weights_t = std::vector<Weights_t>;

/**
 * @brief Handle pairs of source-target entities
 * parts for mismatch fixup.
 *
 * @tparam D                   meshes dimension
 * @tparam onwhat              the entity kind for remap [cell|node]
 * @tparam SourceMesh_Wrapper  the source mesh wrapper to use
 * @tparam SourceState_Wrapper the source state wrapper to use
 * @tparam TargetMesh_Wrapper  the target mesh wrapper to use
 * @tparam TargetState_Wrapper the target state wrapper to use
 */
  template<int D, Entity_kind onwhat,
    class SourceMesh_Wrapper, class SourceState_Wrapper,
    class TargetMesh_Wrapper, class TargetState_Wrapper>
class Parts {
  // shortcut
  using entity_weights_t = std::vector<Weights_t>;

public:
  Parts() = default;
  Parts(
    SourceMesh_Wrapper    const& source_mesh,
    TargetMesh_Wrapper    const& target_mesh,
    SourceState_Wrapper   const& source_state,
    TargetState_Wrapper&         target_state,
    std::vector<int>      const& source_entities,
    std::vector<int>      const& target_entities,
    Wonton::Executor_type const* executor
  ) : source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      target_state_(target_state),
      source_entities_(source_entities),
      target_entities_(target_entities)
  {
#ifdef PORTAGE_ENABLE_MPI
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      distributed_ = true;
      mycomm_ = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm_, &rank_);
      MPI_Comm_size(mycomm_, &nprocs_);
    }
#endif
    // update source/target entities metadata
    source_mesh_size_ = (onwhat == CELL ? source_mesh_.num_owned_cells()
                                        : source_mesh_.num_owned_nodes());
    target_mesh_size_ = (onwhat == CELL ? target_mesh_.num_owned_cells()
                                        : target_mesh_.num_owned_nodes());

    source_part_size_ = source_entities.size();
    target_part_size_ = target_entities.size();
    do_part_by_part_  = (source_part_size_ > 0 and target_part_size_ > 0);

    // Get info about which entities on this processor should be
    // masked out and not accounted for in calculations because they
    // were already encountered on a lower rank processor. We have to
    // do this because in our distributed runs, our source partitions
    // don't form a strict tiling (no overlaps) after redistribution
    source_entities_masks_.resize(source_mesh_size_, 1);
#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {
      get_unique_entity_masks<onwhat, SourceMesh_Wrapper>(
        source_mesh_, &source_ent_masks, mycomm_
      );
    }
#endif

    source_entities_volumes_.reserve(source_part_size_);
    target_entities_volumes_.reserve(target_part_size_);
    intersection_volumes_.reserve(target_part_size);
  }
  ~Parts() = default;

  // has this problem been found to have mismatched mesh boundaries?
  bool has_mismatch() const { return mismatch_; }

  /**
 * @brief Check and fix source/target boundaries mismatch.
 *
 * T is the target mesh, S is the source mesh
 * T_i is the i'th cell of the target mesh and
 * |*| signifies the extent/volume of an entity
 *
 * - if sum_i(|T_i intersect S|) NE |S|, some source cells are not
 *   completely covered by the target mesh.
 * - if sum_i(|T_i intersect S|) NE |T|, some target cells are not
 *   completely covered by the source mesh
 * - if there is a mismatch, adjust values for fields to account so that
 *   integral quantities are conserved and a constant field is readjusted
 *   to a different constant
 *
 * WARNING: 'source_ents_and_weights' is a GLOBAL list defined on the entire
 *          target mesh.
 * @param source... source entities ID and weights for each target entity.
 * @return true if a mismatch has been identified, false otherwise.
 */
  bool test_mismatch(Portage::vector<entity_weights_t> const& source_ents_and_weights) {

    assert(do_part_by_part_);

    // ------------------------------------------
    // COMPUTE VOLUMES ON SOURCE AND TARGET PARTS
    // ------------------------------------------
    // collect volumes of entities that are not masked out and sum them up
    for (auto&& s : source_entities_) {
      source_entities_volumes.emplace_back(
        onwhat == CELL ? source_entities_masks_[s] * source_mesh_.cell_volume(s)
                       : source_entities_masks_[s] * source_mesh_.dual_cell_volume(s)
      );
    }

    for (auto&& t : target_entities_) {
      target_entities_volumes.emplace_back(
        onwhat == CELL ? target_mesh_.cell_volume(t)
                       : target_mesh_.dual_cell_volume(t)
      );
    }

    for (int i = 0; i < target_part_size_; ++i) {
      auto const& t = target_entities_[i];
      // accumulate weights
      entity_weights_t const& weights = source_ents_and_weights[t];
      intersection_volumes_[i] = 0.;
      for (auto&& sw : weights) {
        intersection_volumes_[i] += sw.weights[0];
      }
    }

    double source_volume = std::accumulate(source_entities_volumes_.begin(),
                                           source_entities_volumes_.end(), 0);
    double target_volume = std::accumulate(target_entities_volumes_.begin(),
                                           target_entities_volumes_.end(), 0);
    double intersect_volume = std::accumulate(intersection_volumes_.begin(),
                                              intersection_volumes_.end(), 0);

    global_source_volume_    = source_volume_;
    global_target_volume_    = target_volume_;
    global_intersect_volume_ = intersect_volume;

#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {
      MPI_Allreduce(
        &source_volume_, &global_source_volume_, 1, MPI_DOUBLE, MPI_SUM, mycomm_
      );

      MPI_Allreduce(
        &target_volume_, &global_target_volume_, 1, MPI_DOUBLE, MPI_SUM, mycomm_
      );

      MPI_Allreduce(
        &intersect_volume, &global_intersect_volume_, 1, MPI_DOUBLE, MPI_SUM, mycomm_
      );
    }
#endif

    auto const epsilon = std::numeric_limits<double>::epsilon();
    auto const tolerance_voldiff = 1.E2 * epsilon;

    // In our initial redistribution phase, we will move as many
    // source cells as needed from different partitions to cover the
    // target partition, UNLESS NO SOURCE CELLS COVERING THE TARGET
    // CELLS EXIST GLOBALLY. So, at this stage, we can just check
    // on-rank intersections only

    // Are some source cells not fully covered by target cells?
    // Are some target cells not fully covered by source cells

    double const relative_voldiff_source =
      std::abs(global_intersect_volume_ - global_source_volume_)
      / global_source_volume_;

    double const relative_voldiff_target_ =
      std::abs(global_intersect_volume_ - global_target_volume_)
      / global_target_volume_;

    if (relative_voldiff_source_ > tolerance_voldiff) {
      has_mismatch_ = true;

      if (rank_ == 0) {
        std::fprintf(stderr, "\n** MESH MISMATCH - some source cells ");
        std::fprintf(stderr, "are not fully covered by the target mesh\n");
      }
    }

    if (relative_voldiff_target_ > tolerance_voldiff) {
      has_mismatch_ = true;

      if (rank_ == 0) {
        std::fprintf(stderr, "\n** MESH MISMATCH - some target cells ");
        std::fprintf(stderr, "are not fully covered by the source mesh\n");
      }
    }

    if (not has_mismatch_)
      return false;

    // Discrepancy between intersection volume and source mesh volume PLUS
    // Discrepancy between intersection volume and target mesh volume
    relative_voldiff_ = relative_voldiff_source + relative_voldiff_target;

    // Collect the empty target cells in layers starting from the ones
    // next to partially or fully covered cells. At the end of this
    // section, partially or fully covered cells will all have a layer
    // number of 0 and empty cell layers will have positive layer
    // numbers starting from 1
    std::vector<int> empty_entities;

    is_cell_empty_.resize(target_part_size_, false);
    empty_entities.reserve(target_part_size_);

    for (int i = 0; i < target_part_size_; ++i) {
      auto const& t = target_entities_[i];
      if (std::abs(intersection_volumes_[i]) < epsilon) {
        empty_entities.emplace_back(t);
        is_cell_empty_[i] = true;
      }
    }

    int nb_empty = empty_entities.size();
    int global_nb_empty = nb_empty;

#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {
      int nb_empty_all[nprocs_];
      MPI_Gather(&nb_empty, 1, MPI_INT, nb_empty_all, 1, MPI_INT, 0, mycomm_);
      global_nb_empty = std::accumulate(nb_empty_all, nb_empty_all + nprocs_, 0.);
    }
#endif

    if (global_nb_empty > 0 and rank_ == 0) {
      if (onwhat == CELL) {
        std::fprintf(stderr,
          "One or more target cells are not covered by ANY source cells.\n"
          "Will assign values based on their neighborhood\n"
        );
      }
      else {
        std::fprintf(stderr,
          "One/more target dual cells are not covered by ANY source dual cells.\n"
          "Will assign values based on their neighborhood\n"
        );
      }
    }

    if (nb_empty > 0) {
      layer_num_.resize(target_part_size_, 0);

      int nb_layers = 0;
      int nb_tagged = 0;
      int old_nb_tagged = -1;

      while (nb_tagged < nb_empty and nb_tagged > old_nb_tagged) {
        old_nb_tagged = nb_tagged;

        std::vector<int> current_layer_entities;
        std::vector<int> neighbors;

        for (auto&& current_entity : empty_entities) {
          // skip already set entities
          if (layer_num_[current_entity] == 0) {
            neighbors.clear();
            if (onwhat == CELL) {
              target_mesh_.cell_get_node_adj_cells(current_entity, ALL, &neighbors);
            } else {
              target_mesh_.node_get_cell_adj_nodes(current_entity, ALL, &neighbors);
            }
            for (auto&& neigh : neighbors) {
              if (not is_cell_empty_[neigh] or layer_num_[neigh] != 0) {
                // At least one neighbor has some material or will
                // receive some material (indicated by having a +ve
                // layer number)
                current_layer_entities.push_back(current_entity);
                break;
              }
            }
          }
        }

        // Tag the current layer cells with the next layer number
        for (auto&& current_entity : current_layer_entities) {
          layer_num_[current_entity] = nb_layers + 1;
        }
        nb_tagged += current_layer_entities.size();

        empty_layers_.push_back(current_layer_entities);
        nb_layers++;
      }
    }

    return has_mismatch_;
  }


  /**
   * @brief Repair the remapped field to account for boundary mismatch
   * @param src_var_name        field variable on source mesh
   * @param trg_var_name        field variable on target mesh
   * @param global_lower_bound  lower limit on variable
   * @param global_upper_bound  upper limit on variable
   * @param partial_fixup_type  type of fixup in case of partial mismatch
   * @param empty_fixup_type    type of fixup in empty target entities
   *
   * partial_fixup_type can be one of three types:
   *
   * CONSTANT - Fields will see no perturbations BUT REMAP WILL BE
   *            NON-CONSERVATIVE (constant preserving, not linearity
   *            preserving)
   * LOCALLY_CONSERVATIVE - REMAP WILL BE LOCALLY CONSERVATIVE (target cells
   *                        will preserve the integral quantities received from
   *                        source mesh overlap) but perturbations will
   *                        occur in the field (constant fields may not stay
   *                        constant if there is mismatch)
   * SHIFTED_CONSERVATIVE - REMAP WILL BE CONSERVATIVE and field
   *                        perturbations will be minimum but field
   *                        values may be shifted (Constant fields
   *                        will be shifted to different constant; no
   *                        guarantees on linearity preservation)
   *
   * empty_fixup_type can be one of two types:
   *
   * LEAVE_EMPTY - Leave empty cells as is
   * EXTRAPOLATE - Fill empty cells with extrapolated values
   * FILL        - Fill empty cells with specified values (not yet implemented)
  */
  bool fix_mismatch(std::string src_var_name, std::string trg_var_name,
                    double global_lower_bound = -std::numeric_limits<double>::max(),
                    double global_upper_bound = std::numeric_limits<double>::max(),
                    double conservation_tol = 100*std::numeric_limits<double>::epsilon(),
                    int maxiter = 5,
                    Partial_fixup_type partial_fixup_type =
                      Partial_fixup_type::SHIFTED_CONSERVATIVE,
                    Empty_fixup_type empty_fixup_type =
                      Empty_fixup_type::EXTRAPOLATE) {

    if (source_state_.field_type(onwhat, src_var_name) == Field_type::MESH_FIELD) {
      return fix_mismatch_meshvar(src_var_name, trg_var_name,
                                  global_lower_bound, global_upper_bound,
                                  conservation_tol, maxiter,
                                  partial_fixup_type, empty_fixup_type);
    }
  }

  /**
   * @brief Repair a remapped mesh field to account for boundary mismatch
   */
  bool fix_mismatch_meshvar(std::string const & src_var_name,
                            std::string const & trg_var_name,
                            double global_lower_bound,
                            double global_upper_bound,
                            double conservation_tol = 1e2*std::numeric_limits<double>::epsilon(),
                            int maxiter = 5,
                            Partial_fixup_type partial_fixup_type =
                            Partial_fixup_type::SHIFTED_CONSERVATIVE,
                            Empty_fixup_type empty_fixup_type =
                            Empty_fixup_type::EXTRAPOLATE) {
    // TODO
    return false;
  }


private:

  bool do_part_by_part_ = false;
  bool has_mismatch_    = false;

  // source/target meshes and related states
  SourceMesh_Wrapper  const& source_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetMesh_Wrapper  const& target_mesh_;
  TargetState_Wrapper&       target_state_;

  // entities list and data for part-by-part
  int source_mesh_size_ = 0;
  int target_mesh_size_ = 0;
  int source_part_size_ = 0;
  int target_part_size_ = 0;

  std::vector<int> const& source_entities_;
  std::vector<int> const& target_entities_;
  std::vector<int>        source_entities_masks_   = {};
  std::vector<double>     source_entities_volumes_ = {};
  std::vector<double>     target_entities_volumes_ = {};
  std::vector<double>     intersection_volumes_    = {};

  // data needed for mismatch checks
  double global_source_volume_    = 0.;
  double global_target_volume_    = 0.;
  double global_intersect_volume_ = 0.;
  double relative_voldiff_        = 0.;

  // empty target cells management
  std::vector<std::vector<int>> empty_layers_  = {};
  std::vector<int>              layer_num_     = {};
  std::vector<bool>             is_cell_empty_ = {};

  // MPI
  int rank_         = 0;
  int nprocs_       = 1;
  bool distributed_ = false;
#ifdef PORTAGE_ENABLE_MPI
    MPI_Comm mycomm_ = MPI_COMM_NULL;
#endif
};

} // end namespace Portage
#endif //PORTAGE_PARTS_H
