/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

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

#define DEBUG_PART_BY_PART 0

/**
 * @file  parts.h
 *
 * @brief Manages source and target sub-meshes for part-by-part remap.
*/
namespace Portage {

  using Wonton::Entity_kind;
  using Wonton::Entity_type;
  using entity_weights_t = std::vector<Wonton::Weights_t>;

/**
 * @brief Manages source and target sub-meshes for part-by-part remap.
 *        It detects boundaries mismatch and provides the necessary fixup
 *        for partially filled and empty cells values.
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
    class TargetMesh_Wrapper = SourceMesh_Wrapper,
    class TargetState_Wrapper = SourceState_Wrapper
  >
class PartPair {
  // shortcut
  using entity_weights_t = std::vector<Weights_t>;

public:
  /**
   * @brief Construct a default source-target mesh parts pair.
   */
  PartPair() = default;

  /**
   * @brief Construct a source-target mesh parts pair.
   *
   * @param source_mesh     the source mesh
   * @param target_mesh     the target mesh
   * @param source_state    the source mesh data
   * @param target_state    the target mesh data
   * @param source_entities the list of source entities to remap
   * @param target_entities the list of target entities to remap
   * @param executor        the MPI executor to use
   */
  PartPair(
    SourceMesh_Wrapper    const& source_mesh,
    SourceState_Wrapper   const& source_state,
    TargetMesh_Wrapper    const& target_mesh,
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
    source_mesh_size_ =
      (onwhat == Entity_kind::CELL ? source_mesh_.num_owned_cells()
                                   : source_mesh_.num_owned_nodes());
    target_mesh_size_ =
      (onwhat == Entity_kind::CELL ? target_mesh_.num_owned_cells()
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
        source_mesh_, &source_entities_masks_, mycomm_
      );
    }
#endif

    source_entities_volumes_.resize(source_part_size_);
    target_entities_volumes_.resize(target_part_size_);
    intersection_volumes_.resize(target_part_size_);
    source_lookup_.reserve(source_part_size_);
    target_lookup_.reserve(target_part_size_);

    // set relative indexing and populate lookup hashtables
    for (int i = 0; i < source_part_size_; ++i) {
      auto const& s = source_entities[i];
      source_relative_index_[s] = i;
      source_lookup_.insert(s);
    }

    for (int i = 0; i < target_part_size_; ++i) {
      auto const& t = target_entities[i];
      target_relative_index_[t] = i;
      target_lookup_.insert(t);
    }
  }

  /**
   * @brief Delete a source-target mesh parts pair.
   */
  ~PartPair() = default;

  /**
   * @brief Do source and target meshes have a boundary mismatch?
   *
   * @return true if so, false otherwise.
   */
  bool has_mismatch() const { return has_mismatch_; }

  /**
   * @brief Is mismatch already tested?
   *
   * @return true if so, false otherwise.
   */
  bool is_mismatch_tested() const { return is_mismatch_tested_; }

  /**
   * @brief Get source part size.
   *
   * @return source entities list size
   */
  const int& source_part_size() { return source_part_size_; }

  /**
   * @brief Get target part size.
   *
   * @return target entities list size
   */
  const int& target_part_size() { return target_part_size_; }

  /**
   * @brief Retrieve the neighbors of the given entity on target mesh.
   *
   * @tparam entity_type the entity type [ALL|PARALLEL_OWNED]
   * @param  entity      the given entity
   * @return filtered    the filtered neighboring entities list.
   */
  template<Entity_type entity_type = Entity_type::ALL>
  std::vector<int> get_target_filtered_neighbors(int entity) {
    std::vector<int> neighbors, filtered;
    // retrieve neighbors
    if (onwhat == Entity_kind::CELL) {
      target_mesh_.cell_get_node_adj_cells(entity, entity_type, &neighbors);
    } else {
      target_mesh_.node_get_cell_adj_nodes(entity, entity_type, &neighbors);
    }
    // filter then
    filtered.reserve(neighbors.size());
    for (auto&& neigh : neighbors) {
      if (target_lookup_.count(neigh)) {
        filtered.emplace_back(neigh);
      }
    }
    return std::move(filtered);
  }


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
  bool check_mismatch(Portage::vector<entity_weights_t> const& source_ents_and_weights) {

    assert(do_part_by_part_);

    // ------------------------------------------
    // COMPUTE VOLUMES ON SOURCE AND TARGET PARTS
    // ------------------------------------------
    // collect volumes of entities that are not masked out and sum them up
    for (auto&& s : source_entities_) {
      auto const& i = source_relative_index_[s];
      source_entities_volumes_[i] = (
        onwhat == Entity_kind::CELL
          ? source_entities_masks_[s] * source_mesh_.cell_volume(s)
          : source_entities_masks_[s] * source_mesh_.dual_cell_volume(s)
      );
      #if DEBUG_PART_BY_PART
        std::printf("- source_volume[%02d]: %.3f\n", s, source_entities_volumes_[i]);
        if (i == source_part_size_ - 1)
          std::printf("=======\n");
      #endif
    }

    for (auto&& t : target_entities_) {
      auto const& i = target_relative_index_[t];
      target_entities_volumes_[i] = (
        onwhat == Entity_kind::CELL ? target_mesh_.cell_volume(t)
                                    : target_mesh_.dual_cell_volume(t)
      );
      #if DEBUG_PART_BY_PART
        std::printf("- target_volume[%02d]: %.3f\n", t, target_entities_volumes_[i]);
        if (i == target_part_size_ - 1)
          std::printf("=======\n");
      #endif
    }

    for (auto&& t : target_entities_) {
      auto const& i = target_relative_index_[t];
      // accumulate weights
      entity_weights_t const& weights = source_ents_and_weights[t];
      intersection_volumes_[i] = 0.;
      for (auto&& sw : weights) {
        // matched source cell should be in the source part
        if (source_lookup_.count(sw.entityID))
          intersection_volumes_[i] += sw.weights[0];
        #if DEBUG_PART_BY_PART
          std::printf("\tweights[target:%d][source:%d]: %f\n",
                      t, sw.entityID, sw.weights[0]);
        #endif
      }
      #if DEBUG_PART_BY_PART
        std::printf("intersect_volume[%02d]: %.3f\n", t, intersection_volumes_[i]);
      #endif
    }

    double source_volume = std::accumulate(source_entities_volumes_.begin(),
                                           source_entities_volumes_.end(), 0.);
    double target_volume = std::accumulate(target_entities_volumes_.begin(),
                                           target_entities_volumes_.end(), 0.);
    double intersect_volume = std::accumulate(intersection_volumes_.begin(),
                                              intersection_volumes_.end(), 0.);

    global_source_volume_    = source_volume;
    global_target_volume_    = target_volume;
    global_intersect_volume_ = intersect_volume;

#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {
      MPI_Allreduce(&source_volume, &global_source_volume_, 1, MPI_DOUBLE, MPI_SUM, mycomm_);
      MPI_Allreduce(&target_volume, &global_target_volume_, 1, MPI_DOUBLE, MPI_SUM, mycomm_);
      MPI_Allreduce(&intersect_volume, &global_intersect_volume_, 1, MPI_DOUBLE, MPI_SUM, mycomm_);

      if (rank_ == 0) {
        std::printf("source volume: %.3f\n", global_source_volume_);
        std::printf("target volume: %.3f\n", global_target_volume_);
        std::printf("intersect volume: %.3f\n", global_intersect_volume_);
      }
    }
#else
    std::printf("source volume: %.3f\n", source_volume);
    std::printf("target volume: %.3f\n", target_volume);
    std::printf("intersect volume: %.3f\n", intersect_volume);
#endif

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

    double const relative_voldiff_target =
      std::abs(global_intersect_volume_ - global_target_volume_)
      / global_target_volume_;

    if (relative_voldiff_source > tolerance_) {
      has_mismatch_ = true;

      if (rank_ == 0) {
        std::fprintf(stderr, "\n** MESH MISMATCH - some source cells ");
        std::fprintf(stderr, "are not fully covered by the target mesh\n");
      }
    }

    if (relative_voldiff_target > tolerance_) {
      has_mismatch_ = true;

      if (rank_ == 0) {
        std::fprintf(stderr, "\n** MESH MISMATCH - some target cells ");
        std::fprintf(stderr, "are not fully covered by the source mesh\n");
      }
    }

    if (not has_mismatch_) {
      is_mismatch_tested_ = true;
      return false;
    }

    // Discrepancy between intersection volume and source mesh volume PLUS
    // Discrepancy between intersection volume and target mesh volume
    relative_voldiff_ = relative_voldiff_source + relative_voldiff_target;

    // Collect the empty target cells in layers starting from the ones
    // next to partially or fully covered cells. At the end of this
    // section, partially or fully covered cells will all have a layer
    // number of 0 and empty cell layers will have positive layer
    // numbers starting from 1
    std::vector<int> empty_entities;
    empty_entities.reserve(target_part_size_);

    is_cell_empty_.resize(target_part_size_, false);

    for (auto&& entity : target_entities_) {
      auto const& i = target_relative_index_[entity];
      if (std::abs(intersection_volumes_[i]) < epsilon_) {
        empty_entities.emplace_back(entity);
        is_cell_empty_[i] = true;
      }
    }

    int nb_empty = empty_entities.size();
    int global_nb_empty = nb_empty;

#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {
      global_nb_empty = 0;
      MPI_Reduce(&nb_empty, &global_nb_empty, 1, MPI_INT, MPI_SUM, 0, mycomm_);
    }
#endif

    if (global_nb_empty > 0 and rank_ == 0) {
      if (onwhat == Entity_kind::CELL) {
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

        for (auto&& entity : empty_entities) {
          auto const& i = target_relative_index_[entity];
          // skip already set entities
          if (layer_num_[i] == 0) {

            auto neighbors = get_target_filtered_neighbors<Entity_type::ALL>(entity);

            for (auto&& neigh : neighbors) {
              auto const& j = target_relative_index_[neigh];
              // At least one neighbor has some material or will
              // receive some material (indicated by having a +ve
              // layer number)
              if (not is_cell_empty_[j] or layer_num_[j] != 0) {
                current_layer_entities.push_back(entity);
                break;
              }
            }
          }
        }

        // Tag the current layer cells with the next layer number
        for (auto&& entity : current_layer_entities) {
          auto const& t = target_relative_index_[entity];
          layer_num_[t] = nb_layers + 1;
        }
        nb_tagged += current_layer_entities.size();

        empty_layers_.emplace_back(current_layer_entities);
        nb_layers++;
      }
    }

    is_mismatch_tested_ = true;
    return has_mismatch_;
  }


  /**
   * @brief Repair the remapped field to account for boundary mismatch
   *
   * @param src_var_name        field variable on source mesh
   * @param trg_var_name        field variable on target mesh
   * @param global_lower_bound  lower limit on variable
   * @param global_upper_bound  upper limit on variable
   * @param conservation_tol    conservation tolerance treshold
   * @param maxiter             max number of iterations
   * @param partial_fixup_type  type of fixup in case of partial mismatch
   * @param empty_fixup_type    type of fixup in empty target entities
   *
   * ---------------------------------------------------------------------------
   * 'partial_fixup_type' can be one of three types:
   *
   * CONSTANT             - Fields will see no perturbations BUT REMAP WILL BE
   *                        NON-CONSERVATIVE (constant preserving, not linearity
   *                        preserving)
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
   * ---------------------------------------------------------------------------
   * 'empty_fixup_type' can be one of two types:
   *
   * LEAVE_EMPTY - Leave empty cells as is
   * EXTRAPOLATE - Fill empty cells with extrapolated values
   * FILL        - Fill empty cells with specified values (not yet implemented)
   * ---------------------------------------------------------------------------
  */
  bool fix_mismatch(std::string src_var_name,
                    std::string trg_var_name,
                    double global_lower_bound = -infinity_,
                    double global_upper_bound = infinity_,
                    double conservation_tol = tolerance_,
                    int maxiter = 5,
                    Partial_fixup_type partial_fixup_type = SHIFTED_CONSERVATIVE,
                    Empty_fixup_type empty_fixup_type = EXTRAPOLATE) {

    if (source_state_.field_type(onwhat, src_var_name) == Field_type::MESH_FIELD) {
      return fix_mismatch_meshvar(src_var_name, trg_var_name,
                                  global_lower_bound, global_upper_bound,
                                  conservation_tol, maxiter,
                                  partial_fixup_type, empty_fixup_type);
    }
  }

  /**
   * @brief Repair a remapped mesh field to account for boundary mismatch.
   *
   * @param src_var_name        field variable on source mesh
   * @param trg_var_name        field variable on source mesh
   * @param global_lower_bound  lower limit on variable value
   * @param global_upper_bound  upper limit on variable value
   * @param conservation_tol    conservation tolerance treshold
   * @param maxiter             max number of iterations
   * @param partial_fixup_type  type of fixup in case of partial mismatch
   * @param empty_fixup_type    type of fixup in empty target entities
   * @return true if correctly fixed, false otherwise.
   */
  bool fix_mismatch_meshvar(std::string const & src_var_name,
                            std::string const & trg_var_name,
                            double global_lower_bound,
                            double global_upper_bound,
                            double conservation_tol,
                            int maxiter,
                            Partial_fixup_type partial_fixup_type,
                            Empty_fixup_type empty_fixup_type) {

    // valid only for part-by-part scenario
    assert(do_part_by_part_);
    static bool hit_lower_bound  = false;
    static bool hit_higher_bound = false;

    // Now process remap variables
    // WARNING: absolute indexing
    double const* source_data;
    double*       target_data;

    source_state_.mesh_get_data(onwhat, src_var_name, &source_data);
    target_state_.mesh_get_data(onwhat, trg_var_name, &target_data);

    if (partial_fixup_type == LOCALLY_CONSERVATIVE) {
      // In interpolate step, we divided the accumulated integral (U)
      // in a target cell by the intersection volume (v_i) instead of
      // the target cell volume (v_c) to give a target field of u_t =
      // U/v_i. In partially filled cells, this will preserve a
      // constant source field but fill the cell with too much material
      // (this is the equivalent of requesting Partial_fixup_type::CONSTANT).
      // To restore conservation (as requested by
      // Partial_fixup_type::LOCALLY_CONSERVATIVE), we undo the division by
      // the intersection volume and then divide by the cell volume
      // (u'_t = U/v_c = u_t*v_i/v_c). This does not affect the values
      // in fully filled cells

      for (auto&& entity : target_entities_) {
        auto const& t = target_relative_index_[entity];
        if (not is_cell_empty_[t]) {

          #if DEBUG_PART_BY_PART
            std::printf("fixing target_data[%d] with locally conservative fixup\n", entity);
            std::printf("= before: %.3f", target_data[entity]);
          #endif

          auto const relative_voldiff =
            std::abs(intersection_volumes_[t] - target_entities_volumes_[t])
            / target_entities_volumes_[t];

          if (relative_voldiff > tolerance_) {
            target_data[entity] *= intersection_volumes_[t] / target_entities_volumes_[t];
          }
          #if DEBUG_PART_BY_PART
            std::printf(", after: %.3f\n", target_data[entity]);
          #endif
        }
      }
    }


    if (empty_fixup_type != LEAVE_EMPTY) {
      // Do something here to populate fully uncovered target
      // cells. We have layers of empty cells starting out from fully
      // or partially populated cells. We will assign every empty cell
      // in a layer the average value of all its populated neighbors.
      // IN A DISTRIBUTED MESH IT _IS_ POSSIBLE THAT AN EMPTY ENTITY
      // WILL NOT HAVE ANY OWNED NEIGHBOR IN THIS PARTITION THAT HAS
      // MEANINGFUL DATA TO EXTRAPOLATE FROM (we remap data only to
      // owned entities)
      int current_layer_number = 1;

      for (auto const& current_layer : empty_layers_) {
        for (auto&& entity : current_layer) {

          double averaged_value = 0.;
          int nb_extrapol = 0;
          auto neighbors =
            get_target_filtered_neighbors<Entity_type::PARALLEL_OWNED>(entity);

          for (auto&& neigh : neighbors) {
            auto const& i = target_relative_index_[neigh];
            if (layer_num_[i] < current_layer_number) {
              averaged_value += target_data[neigh];
              nb_extrapol++;
            }
          }
          if (nb_extrapol > 0) {
            averaged_value /= nb_extrapol;
          }
          #if DEBUG_PART_BY_PART
            else {
              std::fprintf(stderr,
                "No owned neighbors of empty entity to extrapolate data from\n"
              );
            }
          #endif
          target_data[entity] = averaged_value;
        }
        current_layer_number++;
      }
    }

    // if the fixup scheme is constant or locally conservative then we're done
    if (partial_fixup_type == CONSTANT or partial_fixup_type == LOCALLY_CONSERVATIVE) {
      return true;
    } else if (partial_fixup_type == SHIFTED_CONSERVATIVE) {

      // At this point assume that all cells have some value in them
      // for the variable
      // Now compute the net discrepancy between integrals over source
      // and target. Excess comes from target cells not fully covered by
      // source cells and deficit from source cells not fully covered by
      // target cells
      double source_sum = 0.;
      double target_sum = 0.;

      for (auto&& s : source_entities_) {
        auto const& i = source_relative_index_[s];
        source_sum += source_data[s] * source_entities_volumes_[i];
      }

      for (auto&& t : target_entities_) {
        auto const& i = target_relative_index_[t];
        target_sum += target_data[t] * target_entities_volumes_[i];
      }

      double global_source_sum = source_sum;
      double global_target_sum = target_sum;

#ifdef PORTAGE_ENABLE_MPI
      if (distributed_) {
        MPI_Allreduce(
          &source_sum, &global_source_sum, 1, MPI_DOUBLE, MPI_SUM, mycomm_
        );

        MPI_Allreduce(
          &target_sum, &global_target_sum, 1, MPI_DOUBLE, MPI_SUM, mycomm_
        );
      }
#endif

      double absolute_diff = global_target_sum - global_source_sum;
      double relative_diff = absolute_diff / global_source_sum;

      if (std::abs(relative_diff) < conservation_tol) {
        return true;  // discrepancy is too small - nothing to do
      }

      // Now redistribute the discrepancy among cells in proportion to
      // their volume. This will restore conservation and if the
      // original distribution was a constant it will make the field a
      // slightly different constant. Do multiple iterations because
      // we may have some leftover quantity that we were not able to
      // add/subtract from some cells in any given iteration. We don't
      // expect that process to take more than two iterations at the
      // most.

      double adj_target_volume;
      double global_adj_target_volume;
      double global_covered_target_volume;

      if (empty_fixup_type == Empty_fixup_type::LEAVE_EMPTY) {
        double covered_target_volume = 0.;
        for (auto&& entity : target_entities_) {
          auto const& t = target_relative_index_[entity];
          if (not is_cell_empty_[t]) {
            covered_target_volume += target_entities_volumes_[t];
          }
        }
        global_covered_target_volume = covered_target_volume;
#ifdef PORTAGE_ENABLE_MPI
        if (distributed_) {
          MPI_Allreduce(
            &covered_target_volume, &global_covered_target_volume,
            1, MPI_DOUBLE, MPI_SUM, mycomm_
          );
        }
#endif
        adj_target_volume = covered_target_volume;
        global_adj_target_volume = global_covered_target_volume;
      }
      else {
        adj_target_volume = std::accumulate(
          target_entities_volumes_.begin(), target_entities_volumes_.end(), 0
        );
        global_adj_target_volume = global_target_volume_;
      }

      // get the right entity type
      auto target_entity_type = [&](int entity) -> Entity_type {
        return (onwhat == Entity_kind::CELL ? target_mesh_.cell_get_type(entity)
                                            : target_mesh_.node_get_type(entity));
      };

      // sort of a "unit" discrepancy or difference per unit volume
      double udiff = absolute_diff / global_adj_target_volume;

      int iter = 0;
      while (std::abs(relative_diff) > conservation_tol and iter < maxiter) {

        for (auto&& entity : target_entities_) {
          auto const& t = target_relative_index_[entity];
          bool is_owned = target_entity_type(entity) == Entity_type::PARALLEL_OWNED;
          bool should_fix = (empty_fixup_type != LEAVE_EMPTY or not is_cell_empty_[t]);

          if (is_owned and should_fix) {
            if ((target_data[entity] - udiff) < global_lower_bound) {
              // Subtracting the full excess will make this cell violate the
              // lower bound. So subtract only as much as will put this cell
              // exactly at the lower bound
              target_data[entity] = global_lower_bound;

              if (not hit_lower_bound) {
                std::fprintf(stderr,
                  "Hit lower bound for cell %d (and maybe other ones) on rank %d\n",
                  t, rank_
                );
                hit_lower_bound = true;
              }
              // this cell is no longer in play for adjustment - so remove its
              // volume from the adjusted target_volume
              adj_target_volume -= target_entities_volumes_[t];

            } else if ((target_data[entity] - udiff) > global_upper_bound) {  // udiff < 0
              // Adding the full deficit will make this cell violate the
              // upper bound. So add only as much as will put this cell
              // exactly at the upper bound
              target_data[entity] = global_upper_bound;

              if (not hit_higher_bound) {
                std::fprintf(stderr,
                  "Hit upper bound for cell %d (and maybe other ones) on rank %d\n",
                  t, rank_
                );
                hit_higher_bound = true;
              }

              // this cell is no longer in play for adjustment - so remove its
              // volume from the adjusted target_volume
              adj_target_volume -= target_entities_volumes_[t];
            } else {
              // This is the equivalent of
              //           [curval*cellvol - diff*cellvol/meshvol]
              // curval = ---------------------------------------
              //                       cellvol
              target_data[entity] -= udiff;
            }
          }  // only non-empty cells
        }  // iterate through mesh cells

        // Compute the new integral over all processors
        target_sum = 0.;
        for (auto&& entity : target_entities_) {
          auto const& t = target_relative_index_[entity];
          target_sum += target_entities_volumes_[t] * target_data[entity];
        }

        global_target_sum = target_sum;
#ifdef PORTAGE_ENABLE_MPI
        if (distributed_) {
          MPI_Allreduce(
            &target_sum, &global_target_sum, 1, MPI_DOUBLE, MPI_SUM, mycomm_
          );
        }
#endif

        // If we did not hit lower or upper bounds, this should be
        // zero after the first iteration.  If we did hit some bounds,
        // then recalculate the discrepancy and discrepancy per unit
        // volume, but only taking into account volume of cells that
        // are not already at the bounds - if we use the entire mesh
        // volume for the recalculation, the convergence slows down
        absolute_diff = global_target_sum - global_source_sum;
        global_adj_target_volume = adj_target_volume;

#ifdef PORTAGE_ENABLE_MPI
        if (distributed_) {
          MPI_Allreduce(
            &adj_target_volume, &global_adj_target_volume,
            1, MPI_DOUBLE, MPI_SUM, mycomm_
          );
        }
#endif

        udiff = absolute_diff / global_adj_target_volume;
        relative_diff = absolute_diff / global_source_sum;

        // Now reset adjusted target mesh volume to be the full volume
        // in preparation for the next iteration
        global_adj_target_volume = (
          empty_fixup_type == LEAVE_EMPTY ? global_covered_target_volume
                                          : global_target_volume_
        );

        iter++;
      }  // while leftover is not zero

      if (std::abs(relative_diff) > conservation_tol) {
        if (rank_ == 0) {
          std::fprintf(stderr,
            "Redistribution not entirely successfully for variable %s\n"
            "Relative conservation error is %f\n"
            "Absolute conservation error is %f\n",
            src_var_name.data(), relative_diff, absolute_diff
          );
          return false;
        }
      }

      return true;
    } else {
      std::fprintf(stderr, "Unknown Partial fixup type\n");
      return false;
    }
  }

public:
  // references to user-provided entities lists
  std::vector<int> const& source_entities_;
  std::vector<int> const& target_entities_;
  // hashtables to have constant-time parts lookup queries in average case.
  // remark: for lookup purposes only, not meant to be iterated.
  std::unordered_set<int> source_lookup_;
  std::unordered_set<int> target_lookup_;

private:

  bool do_part_by_part_    = false;
  bool is_mismatch_tested_ = false;
  bool has_mismatch_       = false;

  // useful constants
  static constexpr double infinity_  = std::numeric_limits<double>::max();
  static constexpr double epsilon_   = std::numeric_limits<double>::epsilon();
  static constexpr double tolerance_ = 1.E2 * epsilon_;

  // source/target meshes and related states
  SourceMesh_Wrapper  const& source_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetMesh_Wrapper  const& target_mesh_;
  TargetState_Wrapper&       target_state_;

  // entities list and data for part-by-part
  // WARNING: use relative indexing
  int source_mesh_size_ = 0;
  int target_mesh_size_ = 0;
  int source_part_size_ = 0;
  int target_part_size_ = 0;

  std::vector<int>    source_entities_masks_   = {};
  std::vector<double> source_entities_volumes_ = {};
  std::vector<double> target_entities_volumes_ = {};
  std::vector<double> intersection_volumes_    = {};

  // data needed for mismatch checks
  double global_source_volume_    = 0.;
  double global_target_volume_    = 0.;
  double global_intersect_volume_ = 0.;
  double relative_voldiff_        = 0.;

  // empty target cells management
  std::map<int,int> source_relative_index_    = {};
  std::map<int,int> target_relative_index_    = {};
  std::vector<int>  layer_num_                = {};
  std::vector<bool> is_cell_empty_            = {};
  std::vector<std::vector<int>> empty_layers_ = {};

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
