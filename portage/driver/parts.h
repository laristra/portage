/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_DRIVER_PARTS_H
#define PORTAGE_DRIVER_PARTS_H

#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include <iterator>
#include <unordered_set>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <limits>

#include "portage/support/portage.h"
#include "portage/driver/fix_mismatch.h"

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

  template<Entity_kind onwhat, class Mesh, class State>
  class Part {
  public:
    /**
     * @brief Default constructor
     */
    Part() = default;

    /**
     * @brief Construct a mesh part object.
     *
     * @param mesh: the mesh wrapper to use
     * @param state: the state wrapper to use
     * @param entities: the list of entities to remap
     */
    Part(Mesh const& mesh, State& state, std::vector<int> const& entities)
      : mesh_(mesh),
        state_(state),
        entities_(entities)  // deep-copy
    {
      size_ = entities.size();
      if (size_ > 0) {
        volumes_.resize(size_);
        lookup_.reserve(size_);

        // set relative indexing and populate lookup hashtables
        for (int i = 0; i < size_; ++i) {
          auto const& s = entities[i];
          index_[s] = i;
          lookup_.insert(s);
        }
      } else { throw std::runtime_error("Error: empty part"); }
    }

    ~Part() = default;

    /**
     * @brief Get a constant reference to the underlying mesh.
     *
     * @return a constant reference to the underlying mesh.
     */
    const Mesh& mesh() const { return mesh_; }

    /**
     * @brief Get a constant reference to the underlying state.
     *
     * @return a constant reference to the underlying state.
     */
    const State& state() const { return state_; }

    /**
     * @brief Get a normal reference to the underlying state.
     *
     * @return a reference to the underlying state.
     */
    State& state() { return const_cast<State&>(state_); }

    /**
     * @brief Get a reference to the entities list.
     *
     * @return a reference to entities list.
     */
    std::vector<int> const& entities() const { return entities_; }

    /**
     * @brief Check if a given entity is in the part.
     *
     * @param id entity ID
     * @return true if so, false otherwise.
     */
    bool contains(int id) const { return lookup_.count(id) == 1; }

    /**
     * @brief Retrieve relative index of given entity.
     *
     * @param id: entity absolute index in mesh.
     * @return entity relative index in part.
     */
    const int& index(int id) const { return index_.at(id); }

    /**
     * @brief Get part size.
     *
     * @return entities list size
     */
    const int& size() const { return size_; }

    /**
     * @brief Retrieve the field data of the part.
     *
     * @param variable: name of the field.
     * @param data: the field data.
     */
    void get_field(std::string variable, double** data) const {
      return state_.mesh_get_data(onwhat, variable, data);
    }

    /**
     * @brief Retrieve the field data of the part.
     *
     * @param variable: name of the field.
     * @param data: the field data.
     */
    void get_field(std::string variable, const double** data) const {
      return state_.mesh_get_data(onwhat, variable, data);
    }

    /**
     * @brief Get the volume of the given entity.
     *
     * @param id: the relative index of the entity.
     * @return its volume.
     */
    const double& volume(int id) const {
      assert(id < size_);
      return volumes_[id];
    }

    /**
     * @brief Retrieve the neighbors of the given entity on mesh part.
     *
     * @tparam entity_type the entity type [ALL|PARALLEL_OWNED]
     * @param  entity      the given entity
     * @return filtered    the filtered neighboring entities list.
     */
    template<Entity_type entity_type = Entity_type::ALL>
    std::vector<int> get_neighbors(int entity) const {
      std::vector<int> neighbors, filtered;
      // first retrieve neighbors
      switch (onwhat) {
        case CELL: mesh_.cell_get_node_adj_cells(entity, entity_type, &neighbors); break;
        case NODE: mesh_.node_get_cell_adj_nodes(entity, entity_type, &neighbors); break;
        default: std::cerr << "Error: unsupported entity kind" << std::endl; break;
      }

      // filter then
      filtered.reserve(neighbors.size());
      for (auto const& current : neighbors) {
        if (lookup_.count(current)) {
          filtered.emplace_back(current);
        }
      }
      return filtered;
    }

    /**
     * @brief Compute the total volume of the part.
     *
     * @return the total volume of the part.
     */
    double compute_total_volume() const {
      return std::accumulate(volumes_.begin(), volumes_.end(), 0.0);
    }

    /**
     * @brief Compute entities volumes within the part.
     *
     * @param masks: entities masks to disable some of them.
     * @return the total volume of the part.
     */
    double compute_entities_volumes(int const* masks = nullptr) {
      if (not cached_volumes) {
        // check if entities mask should be used
        bool const use_masks = masks != nullptr;
        bool const on_cell = onwhat == CELL;
        // kernel to compute the volume of an entity
        auto compute_volume = [&](int s) {
          double volume = (on_cell ? mesh_.cell_volume(s) : mesh_.dual_cell_volume(s));
          volumes_[index_[s]] = (use_masks ? masks[s] * volume : volume);
        };
        // apply kernel on all entities of the part
        Portage::for_each(entities_.begin(), entities_.end(), compute_volume);
        // toggle flag
        cached_volumes = true;
      }
      // finally accumulate them to retrieve the global volume
      return compute_total_volume();
    }

  private:
    // mesh and state
    Mesh const& mesh_;
    State& state_;

    int size_ = 0;
    bool cached_volumes = false;

    /* Part meta-data:
     * - entities list, related volumes and relative indices.
     * - a hashtable to have constant-time parts lookup queries in average case.
     *   remark: for lookup purposes only, not meant to be iterated. */
    std::vector<int>        entities_ {};
    std::vector<double>     volumes_  {};
    std::map<int, int>      index_    {};
    std::unordered_set<int> lookup_   {};

    // get rid of long namespaces
    constexpr static auto const CELL = Wonton::Entity_kind::CELL;
    constexpr static auto const NODE = Wonton::Entity_kind::NODE;
  };


/**
 * @brief Manages source and target sub-meshes for part-by-part remap.
 *        It detects boundaries mismatch and provides the necessary fixup
 *        for partially filled and empty cells values.
 *
 * @tparam D           meshes dimension
 * @tparam onwhat      the entity kind for remap [cell|node]
 * @tparam SourceMesh  the source mesh wrapper to use
 * @tparam SourceState the source state wrapper to use
 * @tparam TargetMesh  the target mesh wrapper to use
 * @tparam TargetState the target state wrapper to use
 */
  template<int D, Entity_kind onwhat,
    class SourceMesh, class SourceState,
    class TargetMesh = SourceMesh,
    class TargetState = SourceState
  >
class PartPair {

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
    SourceMesh const& source_mesh, SourceState& source_state,
    TargetMesh const& target_mesh, TargetState& target_state,
    std::vector<int> const& source_entities,
    std::vector<int> const& target_entities,
    Wonton::Executor_type const* executor
  ) : source_(source_mesh, source_state, source_entities),
      target_(target_mesh, target_state, target_entities)
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

    // Get info about which entities on this processor should be
    // masked out and not accounted for in calculations because they
    // were already encountered on a lower rank processor. We have to
    // do this because in our distributed runs, our source partitions
    // don't form a strict tiling (no overlaps) after redistribution
    int nb_masks = (onwhat == Entity_kind::CELL ? source_.mesh().num_owned_cells()
                                                : source_.mesh().num_owned_nodes());

    source_masks_.resize(nb_masks, 1);
#ifdef PORTAGE_ENABLE_MPI
    if (distributed_) {
      get_unique_entity_masks<onwhat, SourceMesh>(
        source_.mesh(), &source_masks_, mycomm_
      );
    }
#endif

    intersect_volumes_.resize(target_.size());
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
   * @brief Retrieve a constant pointer to source mesh part.
   *
   * @return a constant pointer to source mesh part
   */
  Part<onwhat, SourceMesh, SourceState> const* get_source() const { return &source_; }

  /**
   * @brief Retrieve a constant pointer to target mesh part.
   *
   * @return a constant pointer to target mesh part
   */
  Part<onwhat, TargetMesh, TargetState> const* get_target() const { return &target_; }

  /**
   * @brief Get a reference to source entities list.
   *
   * @return a reference to source entities list.
   */
  std::vector<int> const& get_source_entities() const { return source_.entities(); }

  /**
   * @brief Get a reference to target entities list.
   *
   * @return a reference to source entities list.
   */
  std::vector<int> const& get_target_entities() const { return target_.entities(); }

  /**
   * @brief Check if a given entity is in source part list.
   *
   * @param id entity ID
   * @return true if so, false otherwise.
   */
  bool is_source_entity(int id) const { return source_.contains(id); }

  /**
   * @brief Check if a given entity is in target part list.
   *
   * @param id entity ID
   * @return true if so, false otherwise.
   */
  bool is_target_entity(int id) const { return target_.contains(id); }

  /**
   * @brief Get source part size.
   *
   * @return source entities list size
   */
  const int& source_part_size() const { return source_.size(); }

  /**
   * @brief Get target part size.
   *
   * @return target entities list size
   */
  const int& target_part_size() const { return target_.size(); }

  /**
   * @brief Retrieve the neighbors of the given entity on source mesh.
   *
   * @tparam entity_type the entity type [ALL|PARALLEL_OWNED]
   * @param  entity      the given entity
   * @return filtered    the filtered neighboring entities list.
   */
  template<Entity_type entity_type = Entity_type::ALL>
  std::vector<int> get_source_filtered_neighbors(int entity) const {
    return source_.get_neighbors(entity);
  }

  /**
   * @brief Retrieve the neighbors of the given entity on target mesh.
   *
   * @tparam entity_type the entity type [ALL|PARALLEL_OWNED]
   * @param  entity      the given entity
   * @return filtered    the filtered neighboring entities list.
   */
  template<Entity_type entity_type = Entity_type::ALL>
  std::vector<int> get_target_filtered_neighbors(int entity) const {
    return target_.get_neighbors(entity);
  }

  /**
   * @brief Compute source and target parts intersection volume.
   *
   * @param source_weights: candidate source cells and their intersection moments.
   * @return the total intersection volume.
   */
  double compute_intersect_volumes
    (Portage::vector<entity_weights_t> const& source_weights) {
    // retrieve target entities list
    auto const& target_entities = target_.entities();

    // compute the intersected volume of each target part entity
    Portage::for_each(target_entities.begin(), target_entities.end(), [&](int t) {
      auto const& i = target_.index(t);
      // accumulate moments
      entity_weights_t const& moments = source_weights[t];
      intersect_volumes_[i] = 0.;
      for (auto const& current : moments) {
        // matched source cell should be in the source part
        if (source_.contains(current.entityID))
          intersect_volumes_[i] += current.weights[0];
        #if DEBUG_PART_BY_PART
          std::printf("\tmoments[target:%d][source:%d]: %f\n"
                      , t, current.entityID, current.weights[0]);
        #endif
      }
      #if DEBUG_PART_BY_PART
        std::printf("intersect_volume[%02d]: %.3f\n", t, intersect_volumes_[i]);
      #endif
    });

    // accumulate values to retrieve total intersected volume
    return std::accumulate(intersect_volumes_.begin(), intersect_volumes_.end(), 0.);
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
  bool check_mismatch(Portage::vector<entity_weights_t> const& source_weights) {

    // ------------------------------------------
    // COMPUTE VOLUMES ON SOURCE AND TARGET PARTS
    // ------------------------------------------
    // collect volumes of entities that are not masked out and sum them up
    double source_volume = source_.compute_entities_volumes(source_masks_.data());
    double target_volume = target_.compute_entities_volumes();
    double intersect_volume = compute_intersect_volumes(source_weights);

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
    empty_entities.reserve(target_.size());

    is_cell_empty_.resize(target_.size(), false);

    for (auto&& entity : target_.entities()) {
      auto const& i = target_.index(entity);
      if (std::abs(intersect_volumes_[i]) < epsilon_) {
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
      layer_num_.resize(target_.size(), 0);

      int nb_layers = 0;
      int nb_tagged = 0;
      int old_nb_tagged = -1;

      while (nb_tagged < nb_empty and nb_tagged > old_nb_tagged) {
        old_nb_tagged = nb_tagged;

        std::vector<int> current_layer_entities;

        for (auto&& entity : empty_entities) {
          auto const& i = target_.index(entity);
          // skip already set entities
          if (layer_num_[i] == 0) {

            auto neighbors = get_target_filtered_neighbors<Entity_type::ALL>(entity);

            for (auto&& neigh : neighbors) {
              auto const& j = target_.index(neigh);
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
          auto const& t = target_.index(entity);
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
                    Empty_fixup_type empty_fixup_type = EXTRAPOLATE) const {

    if (source_.state().field_type(onwhat, src_var_name) == Field_type::MESH_FIELD) {
      return fix_mismatch_meshvar(src_var_name, trg_var_name,
                                  global_lower_bound, global_upper_bound,
                                  conservation_tol, maxiter,
                                  partial_fixup_type, empty_fixup_type);
    } else
      return false;
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
                            Empty_fixup_type empty_fixup_type) const {

    // valid only for part-by-part scenario
    static bool hit_lower_bound  = false;
    static bool hit_higher_bound = false;

    // use aliases
    auto const& source_entities = source_.entities();
    auto const& target_entities = target_.entities();

    // Now process remap variables
    // WARNING: absolute indexing
    double const* source_data;
    double*       target_data;
    source_.get_field(src_var_name, &source_data);
    target_.get_field(trg_var_name, &target_data);

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

      for (auto&& entity : target_entities) {
        auto const& t = target_.index(entity);
        if (not is_cell_empty_[t]) {

          #if DEBUG_PART_BY_PART
            std::printf("fixing target_data[%d] with locally conservative fixup\n", entity);
            std::printf("= before: %.3f", target_data[entity]);
          #endif

          auto const relative_voldiff =
            std::abs(intersect_volumes_[t] - target_.volume(t)) / target_.volume(t);

          if (relative_voldiff > tolerance_) {
            target_data[entity] *= intersect_volumes_[t] / target_.volume(t);
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
            auto const& i = target_.index(neigh);
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

      for (auto&& s : source_entities) {
        auto const& i = source_.index(s);
        source_sum += source_data[s] * source_.volume(i);
      }

      for (auto&& t : target_entities) {
        auto const& i = target_.index(t);
        target_sum += target_data[t] * target_.volume(i);
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
        for (auto&& entity : target_entities) {
          auto const& t = target_.index(entity);
          if (not is_cell_empty_[t]) {
            covered_target_volume += target_.volume(t);
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
        adj_target_volume = target_.compute_total_volume();
        global_adj_target_volume = global_target_volume_;
      }

      // get the right entity type
      auto target_entity_type = [&](int entity) -> Entity_type {
        return (onwhat == Entity_kind::CELL ? target_.mesh().cell_get_type(entity)
                                            : target_.mesh().node_get_type(entity));
      };

      // sort of a "unit" discrepancy or difference per unit volume
      double udiff = absolute_diff / global_adj_target_volume;

      int iter = 0;
      while (std::abs(relative_diff) > conservation_tol and iter < maxiter) {

        for (auto&& entity : target_entities) {
          auto const& t = target_.index(entity);
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
              adj_target_volume -= target_.volume(t);

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
              adj_target_volume -= target_.volume(t);
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
        for (auto&& entity : target_entities) {
          auto const& t = target_.index(entity);
          target_sum += target_.volume(t) * target_data[entity];
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


private:
  // source and target mesh parts
  Part<onwhat, SourceMesh, SourceState> source_;
  Part<onwhat, TargetMesh, TargetState> target_;

  bool is_mismatch_tested_ = false;
  bool has_mismatch_       = false;

  // useful constants
  static constexpr double infinity_  = std::numeric_limits<double>::max();
  static constexpr double epsilon_   = std::numeric_limits<double>::epsilon();
  static constexpr double tolerance_ = 1.E2 * epsilon_;

  // data needed for mismatch checks
  double global_source_volume_    = 0.;
  double global_target_volume_    = 0.;
  double global_intersect_volume_ = 0.;
  double relative_voldiff_        = 0.;

  std::vector<int>    source_masks_   = {};
  std::vector<double> intersect_volumes_ = {};
  // empty target cells management
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
