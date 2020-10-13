/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

#include "params.h"

using T = double;
using Portage::Interpolate_1stOrder;
using Portage::Interpolate_2ndOrder;

/**
 * @brief Process part-by-part remapping for the given field.
 *
 * @tparam dim: source/target meshes dimension.
 * @param field: a string expression of the numerical field.
 * @param nb_parts: the number of source-target parts couples.
 * @param source_mesh: a pointer to the source mesh.
 * @param target_mesh: a pointer to the target mesh.
 * @param source_mesh_wrapper: a wrapper to access source mesh data.
 * @param target_mesh_wrapper: a wrapper to access target mesh data.
 * @param source_state_wrapper: a wrapper to access source state data.
 * @param target_state_wrapper: a wrapper to access source state data.
 * @param executor: a pointer to the MPI executor.
 * @param source_cells: list of source cells for each part.
 * @param target_cells: list of target cells for each part.
 * @param params: input parameters
 */
template <int dim>
void remap(std::string const& field, int nb_parts,
           std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh,
           Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
           Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
           Wonton::Jali_State_Wrapper& source_state_wrapper,
           Wonton::Jali_State_Wrapper& target_state_wrapper,
           Wonton::Executor_type* executor,
           std::vector<std::vector<int>> const& source_cells,
           std::vector<std::vector<int>> const& target_cells,
           Params const& params);

/**
 * @brief Process part-by-part remapping for the given 2D field.
 *
 * @param field: a string expression of the numerical field.
 * @param nb_parts: the number of source-target parts couples.
 * @param source_mesh: a pointer to the source mesh.
 * @param target_mesh: a pointer to the target mesh.
 * @param source_mesh_wrapper: a wrapper to access source mesh data.
 * @param target_mesh_wrapper: a wrapper to access target mesh data.
 * @param source_state_wrapper: a wrapper to access source state data.
 * @param target_state_wrapper: a wrapper to access source state data.
 * @param executor: a pointer to the MPI executor.
 * @param source_cells: list of source cells for each part.
 * @param target_cells: list of target cells for each part.
 * @param params: input parameters
 */
template<>
void remap<2>(std::string const& field, int nb_parts,
              std::shared_ptr<Jali::Mesh> source_mesh,
              std::shared_ptr<Jali::Mesh> target_mesh,
              Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
              Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
              Wonton::Jali_State_Wrapper& source_state_wrapper,
              Wonton::Jali_State_Wrapper& target_state_wrapper,
              Wonton::Executor_type* executor,
              std::vector<std::vector<int>> const& source_cells,
              std::vector<std::vector<int>> const& target_cells,
              Params const& params) {

  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Wonton::Jali_State_Wrapper>;

  using PartPair = Portage::PartPair<2, Wonton::Jali_Mesh_Wrapper,
                                     Wonton::Jali_State_Wrapper>;


  std::vector<PartPair> parts_registry;
  parts_registry.reserve(nb_parts);

  for (int i = 0; i < nb_parts; ++i) {
    // create source-target mesh parts manager and
    // populate cell lists for the current part.
    parts_registry.emplace_back(source_mesh_wrapper, source_state_wrapper,
                                target_mesh_wrapper, target_state_wrapper,
                                source_cells[i], target_cells[i], executor);
  }

  // perform remap kernels.
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectRnD>(candidates);

  for (auto&& part : parts_registry) {

    if (params.order == 2) {
      auto source_part = part.source();
      auto gradients = remapper.compute_source_gradient(field, params.limiter,
                                                        params.bnd_limiter,0,
                                                        &source_part);
      remapper.interpolate_mesh_var<T, Interpolate_2ndOrder>(field, field, weights,
                                                             &part, &gradients);
    } else {
      remapper.interpolate_mesh_var<T, Interpolate_1stOrder>(field, field, weights, &part);
    }

    if (part.check_mismatch(weights)) {
      part.fix_mismatch(field, field, params.lower_bound, params.upper_bound,
                        params.tolerance, params.fix_iter,
                        params.partial_fixup, params.empty_fixup);
    }
  }
}

/**
 * Process part-by-part remapping for the given 3D field.
 *
 * @param field: a string expression of the numerical field.
 * @param nb_parts: the number of source-target parts couples.
 * @param source_mesh: a pointer to the source mesh.
 * @param target_mesh: a pointer to the target mesh.
 * @param source_mesh_wrapper: a wrapper to access source mesh data.
 * @param target_mesh_wrapper: a wrapper to access target mesh data.
 * @param source_state_wrapper: a wrapper to access source state data.
 * @param target_state_wrapper: a wrapper to access source state data.
 * @param executor: a pointer to the MPI executor.
 * @param source_cells: list of source cells for each part.
 * @param target_cells: list of target cells for each part.
 * @param params: input parameters
 */
template<>
void remap<3>(std::string const& field, int nb_parts,
              std::shared_ptr<Jali::Mesh> source_mesh,
              std::shared_ptr<Jali::Mesh> target_mesh,
              Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
              Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
              Wonton::Jali_State_Wrapper& source_state_wrapper,
              Wonton::Jali_State_Wrapper& target_state_wrapper,
              Wonton::Executor_type* executor,
              std::vector<std::vector<int>> const& source_cells,
              std::vector<std::vector<int>> const& target_cells,
              Params const& params) {

  using Remapper = Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
    Wonton::Jali_Mesh_Wrapper,
    Wonton::Jali_State_Wrapper>;

  using PartPair = Portage::PartPair<3, Wonton::Jali_Mesh_Wrapper,
    Wonton::Jali_State_Wrapper>;

  std::vector<PartPair> parts_registry;
  parts_registry.reserve(nb_parts);

  // filter cells and populate lists
  for (int i = 0; i < nb_parts; ++i) {
    parts_registry.emplace_back(source_mesh_wrapper, source_state_wrapper,
                                target_mesh_wrapper, target_state_wrapper,
                                source_cells[i], target_cells[i], executor);
  }

  // perform remap kernels.
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectRnD>(candidates);

  for (auto&& part : parts_registry) {

    if (params.order == 2) {
      auto source_part = part.source();
      auto gradients = remapper.compute_source_gradient(field, params.limiter,
                                                        params.bnd_limiter,0,
                                                        &source_part);
      remapper.interpolate_mesh_var<T, Interpolate_2ndOrder>(field, field, weights,
                                                             &part, &gradients);
    } else {
      remapper.interpolate_mesh_var<T, Interpolate_1stOrder>(field, field, weights, &part);
    }

    if (part.check_mismatch(weights)) {
      part.fix_mismatch(field, field, params.lower_bound, params.upper_bound,
                        params.tolerance, params.fix_iter,
                        params.partial_fixup, params.empty_fixup);
    }
  }
}
