/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "wonton/support/wonton.h"
#include "wonton/swarm/swarm.h"
#include "wonton/swarm/swarm_state.h"
#include "portage/support/portage.h"
#include "portage/driver/driver_swarm.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/search/search_points_bins.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

#include "particles.h"
#include "analysis.h"
#include "steps.h"

// aliases
constexpr int dim = 2;
using Matrix = std::vector<std::vector<double>>;

// driver using legacy search kernel for both weight forms.
using RemapLegacy = Portage::SwarmDriver<Portage::SearchPointsByCells,
                                         Portage::Accumulate,
                                         Portage::Estimate,
                                         dim, Wonton::Swarm<dim>,
                                         Wonton::SwarmState<dim>>;

// driver using lightweight search kernel for gather form.
using RemapGather = Portage::SwarmDriver<Portage::SearchPointsBins,
                                         Portage::Accumulate,
                                         Portage::Estimate,
                                         dim, Wonton::Swarm<dim>,
                                         Wonton::SwarmState<dim>>;

/**
 * @brief Run the particle remap analysis application.
 *
 * It aims to ease numerical testing on particle remap by providing means
 * to run specific remap scenarii and to post-process the results using
 * the notebook. The input parameters are given in a JSON file. The source
 * and target points can be imported or generated using a cartesian, radial
 * or uniformly random pattern. The source field is given as an analytical
 * function to be parsed. For now, only the gather form is supported for
 * the remap weights which means that they are assigned to target points.
 * For kinetic beam tests, coordinates and field are scaled on a power of the
 * Lorentz factor gamma, and the smoothing lengths can be rescaled by
 * specifying related factors.
 *
 * @param argc: arguments count
 * @param argv: arguments values
 * @return status code
 */
int main(int argc, char* argv[]) {

  Params params;
  if (not params.parse(argc, argv))
    return EXIT_FAILURE;

  Step step(params);
  /* ------------------------------------------------------------------------ */
  step.start("Initialize");

  auto executor = params.get_executor();

  auto source_swarm = particles::init(params.source.file,
                                      params.source.size, params.source.distrib,
                                      params.source.min, params.source.max,
                                      params.source.span, params.source.width,
                                      params.source.radius, params.source.center,
                                      params.scale.coords, params.source.frame);

  auto target_swarm = particles::init(params.target.file,
                                      params.target.size, params.target.distrib,
                                      params.target.min, params.target.max,
                                      params.target.span, params.target.width,
                                      params.target.radius, params.target.center,
                                      params.scale.coords, params.target.frame,
                                      0, params.target.output);

#ifdef DEBUG
  std::cout << " ===== source ===== " << std::endl;
  for (int i = 0; i < 20; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    std::printf("[%f, %f]\n", p[0], p[1]);
  }

  std::cout << " ===== target ===== " << std::endl;
  for (int i = 0; i < 20; ++i) {
    auto p = target_swarm.get_particle_coordinates(i);
    std::printf("[%f, %f]\n", p[0], p[1]);
  }
#endif

  Wonton::SwarmState<dim> source_state(source_swarm);
  Wonton::SwarmState<dim> target_state(target_swarm);

  // add source field
  auto source_field = analysis::eval(source_swarm, params.source.field,
                                     params.scale.coords, params.scale.field,
                                     params.output.source);

  source_state.add_field("field", source_field);
  target_state.add_field<double>("field", 0.0);

  // set smoothing lengths
  int const num_source = source_swarm.num_particles();
  int const num_target = target_swarm.num_particles();
  bool const gather_form = (params.remap.weight == Portage::WeightCenter::Gather);

  Wonton::vector<Matrix> smoothing_lengths(gather_form ? num_target : num_source);

  auto const h = particles::mean_distance(target_swarm,
                                          params.scale.smooth, params.target.distrib,
                                          params.target.min, params.target.max,
                                          params.scale.coords, params.target.radius,
                                          params.target.span, params.target.width,
                                          params.target.size);

  Wonton::transform(smoothing_lengths.begin(),
                    smoothing_lengths.end(),
                    smoothing_lengths.begin(),
                    [&](auto& m) { return Matrix(1, h); });

  /* ------------------------------------------------------------------------ */
  step.start("Remap");

  if (gather_form) {
    RemapGather remapper(source_swarm, source_state,
                         target_swarm, target_state,
                         smoothing_lengths, params.remap.kernel,
                         params.remap.support, params.remap.weight);

    remapper.set_remap_var_names({"field"}, {"field"},
                                 params.remap.estimator, params.remap.basis);
    remapper.run(executor);
  } else {
    RemapLegacy remapper(source_swarm, source_state,
                         target_swarm, target_state,
                         smoothing_lengths, params.remap.kernel,
                         params.remap.support, params.remap.weight);

    remapper.set_remap_var_names({"field"}, {"field"},
                                 params.remap.estimator, params.remap.basis);
    remapper.run(executor);
  }

  /* ------------------------------------------------------------------------ */
  step.start("Assess");

  int const norm = 2;
  auto target_remap = target_state.get_field("field");
  auto target_exact = analysis::eval(target_swarm, params.source.field,
                                     params.scale.coords, params.scale.field);
  auto target_error = analysis::error(target_exact, target_remap, norm,
                                      params.output.error,
                                      params.output.exact,
                                      params.output.remap);

  step.stop();

  if (params.mpi.rank == 0) {
    std::printf("unscaled error L%d: %.2f\n", norm, target_error / params.scale.field);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
