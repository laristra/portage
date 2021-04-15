/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#pragma once

#include <iomanip>
#include "params.h"

namespace analysis {

/**
 * @brief Assess analytical field values on given points.
 *
 * @param points: the given set of points.
 * @param function: the analytical expression to evaluate.
 * @param coord_scale: coordinates scaling factors.
 * @param field_scale: field scaling factor.
 * @param output: the output file path.
 * @return field values.
 */
Wonton::vector<double> eval(Wonton::Swarm<2> const& points,
                            std::string const& function,
                            Point<2> const& coord_scale = { 1, 1 },
                            double field_scale = 1,
                            std::string const& output = "") {
  Formula formula;
  formula.initialize(function);

  int const num_points = points.num_particles();
  Wonton::vector<double> values(num_points);
  double const s[] = { coord_scale[0],
                       coord_scale[1],
                       field_scale };

  for (int i = 0; i < num_points; ++i) {
    auto const p = points.get_particle_coordinates(i);
    values[i] = formula({ p[0] * s[0], p[1] * s[1]}) * s[2];
  }

  if (not output.empty()) {
    std::ofstream file(output);
    if (file.good()) {
      for (int i = 0; i < num_points; ++i) {
        file << std::setprecision(12) << values[i] << std::endl;
      }
    } else {
      throw std::runtime_error("failed to open '"+ output +"'");
    }
  }

  return values;
}

/**
 * @brief Assess error of approximated field.
 *
 * @param exact: the exact field values.
 * @param approx: the approximateed field values.
 * @param norm: the norn to consider.
 * @param output: the output file path.
 * @return error field.
 */
double error(Wonton::vector<double> const& exact,
             Wonton::vector<double> const& approx,
             int norm,
             std::string const& output_error,
             std::string const& output_exact,
             std::string const& output_approx) {

  assert(exact.size() == approx.size());

  // step 1: compute error
  int const num_values = exact.size();
  double error[num_values];
  double error_norm = 0.0;

  // store pointwise error
  for (int i = 0; i < num_values; ++i) {
    double const delta = std::abs(exact[i] - approx[i]);
    error[i] = (not std::isnan(delta) ? delta : 0.);
  }

  // compute error in specified norm
  switch (norm) {
    case 0:  for (auto&& value : error) { error_norm = std::max(value, error_norm); } break;
    case 1:  for (auto&& value : error) { error_norm += value; } break;
    case 2:  for (auto&& value : error) { error_norm += value * value; } error_norm = sqrt(error_norm); break;
    default: for (auto&& value : error) { error_norm += std::pow(value, norm); } error_norm = std::pow(error_norm, 1./norm); break;
  }

  // step 2: export if requested
  if (not output_error.empty()) {
    std::ofstream file(output_error);
    if (file.good()) {
      for (int i = 0; i < num_values; ++i) {
        file << std::setprecision(12) << error[i] << std::endl;
      }
    } else {
      throw std::runtime_error("failed to open '"+ output_error +"'");
    }
  }

  if (not output_exact.empty()) {
    std::ofstream file(output_exact);
    if (file.good()) {
      for (int i = 0; i < num_values; ++i) {
        file << std::setprecision(12) << exact[i] << std::endl;
      }
    } else {
      throw std::runtime_error("failed to open '"+ output_exact +"'");
    }
  }

  if (not output_approx.empty()) {
    std::ofstream file(output_approx);
    if (file.good()) {
      for (int i = 0; i < num_values; ++i) {
        file << std::setprecision(12) << approx[i] << std::endl;
      }
    } else {
      throw std::runtime_error("failed to open '"+ output_approx +"'");
    }
  }

  return error_norm;
}

} // namespace 'analysis'
