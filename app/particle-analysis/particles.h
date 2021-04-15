/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

#include <random>
#include "params.h"
#include "formula.h"

namespace particles {

/**
 * @brief Generate points.
 *
 * @param num: number of points.
 * @param distrib: statistical distribution.
 * @param p_min: lower bound
 * @param p_max: upper bound
 * @param scaled_span: scaled span angle for radial grid.
 * @param scaled_width: scaled width for radial grid.
 * @param radius: distance from origin for radial grid.
 * @param center: center of frame for radial grid.
 * @param scaling: scaling factors for coordinates.
 * @param frame: coordinate system.
 * @param seed: user seed for random points.
 * @param output: output file if desired.
 * @return a swarm.
 */
Wonton::Swarm<2> generate(const int num[2], Distrib distrib,
                          Wonton::Point<2> const& p_min,
                          Wonton::Point<2> const& p_max,
                          double scaled_span = 0,
                          double scaled_width = 0,
                          double radius = 0,
                          Wonton::Point<2> const& center = { 0, 0 },
                          Wonton::Point<2> const& scaling = { 1, 1 },
                          Frame frame = LOCAL,
                          unsigned seed = 0,
                          std::string const& output = "") {

  std::vector<Wonton::Point<2>> points;
  int const num_points = num[0] * num[1];

  // generate a radial mesh
  if (distrib == RADIAL) {
    double center_angle = 0.;
    Wonton::Point<2> reference(0.0, 0.0);

    // NOTE: this algorithm works by determining the outer mesh point for an angle,
    // move radially inward to determine each inner mesh point, then repeat for the next angle.
    // locate the topleft grid point from the center point of mesh
    // for local (x',y') frame, use (x,y) = (0,1) as reference
    switch (frame) {
      case LOCAL: center_angle = 0.5 * M_PI; reference[1] = radius; break;
      case GLOBAL: center_angle = std::atan2(center[1], center[0]); break;
      default: throw std::runtime_error("invalid coords frame");
    }

    double const span = scaled_span / scaling[0];
    double const width = scaled_width / scaling[1];
    double const unit_angle = span / (num[0] - 1);
    double const unit_width = width / (num[1] - 1);
    double const half_width = 0.5 * width;
    double const half_angle = 0.5 * span;
    double const cos_unit_angle = std::cos(unit_angle);
    double const sin_unit_angle = std::sin(unit_angle);

    double top_x = (radius + half_width) * cos(center_angle);
    double top_y = (radius + half_width) * sin(center_angle);
    double topleft_x = top_x * cos(half_angle) - top_y * sin(half_angle);
    double topleft_y = top_x * sin(half_angle) + top_y * cos(half_angle);

    points.resize(num_points);

    // locate the grid points from 'topleft' vertically
    // and rotate the 'topleft' grid point
    // move from small angle to large angle
    for (int i = 0; i < num[0]; ++i) {
      int k = i * num[1];
      // update dx and dy for current radial direction
      double const dx = unit_width * topleft_x / (radius + half_width);
      double const dy = unit_width * topleft_y / (radius + half_width);
      // set top mesh point for each angle
      points[k][0] = topleft_x - reference[0];
      points[k][1] = topleft_y - reference[1];

      // move radially inward for a specific angle from top to bottom
      for (int j = k; j < num[1] + k; ++j) {
        points[j + 1][0] = points[j][0] - dx;
        points[j + 1][1] = points[j][1] - dy;
      }

      top_x = topleft_x;
      top_y = topleft_y;
      // rotate to calculate top mesh point for next angle
      topleft_x =  top_x * cos_unit_angle + top_y * sin_unit_angle;
      topleft_y = -top_x * sin_unit_angle + top_y * cos_unit_angle;
    }
  } else {

    double const& x_min = p_min[0];
    double const& y_min = p_min[1];
    double const& x_max = p_max[0];
    double const& y_max = p_max[1];

    if (distrib == CARTESIAN) {
      double const dx = (x_max - x_min) / (num[0] - 1);
      double const dy = (y_max - y_min) / (num[1] - 1);

      for (int i = 0; i < num[0]; ++i) {
        double const x = x_min + i * dx;
        for (int j = 0; j < num[1]; ++j) {
          double const y = y_min + j * dy;
          points.emplace_back(x, y);
        }
      }
    } else if (distrib == RANDOM) {
      // set the random engine and generator
      std::random_device device;
      std::mt19937 engine{seed ? seed : device()};
      std::uniform_real_distribution<double> x_distrib(x_min, x_max);
      std::uniform_real_distribution<double> y_distrib(y_min, y_max);

      for (int i = 0; i < num_points; ++i) {
        points.emplace_back(x_distrib(engine), y_distrib(engine));
      }
    }

    // rescale
    bool const rescale = std::abs(scaling[0] - 1.) > 0
                      or std::abs(scaling[1] - 1.) > 0;
    if (rescale) {
      for (auto& p : points) {
        for (int d = 0; d < 2; ++d) {
          p[d] /= scaling[d];
        }
      }
    }
  }

  if (not output.empty()) {
    std::ofstream file(output);
    if (file.good()) {
      for (auto&& p : points) {
        file << p[0] << "," << p[1] << std::endl;
      }
    } else {
      throw std::runtime_error("failed to open '"+ output +"'");
    }
  }

  return Wonton::Swarm<2>(points);
}

/**
 * @brief Compute mean distance between points.
 *
 * @param swarm: the set of points
 * @param scale: coordinates scaling factors
 * @param distrib: statistical distribution
 * @param p_min: lower bound
 * @param p_max: upper bound
 * @param scaling: scaling factors for coordinates.
 * @param radius: offset from origin for radial grid
 * @param scaled_span: scaled span angle for radial grid.
 * @param scaled_width: scaled width for radial grid.
 * @param size: number of points in each direction.
 * @return the mean distance in each direction.
 */
std::vector<double> mean_distance(Wonton::Swarm<2> const& swarm,
                                  Wonton::Point<2> const& scale = { 1, 1 },
                                  Distrib distrib = Distrib::RADIAL,
                                  Wonton::Point<2> const& p_min = { 0, 0 },
                                  Wonton::Point<2> const& p_max = { 0, 0 },
                                  Wonton::Point<2> const& scaling = { 1, 1 },
                                  double radius = 1,
                                  double scaled_span = 0,
                                  double scaled_width = 0,
                                  const int size[2] = nullptr) {

  if (size == nullptr or size[0] < 1 or size[1] < 1) {
    throw std::runtime_error("invalid swarm size");
  }

  std::vector<double> h(2, 0.);

  if (distrib == RADIAL) {
    double const span = scaled_span / scaling[0];
    double const width = scaled_width / scaling[1];
    h[0] = scale[0] * radius * span / (size[0] - 1);
    h[1] = scale[1] * width / (size[1] - 1);
  } else {
    for (int d = 0; d < 2; ++d) {
      h[d] = scale[d] * (p_max[d] - p_min[d]) / (size[d] - 1);
    }
  }

  return h;
}


/**
 * @brief Import points.
 *
 * @param path: file to load
 * @param scaling: scaling factors for coordinates.
 * @param unscale: unscale coordinates instead.
 * @return a swarm
 */
Wonton::Swarm<2> import(std::string const& path,
                        Wonton::Point<2>& p_min,
                        Wonton::Point<2>& p_max,
                        int size[2],
                        Wonton::Point<2> const& scaling = { 1, 1 },
                        bool unscale = true) {

  std::ifstream file(path);
  if (not file.good()) { throw std::runtime_error("cannot open "+ path); }

  assert(size != nullptr);

  std::string line;
  double p[2];
  double const s[] = { unscale ? 1./scaling[0] : scaling[0],
                       unscale ? 1./scaling[1] : scaling[1] };

  std::vector<Wonton::Point<2>> points;

  char comma;
  while (std::getline(file, line)) {
    std::istringstream buffer(line);
    buffer >> p[0] >> comma >> p[1];
    points.emplace_back(p[0] * s[0], p[1] * s[1]);
  }

  // deduce coordinates extents
  for (auto&& p_cur : points) {
    for (int d = 0; d < 2; ++d) {
      if (p_cur[d] < p_min[d]) { p_min[d] = p_cur[d]; }
      if (p_cur[d] > p_max[d]) { p_max[d] = p_cur[d]; }
    }
  }

  // deduce number of points per axis
  int const num_points = points.size();
  size[0] = size[1] = static_cast<int>(std::sqrt(num_points));

  return Wonton::Swarm<2>(points);
}

/**
 * @brief Initialize a swarm.
 *
 * @param path: file to load if any
 * @param num: number of points
 * @param distrib: their distribution
 * @param p_min: lower bound
 * @param p_max: upper bound
 * @param span: scaled span angle for radial grid.
 * @param width: scaled width for radial grid.
 * @param radius: distance from origin for radial grid.
 * @param center: center of frame for radial grid.
 * @param scaling: scaling factors for coordinates.
 * @param frame: coordinate system.
 * @param seed: user seed for random points.
 * @param output: output file if desired.
 * @return a swarm.
 */
Wonton::Swarm<2> init(std::string const& path,
                      int num[2],
                      Distrib distrib,
                      Wonton::Point<2>& p_min,
                      Wonton::Point<2>& p_max,
                      double span = 0,
                      double width = 0,
                      double radius = 0,
                      Wonton::Point<2> const& center = { 0,0 },
                      Wonton::Point<2> const& scaling = { 1,1 },
                      Frame frame = LOCAL,
                      unsigned seed = 0,
                      std::string const& output = "") {

  return not path.empty() ? import(path, p_min, p_max, num, scaling)
                          : generate(num, distrib, p_min, p_max, span, width,
                                     radius, center, scaling, frame, seed, output);
}


} // namespace 'particles'