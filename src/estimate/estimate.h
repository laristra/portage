/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef ESTIMATOR_H_INC_
#define ESTIMATOR_H_INC_

#include <vector>
#include <tr1/memory>
#include <string>
#include <cmath>
#include <array>

#include "Point.h"
#include "weight.h"
#include "swarm.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::tr1::shared_ptr;
using std::map;

enum EstimatorType {
  Min,
  Max,
  KernelDensity,
  LocalRegression
};

enum WeightCenter {
  Gather,
  Scatter
};

/*!
 * @brief
 */
template<size_t dim>
class Estimator {
 public:
  Estimator(shared_ptr<Swarm<dim>> source, shared_ptr<Swarm<dim>> target,
            EstimatorType estimator, WeightCenter center,
            shared_ptr<vector<Weight::Type>> weights,
            shared_ptr<vector<vector<double>>> smoothing,
            Basis::Type basis)
      : source_(source),
        target_(target),
        weights_(weights),
        basis_(basis),
        estimator_(estimator),
        center_(center)
  {}

  double operator()(int const target_index,
                    std::vector<Weights_t> const & sources_and_weights) const
  {}

  void set_interpolation_variable(std::string const & interp_var_name) {}

 private:
  shared_ptr<Swarm<dim>> source_;
  shared_ptr<Swarm<dim>> target_;
  shared_ptr<vector<Weight::Type>> weights_;
  Basis::Type basis_;
  EstimatorType estimator_;
  WeightCenter center_;
};

}
}

#endif ESTIMATOR_H_INC
