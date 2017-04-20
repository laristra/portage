/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef ESTIMATE_H_INC_
#define ESTIMATE_H_INC_

#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <array>
#include <cassert>

#include "portage/support/Point.h"
#include "weight.h"
#include "basis.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::shared_ptr;
using std::map;

enum EstimateType {
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
class Estimate {
 public:
  Estimate(shared_ptr<Swarm<dim>> source, shared_ptr<Swarm<dim>> target,
	   shared_ptr<SwarmState<dim>> source_state, shared_ptr<SwarmState<dim>> target_state,
            EstimateType estimate, WeightCenter center,
            shared_ptr<vector<Weight::Kernel>> kernels,
            shared_ptr<vector<Weight::Geometry>> geometries,
            shared_ptr<vector<vector<double>>> smoothing,
            Basis::Type basis)
      : source_(source),
        target_(target),
	source_state_(source_state),
	target_state_(target_state),
        estimate_(estimate),
        center_(center)
        kernels_(kernels),
        geometries_(geometries),
        smoothing_(smoothing),
        basis_(basis)
  {}

  double operator()(int const target_index,
                    vector<vector<vector<double>>> const & moment_matrix_bit) const
  {}

  void set_interpolation_variable(std::string const & interp_var_name) {}

 private:
  shared_ptr<Swarm<dim>> source_;
  shared_ptr<Swarm<dim>> target_;
  shared_ptr<SwarmState<dim>> source_state_;
  shared_ptr<SwarmState<dim>> target_state_;
  EstimateType estimate_;
  WeightCenter center_;
  shared_ptr<vector<Weight::Kernel>> kernels_;
  shared_ptr<vector<Weight::Geometry>> geometries_;
  shared_ptr<vector<vector<double>>> smoothing_;
  Basis::Type basis_;
};

}
}

#endif // ESTIMATE_H_INC
