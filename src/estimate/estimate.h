/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef ESTIMATE_H_INC_
#define ESTIMATE_H_INC_

#include <vector>
#include <string>
#include <cassert>

#include "portage/swarm/swarm_state.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;

template<size_t dim>
class Estimate {
 public:
  Estimate(shared_ptr<SwarmState<dim>> source_state):
	source_state_(source_state)
  {}

  double operator()(int const target_index,
                    vector<vector<double>> const & shape_vec,
		    vector<size_t> const& source_particles)
  {
    assert(shape_vec.size() == source_particles.size());
    assert(derivative_ < shape_vec[0].size());

    typename SwarmState<dim>::DblVecPtr source_field_ptr;
    source_state_->get_field(var_name_, source_field_ptr);
    vector<double> &source_field(*source_field_ptr);

    double result=0.;
    for (size_t i=0; i<source_particles.size(); i++) {
      result += source_field[source_particles[i]] * shape_vec[i][derivative_];
    }
    return result;
  }

  void set_interpolation_variable(std::string const & interp_var_name, size_t derivin=0) {
    var_name_ = interp_var_name;
    derivative_ = derivin;
  }

 private:
  shared_ptr<SwarmState<dim>> source_state_;
  std::string var_name_;
  size_t derivative_;
};

}
}

#endif // ESTIMATE_H_INC
