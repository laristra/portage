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

template<size_t dim, class TargetSwarmState>
class Estimate {
 public:
  Estimate(SwarmState<dim> const& source_state):
      source_state_(source_state)
  {}

  double operator()(int const target_index,
                    vector<Weights_t> const &sources_and_mults) const
  {    
    int nsrc = sources_and_mults.size();
    assert(nsrc > 0);

    assert(derivative_ < sources_and_mults[0].weights.size());

    double result=0.;
    for (size_t i=0; i<nsrc; i++) {
      Weights_t const& wt = sources_and_mults[i];
      int p = wt.entityID;
      vector<double> const& shape_vec = wt.weights;
      result += source_vals_[p] * shape_vec[derivative_];
    }
    return result;
  }

  void set_variable(std::string const & var_name, size_t derivin=0) {
    var_name_ = var_name;
    derivative_ = derivin;
    typename SwarmState<dim>::DblVecPtr source_field_ptr;
    source_state_.get_field(var_name_, source_field_ptr);
    source_vals_ = &((*source_field_ptr)[0]);
  }

 private:
  SwarmState<dim> const& source_state_;
  std::string var_name_;
  size_t derivative_;
  double const *source_vals_;
};

}
}

#endif // ESTIMATE_H_INC
