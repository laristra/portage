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

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/weight.h"
#include "portage/support/basis.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::shared_ptr;
using std::map;


/*!
 * @brief
 */
template<size_t dim,
         template <size_t> class /* SourceSwarmState */>
class Estimate {
 public:
  Estimate(SwarmState<dim> const& source_state)
      : source_state_(source_state)
  {}

  void set_variable(std::string const & var_name) {
    var_name_ = var_name;
    source_state_.get_data(PARTICLE, var_name, &source_vals_);
  }

  double operator()(int const target_index,
                    vector<Weights_t> const & sources_and_weights) const
  {

    // contribution of the source cell is its field value weighted by
    // its "weight" 
    double val = 0.0;
    for (auto const& wt : sources_and_weights) {
      int srccell = wt.entityID;
      std::vector<double> pair_weights = wt.weights;
      val += source_vals_[srccell] * pair_weights[0];
    }
    
    return val;
  }

 private:
  SwarmState<dim> const& source_state_;
  std::string var_name_;
  double const * source_vals_;
};

}
}

#endif // ESTIMATE_H_INC
