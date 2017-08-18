/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#define ESTIMATE_H_INC_

#include <vector>
#include <string>
#include <cassert>

#include "portage/swarm/swarm_state.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;

/**
 * @class Estimate portage/estimate/estimate.h
 * @brief This functional will take the result of 
 * @code Accumulate::operator()@endcode and use it to
 *  estimate derivatives of the source state data.
 */
template<size_t dim, class TargetSwarmState>
class Estimate {
 public:
  /** 
   * @brief Constructor
   * @param source_state the very same source state provided to @code Accumulate::Accumulate()@endcode
   */
  Estimate(SwarmState<dim> const& source_state):
	source_state_(source_state)
  {}

  /** 
   * @brief Set the interpolation variable and derivative
   * @param interp_var_name the name of the variable to estimate
   * @param derivative the derivative specification in the range 0...n-1 where n is the number of 
   *        basis components used in @code Accumulate::operator()@endcode
   */
  void set_interpolation_variable(std::string const & interp_var_name, size_t derivative=0) {
    var_name_ = interp_var_name;
    derivative_ = derivative;
  }

  /**  
   * @brief apply corrected weight function to compute derivatives
   *
   * @param target_index the index of the target swarm particle at which to estimate the derivative
   * @param shape_vec the result of @code Accumulate::operator()@endcode
   * @param source_particles the same source particles provided to @code Accumulate::operator()@endcode
   * @return the derivative of the data in @code source_state_@endcode specified by @code derivative_@endcode
   */
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
