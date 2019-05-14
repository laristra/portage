/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef PORTAGE_SUPPORT_MOMENT_INDEX_
#define PORTAGE_SUPPORT_MOMENT_INDEX_

#include <array>
#include <cassert>
#include <cmath>
#include <utility>

//! This file contains convenience functions to help extract moments from the
//! moments list by converting between a moment specification and a moments
//! list index.

namespace Portage {

  // ==========================================================================

  //! Convert from moment specification to index
  //! @param[in]  order  Order of the moment
  //! @param[in]  exponents  Powers for the coordinate axes
  template<int D>
  int moment_to_index(int const order, std::array<int,D> const exponents) {
    static_assert(false, "\"moment_to_index\" only works in 1, 2, or 3D");
  }

  template<>
  int moment_to_index<3>(int const order, std::array<int,3> const exponents) {
    int constexpr D = 3;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    int index = 0;
    index += (order+2) * (order+1) * order / 6;
    int order1 = order - exponents[0];
    index += (order1+1) * order1 / 2;
    index += order1 - exponents[1];
    return(index);
  }

  template<>
  int moment_to_index<2>(int const order, std::array<int,2> const exponents) {
    int constexpr D = 2;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    int index = 0;
    index += (order+1) * order / 2;
    index += order - exponents[0];
    return(index);
  }

  template<>
  int moment_to_index<1>(int const order, std::array<int,1> const exponents) {
    int constexpr D = 1;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    int index = 0;
    index += order;
    return(index);
  }

  // ==========================================================================

  //! Convert from index to moment specification
  //! @param[in]  index  Index in moment list
  template<int D>
  std::pair<int,std::array<int,D>> index_to_moment(int const index) {
    static_assert(false, "\"index_to_moment\" only works in 1, 2, or 3D");
  }

  template<>
  std::pair<int,std::array<int,3>> index_to_moment<3>(int const index) {
    int constexpr D = 3;
    int n = index;
    // Declare the output variables
    int order;
    std::array<int,D> exponents;
    // Find order
    //    While this could be done analytically, the expression is ugly and
    // involves cube roots and complex numbers.  With only a small number of
    // orders, this should be plenty fast enough (and possibly faster than the
    // ugly expression).
    order = 0;
    while (true) {
      if (n < (order+2)*(order+1)*order/6)
        break;
      order++;
    }
    order--;
    // Skip the lower-order moments
    n -= (order+2)*(order+1)*order/6;
    // Find the first exponent
    exponents[0] = order - (int) floor((sqrt((double)8*n+1) - 1) / 2);
    auto order1 = order - exponents[0];
    // Skip moments with wrong first exponent
    n -= (order1+1) * order1 / 2;
    // Find the second exponent
    exponents[1] = order1 - n;
    // Find the third exponent
    exponents[2] = n;
    // Verify results
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    // Return
    return(std::move(std::make_pair(order,exponents)));
  }

  template<>
  std::pair<int,std::array<int,2>> index_to_moment<2>(int const index) {
    int constexpr D = 2;
    int n = index;
    // Declare the output variables
    int order;
    std::array<int,D> exponents;
    // Find order
    order = (int) floor((sqrt((double)8*n+1) - 1) / 2);
    // Skip the lower-order moments
    n -= (order+1) * order / 2;
    // Find the first exponent
    exponents[0] = order - n;
    // Find the second exponent
    exponents[1] = n;
    // Verify results
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    // Return
    return(std::move(std::make_pair(order,exponents)));
  }

  template<>
  std::pair<int,std::array<int,1>> index_to_moment<1>(int const index) {
    int constexpr D = 1;
    // Declare the output variables
    int order;
    std::array<int,D> exponents;
    // Find order
    order = index;
    // Find the first exponent
    exponents[0] = index;
    // Verify results
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    // Return
    return(std::move(std::make_pair(order,exponents)));
  }

}  // namespace Portage

#endif  // PORTAGE_SUPPORT_MOMENT_INDEX_
