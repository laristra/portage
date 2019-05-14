/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

#include "portage/support/moment_index.h"

// ============================================================================
// Test core routine

template<int D>
void run_moment_index_test() {
  int const NMAX = 100;
  for (int index = 0; index < NMAX; ++index) {
    // Convert index to moment specification
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = Portage::index_to_moment<D>(index);
    // Verify moment specification is internally consistent
    ASSERT_GE(order, 0);
    for (int d = 0; d < D; ++d) {
      ASSERT_GE(exponents[d], 0);
    }
    ASSERT_EQ(order, std::accumulate(exponents.begin(), exponents.end(), 0));
    // Convert moment specification back to index
    auto index2 = Portage::moment_to_index<D>(order, exponents);
    // Ensure the round trip is self-consistent
    ASSERT_EQ(index, index2);
  }
}

// ============================================================================

// 1D
TEST(Moment_Index_test, Test1D) {
  run_moment_index_test<1>();
}

// 2D
TEST(Moment_Index_test, Test2D) {
  run_moment_index_test<2>();
}

// 3D
TEST(Moment_Index_test, Test3D) {
  run_moment_index_test<3>();
}

