/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <string>
#include <limits>

#include "gtest/gtest.h"

#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

// portage includes
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/operator.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

namespace {

double TOL = 1e-12;

/**
 * @brief Compute field value based on point coordinates.
 *
 * @tparam dim: spatial dimension.
 * @param p: current point coordinates.
 * @return related field value.
 */
template<int dim>
double exact_value(Wonton::Point<dim> const& p) {
  switch (dim) {
    case 2: return (p[0] < 0. and p[1] < 0.) or (p[0] > 0. and p[1] > 0.) ? 2.0 : 1.0;
    case 3: return (p[0] < 0. and p[1] < 0. and p[2] < 0.)
                or (p[0] > 0. and p[1] > 0. and p[2] > 0.) ? 2.0 : 1.0;
    default: return 0.0;
  }
}


TEST(PartByParticle, 2D) {

  using namespace Portage::swarm;

  int const ncells = 4;

  Wonton::Simple_Mesh source_mesh(-1., -1., 1., 1., ncells, ncells);
  Wonton::Simple_Mesh target_mesh(-1., -1., 1., 1., 2 * ncells + 2, 2 * ncells + 2);
  Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(source_mesh);
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(target_mesh);

  int const nb_source = source_mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  int const nb_target = target_mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  assert(nb_source == ncells * ncells);

  double source_values[nb_source];
  double target_values[nb_target];

  // set source field values
  for (int i = 0; i < nb_source; i++) {
    Wonton::Point<2> p;
    source_mesh_wrapper.cell_centroid(i, &p);
    source_values[i] = exact_value<2>(p);
  }

  Wonton::Simple_State source_state(std::make_shared<Wonton::Simple_Mesh>(source_mesh));
  Wonton::Simple_State target_state(std::make_shared<Wonton::Simple_Mesh>(target_mesh));
  Wonton::Simple_State_Wrapper source_state_wrapper(source_state);
  Wonton::Simple_State_Wrapper target_state_wrapper(target_state);
  source_state.add("indicate", Wonton::CELL, source_values);
  target_state.add("indicate", Wonton::CELL, target_values);

  using Remapper = Portage::MSM_Driver<Portage::SearchPointsByCells,
                                       Accumulate, Estimate, 2,
                                       Wonton::Simple_Mesh_Wrapper,
                                       Wonton::Simple_State_Wrapper>;

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper,
                    1.0, 1.0, Weight::FACETED, Weight::POLYRAMP,
                    Scatter, "indicate", 0.25);

  remapper.set_remap_var_names({"indicate"}, {"indicate"});
  remapper.run();

  double* remapped;
  target_state_wrapper.mesh_get_data(Wonton::CELL, "indicate", &remapped);

  for (int i = 0; i < nb_target; i++) {
    Wonton::Point<2> p;
    target_mesh_wrapper.cell_centroid(i, &p);
    double expected = exact_value<2>(p);
    ASSERT_NEAR(expected, remapped[i], TOL);
  }
}


TEST(PartByParticle, 3D) {

  using namespace Portage::swarm;

  int const ncells = 4;

  Wonton::Simple_Mesh source_mesh(-1., -1., -1., 1., 1., 1., ncells, ncells, ncells);
  Wonton::Simple_Mesh target_mesh(-1., -1., -1., 1., 1., 1., 2 * ncells + 2, 2 * ncells + 2, 2 * ncells + 2);
  Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(source_mesh);
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(target_mesh);

  int const nb_source = source_mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  int const nb_target = target_mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  assert(nb_source == ncells * ncells * ncells);

  double source_values[nb_source];
  double target_values[nb_target];

  // set source field values
  for (int i = 0; i < nb_source; i++) {
    Wonton::Point<3> p;
    source_mesh_wrapper.cell_centroid(i, &p);
    source_values[i] = exact_value<3>(p);
  }

  Wonton::Simple_State source_state(std::make_shared<Wonton::Simple_Mesh>(source_mesh));
  Wonton::Simple_State target_state(std::make_shared<Wonton::Simple_Mesh>(target_mesh));
  Wonton::Simple_State_Wrapper source_state_wrapper(source_state);
  Wonton::Simple_State_Wrapper target_state_wrapper(target_state);
  source_state.add("indicate", Wonton::CELL, source_values);
  target_state.add("indicate", Wonton::CELL, target_values);

  using Remapper = Portage::MSM_Driver<Portage::SearchPointsByCells,
                                       Accumulate, Estimate, 3,
                                       Wonton::Simple_Mesh_Wrapper,
                                       Wonton::Simple_State_Wrapper>;

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper,
                    1.0, 1.0, Weight::FACETED, Weight::POLYRAMP,
                    Scatter, "indicate", 0.25);

  remapper.set_remap_var_names({"indicate"}, {"indicate"});
  remapper.run();

  double* remapped;
  target_state_wrapper.mesh_get_data(Wonton::CELL, "indicate", &remapped);

  for (int i = 0; i < nb_target; i++) {
    Wonton::Point<3> p;
    target_mesh_wrapper.cell_centroid(i, &p);
    double expected = exact_value<3>(p);
    ASSERT_NEAR(expected, remapped[i], TOL);
  }
}

}  // namespace
