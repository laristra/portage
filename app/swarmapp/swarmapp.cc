/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <cassert>
#include <cmath>
#ifdef PORTAGE_ENABLE_MPI
  #include <mpi.h>
#endif

#include "portage/support/portage.h"
#include "portage/driver/driver_swarm.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/search/search_simple_points.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "wonton/support/Point.h"
#ifdef HAVE_NANOFLANN
  #include "portage/search/search_kdtree_nanoflann.h"
#endif

using namespace Portage::Meshfree;

//////////////////////////////////////////////////////////////////////
// Helper routines and data structures
//////////////////////////////////////////////////////////////////////

/**
 * @brief Compute the field value for the given point.
 *
 * @tparam dim: spatial dimension.
 * @param field_order: the field order.
 * @param p: the given point.
 * @return the value at p.
 */
template<int dim>
double field_func(int field_order, Wonton::Point<dim> const& p) {
  double value = 0.0;
  switch (field_order) {
    case -1: {
      double rsqr = 0.0;
      for (int i = 0; i < dim; i++) { rsqr += p[i] * p[i]; }
      value = 1.0;
      for (int i = 0; i < dim; i++) { value *= std::sin(0.9 * 2 * M_PI * p[i]); }
      value *= std::exp(-1.5 * std::sqrt(rsqr));
      break;
    }
    case 0: value = 25.3; break;
    case 1: for (int i = 0; i < dim; i++) value += p[i]; break;
    case 2: for (int i = 0; i < dim; i++) value += p[i] * p[i]; break;
    case 3: for (int i = 0; i < dim; i++) value += p[i] * p[i] * p[i]; break;
    default: throw std::runtime_error("invalid field order");
  }

  return value;
}

/**
 * @brief Example properties
 *
 */
class ExampleProperties {
public:
  ExampleProperties(int in_field_order,
                    int in_estimation_order,
                    int in_dimension)
    : field_order(in_field_order),
      estimation_order(in_estimation_order),
      dimension(in_dimension)
  {}

  int estimation_order = 0;   // estimation order in example
  int field_order = 0;        // order of the field to map
  int dimension = 0;          // spatial dimension of problem
};

/**
 * @brief Example properties list.
 *
 * Use this to add new problems. If needed, the related struct
 * can be extented to contain more information.
 */
static ExampleProperties const examples[20] = {
  { 0, 0, 2}, // 2D const func, unitary basis
  { 1, 0, 2}, // 2D linear func, unitary basis
  { 1, 1, 2}, // 2D linear func, linear basis
  { 2, 0, 2}, // 2D quadratic func, unitary basis
  { 2, 1, 2}, // 2D quadratic func, linear basis
  { 2, 2, 2}, // 2D quadratic func, quadratic basis
  { 3, 0, 2}, // 2D cubic func, unitary basis
  { 3, 1, 2}, // 2D cubic func, linear basis
  { 3, 2, 2}, // 2D cubic func, quadratic basis
  {-1, 2, 2}, // 2D exp(-1.5*r)*sin(0.9*2*pi*x)*sin(0.9*2*pi*y)
  { 0, 0, 3}, // 3D const func, unitary basis
  { 1, 0, 3}, // 3D linear func, unitary basis
  { 1, 1, 3}, // 3D linear func, linear basis
  { 2, 0, 3}, // 3D quadratic func, unitary basis
  { 2, 1, 3}, // 3D quadratic func, linear basis
  { 2, 2, 3}, // 3D quadratic func, quadratic basis
  { 3, 0, 3}, // 3D cubic func, unitary basis
  { 3, 1, 3}, // 3D cubic func, linear basis
  { 3, 2, 3}, // 3D cubic func, quadratic basis
  {-1, 2, 3}  // 3D exp(-1.5*r)*sin(0.9*2*pi*x)*sin(0.9*2*pi*y)*sin(0.9*2*pi*z)
};


/**
 * @brief Print command usage and options.
 *
 */
void print_usage() {
  std::cout << "Usage: swarmapp example-number nsourcepts ntargetpts distribution seed hfactor center" << std::endl;
  std::cout << "- example-number must be between 0-9 for 2D and 10-19 for 3D" << std::endl;
  std::cout << "- nsourcepts is the number of source points per axis" << std::endl;
  std::cout << "- ntargetpts is the number of target points per axis" << std::endl;
  std::cout << "- distribution = 0 for random, 1 for regular grid, 2 for perturbed grid" << std::endl;
  std::cout << "- seed sets the random number seed to get reproducible results" << std::endl;
  std::cout << "- hfactor scales the smoothing length in multiples of the average inter-particle spacing" << std::endl;
  std::cout << "- center governs where weights are centered: 0 for target, 1 for source " << std::endl;

  std::cout << "List of example numbers:" << std::endl;
  for (int i = 0; i < 20; ++i) {
    std::printf("  %d: order %d remap of order %d function\n",
                i, examples[i].estimation_order, examples[i].field_order);
  }
}

/**
 * @brief Run the particle remap.
 *
 * @tparam dim: spatial dimension.
 * @param example_num: the example number.
 * @param n_source: number of source points per axis.
 * @param n_target: number of target points per axis.
 * @param distribution: 0 for random, 1 for regular grid, 2 for perturbed grid.
 * @param seed: random number seed to get reproducible results.
 * @param hfactor: smoothing lengths scaling factor.
 * @param center: governs where weights are centered: 0 for target, 1 for source.
 */
template<int dim>
void run(int example_num, int n_source, int n_target,
         int distribution, unsigned seed, double hfactor, int center);

/**
 * @brief Run the 2D particle remap.
 *
 * @param example_num: the example number.
 * @param n_source: number of source points per axis.
 * @param n_target: number of target points per axis.
 * @param distribution: 0 for random, 1 for regular grid, 2 for perturbed grid.
 * @param seed: random number seed to get reproducible results.
 * @param hfactor: smoothing lengths scaling factor.
 * @param center: governs where weights are centered: 0 for target, 1 for source.
 *
 * NOTE: gather form tests.
 */
template<>
void run<2>(int example_num, int n_source, int n_target,
            int distribution, unsigned seed, double hfactor, int center) {

  std::cout << "starting swarm app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  auto const& example = examples[example_num];
  std::cout << "estimate order: " << example.estimation_order << std::endl;

  // Regularly ordered input swarm; randomly ordered output swarm
  Swarm<2> source_swarm(n_source * n_source, distribution, seed,
                        -1.1, 1.1, -1.1, 1.1);
  Swarm<2> target_swarm(n_target * n_target, distribution, seed,
                        -1.0, 1.0, -1.0, 1.0);

  SwarmState<2> source_state(source_swarm);
  SwarmState<2> target_state(target_swarm);

  int const num_source_particles = source_swarm.num_particles();
  int const num_target_particles = target_swarm.num_particles();

  Portage::vector<double> source_data(num_source_particles, 0.);
  Portage::vector<double> target_data(num_target_particles, 0.);
  std::vector<std::string> remap_fields;

  for (int i = 0; i < num_source_particles; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data[i] = field_func<2>(example.field_order, p);
  }
  
  source_state.add_field("remapdata", source_data);
  target_state.add_field("remapdata", target_data);
  remap_fields.emplace_back("remapdata");

  // Smoothing lengths (or in other words "size" of the support function) in each dimension
  double h = 0.;
  int nsmooth = 0;
  WeightCenter center_type;

  if (center == 0) {
    center_type = Gather;
    h = 2.0 * hfactor / n_target;
    nsmooth = num_target_particles;
  } else {
    center_type = Scatter;
    h = 2.2 * hfactor / n_source;
    nsmooth = num_source_particles;
  }

  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(2, h));
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths(nsmooth, default_lengths);

#ifdef HAVE_NANOFLANN  // Search by kdtree
  using Remapper = SwarmDriver<Portage::Search_KDTree_Nanoflann,
                               Accumulate, Estimate, 2,
                               Swarm<2>, SwarmState<2>>;
#else
  using Remapper = SwarmDriver<Portage::SearchPointsByCells,
                               Accumulate, Estimate, 2,
                               Swarm<2>, SwarmState<2>>;
#endif

  Remapper remapper(source_swarm, source_state,
                    target_swarm, target_state,
                    smoothing_lengths, Weight::B4,
                    Weight::ELLIPTIC, center_type);

  switch (example.estimation_order) {
    case 0: remapper.set_remap_var_names(remap_fields, remap_fields,
                                         LocalRegression,
                                         Basis::Unitary); break;
    case 1: remapper.set_remap_var_names(remap_fields, remap_fields,
                                        LocalRegression,
                                        Basis::Linear); break;
    case 2: remapper.set_remap_var_names(remap_fields, remap_fields,
                                         LocalRegression,
                                         Basis::Quadratic); break;
    default: break;
  }

#ifdef PORTAGE_ENABLE_MPI
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  remapper.run(&mpiexecutor, true);
#else
  remapper.run();
#endif

  auto& target_field = target_state.get_field("remapdata");
  std::vector<double> expected_value(num_target_particles, 0.0);
  std::vector<double> total_error(3, 0.);

  for (int i = 0; i < num_target_particles; ++i) {
    auto p = target_swarm.get_particle_coordinates(i);
    expected_value[i] = field_func<2>(example.field_order, p);
    double error = std::abs(expected_value[i] - target_field[i]);
    total_error[0] = std::max(error, total_error[0]);
    total_error[1] += error;
    total_error[2] += error * error;
  }

  total_error[1] /= num_source_particles;
  total_error[2] /= num_source_particles;
  total_error[2] = std::sqrt(total_error[2]);

  std::cout << std::endl;
  std::cout << std::scientific;
  std::cout.precision(17);
  std::cout << "weight center = " << center << std::endl;
  std::cout << "distribution = " << distribution << std::endl;
  std::cout << "dimension = 2" << std::endl;
  std::cout << "seed = " << seed << std::endl;
  std::cout << std::endl;
  std::cout << "smoothing length = " << h << std::endl;
  std::cout << "Linf NORM OF ERROR = " << total_error[0] << std::endl;
  std::cout << "L1 NORM OF ERROR = " << total_error[1] << std::endl;
  std::cout << "L2 NORM OF ERROR = " << total_error[2] << std::endl;
  
  // Create the input file for comparison.
  std::string source_values_csv = "infield" + std::to_string(example_num) + ".csv";
  std::string target_values_csv = "outfield-"+ std::to_string(example_num) +"-"+ std::to_string(center) + ".csv";

  std::ofstream file[2];
  file[0].open(source_values_csv);
  file[0] << std::scientific;
  file[0].precision(17);
  file[0] << " X, Y, Value\n";

  for (int i = 0; i < num_source_particles; ++i) {
    auto const p = source_swarm.get_particle_coordinates(i);
    file[0] << p[0] << ", " << p[1] << ", " << source_data[i] << std::endl;
  }

  // Create the output file for comparison.
  file[1].open(target_values_csv);
  file[1] << std::scientific;
  file[1].precision(17);
  file[1] << " X, Y, Value\n";

  for (int i = 0; i < num_target_particles; ++i) {
    auto const p = target_swarm.get_particle_coordinates(i);
    std::cout << p << std::endl;
    file[1] << p[0] << ", " << p[1] << ", " << target_field[i] << std::endl;
  }

  std::cout << "finishing swarm app..." << std::endl;
}

/**
 * @brief Run the 3D particle remap.
 *
 * @param example_num: the example number.
 * @param n_source: number of source points per axis.
 * @param n_target: number of target points per axis.
 * @param distribution: 0 for random, 1 for regular grid, 2 for perturbed grid.
 * @param seed: random number seed to get reproducible results.
 * @param hfactor: smoothing lengths scaling factor.
 * @param center: governs where weights are centered: 0 for target, 1 for source.
 *
 * NOTE: gather form tests.
 */
template<>
void run<3>(int example_num, int n_source, int n_target,
            int distribution, unsigned seed, double hfactor, int center) {

  std::cout << "starting swarm app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  auto const& example = examples[example_num];

  // Regularly ordered input swarm; randomly ordered output swarm
  Swarm<3> source_swarm(n_source * n_source * n_source, distribution, seed,
                        -1.1, 1.1, -1.1, 1.1, -1.1, 1.1);
  Swarm<3> target_swarm(n_target * n_target * n_target, distribution, seed,
                        -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

  SwarmState<3> source_state(source_swarm);
  SwarmState<3> target_state(target_swarm);

  int const num_source_particles = source_swarm.num_particles();
  int const num_target_particles = target_swarm.num_particles();

  Portage::vector<double> source_data(num_source_particles, 0.);
  Portage::vector<double> target_data(num_target_particles, 0.);
  std::vector<std::string> remap_fields;

  for (int i = 0; i < num_source_particles; ++i) {
    auto coord = source_swarm.get_particle_coordinates(i);
    source_data[i] = field_func<3>(example.field_order, coord);
  }

  source_state.add_field("remapdata", source_data);
  target_state.add_field("remapdata", target_data);
  remap_fields.emplace_back("remapdata");

  // Smoothing lengths (or in other words "size" of the support function) in each dimension
  double h = 0.;
  int nsmooth = 0;
  WeightCenter center_type {};

  if (center == 0) {
    center_type = Gather;
    h = 2.0 * hfactor / n_target;
    nsmooth = num_target_particles;
  } else {
    center_type = Scatter;
    h = 2.2 * hfactor / n_source;
    nsmooth = num_source_particles;
  }

  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(3, h));
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths(nsmooth, default_lengths);

#ifdef HAVE_NANOFLANN  // Search by kdtree
  using Remapper = SwarmDriver<Portage::Search_KDTree_Nanoflann,
                               Accumulate, Estimate, 3,
                               Swarm<3>, SwarmState<3>>;
#else
  using Remapper = SwarmDriver<Portage::SearchPointsByCells,
                               Accumulate, Estimate, 3,
                               Swarm<3>, SwarmState<3>>;
#endif

  Remapper remapper(source_swarm, source_state,
                    target_swarm, target_state,
                    smoothing_lengths, Weight::B4,
                    Weight::ELLIPTIC, center_type);


  switch (example.estimation_order) {
    case 0: remapper.set_remap_var_names(remap_fields, remap_fields,
                                         LocalRegression,
                                         Basis::Unitary); break;
    case 1: remapper.set_remap_var_names(remap_fields, remap_fields,
                                         LocalRegression,
                                         Basis::Linear); break;
    case 2: remapper.set_remap_var_names(remap_fields, remap_fields,
                                         LocalRegression,
                                         Basis::Quadratic); break;
    default: break;
  }

#ifdef PORTAGE_ENABLE_MPI
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  remapper.run(&mpiexecutor, true);
#else
  remapper.run();
#endif

  auto& target_field = target_state.get_field("remapdata");
  std::vector<double> expected_value(num_target_particles, 0.0);
  std::vector<double> total_error(3, 0.);

  for (int i = 0; i < num_target_particles; ++i) {
    auto p = target_swarm.get_particle_coordinates(i);
    expected_value[i] = field_func<3>(example.field_order, p);
    double error = std::abs(expected_value[i] - target_field[i]);
    total_error[0] = std::max(error, total_error[0]);
    total_error[1] += error;
    total_error[2] += error * error;
  }

  total_error[1] /= num_source_particles;
  total_error[2] /= num_source_particles;
  total_error[2] = std::sqrt(total_error[2]);

  std::cout << std::endl;
  std::cout << std::scientific;
  std::cout.precision(17);
  std::cout << "weight center = " << center << std::endl;
  std::cout << "distribution = " << distribution << std::endl;
  std::cout << "dimension = 3" << std::endl;
  std::cout << "seed = " << seed << std::endl;
  std::cout << std::endl;
  std::cout << "smoothing length = " << h << std::endl;
  std::cout << "Linf NORM OF ERROR = " << total_error[0] << std::endl;
  std::cout << "L1 NORM OF ERROR = " << total_error[1] << std::endl;
  std::cout << "L2 NORM OF ERROR = " << total_error[2] << std::endl;

  // Create the input file for comparison.
  std::string source_values_csv = "infield" + std::to_string(example_num) + ".csv";
  std::string target_values_csv = "outfield-"+ std::to_string(example_num) +"-"+ std::to_string(center) + ".csv";

  std::ofstream file[2];
  file[0].open(source_values_csv);
  file[0] << std::scientific;
  file[0].precision(17);
  file[0] << " X, Y, Value\n";

  for (int i = 0; i < num_source_particles; ++i) {
    auto const p = source_swarm.get_particle_coordinates(i);
    file[0] << p[0] << ", " << p[1] << ", " << p[2] << ", "<< source_data[i] << std::endl;
  }

  // Create the output file for comparison.
  file[1].open(target_values_csv);
  file[1] << std::scientific;
  file[1].precision(17);
  file[1] << " X, Y, Z, Value\n";

  for (int i = 0; i < num_target_particles; ++i) {
    auto const p = target_swarm.get_particle_coordinates(i);
    std::cout << p << std::endl;
    file[1] << p[0] << ", " << p[1] << ", " << p[2] << ", " << target_field[i] << std::endl;
  }

  std::cout << "finishing swarm app..." << std::endl;
}

/**
 * @brief Main method.
 *
 * @param argc: number of arguments
 * @param argv: array of arguments.
 * @return status code: 0 if ok
 */
int main(int argc, char** argv) {

#ifdef PORTAGE_ENABLE_MPI
  int nb_ranks = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_ranks);
#endif

  if (argc < 7) {
    print_usage();
#ifdef PORTAGE_ENABLE_MPI
    MPI_Finalize();
#endif
    return EXIT_FAILURE;
  }

  int example    = std::atoi(argv[1]);
  int n_source   = std::atoi(argv[2]);
  int n_target   = std::atoi(argv[3]);
  int distrib    = std::atoi(argv[4]);
  unsigned seed  = std::atoi(argv[5]);
  int dimension  = (example < 10 ? 2 : 3);
  double hfactor = std::stod(argv[6]);
  int center     = std::atoi(argv[7]);

  // check inputs
  assert(n_source > 0);
  assert(n_target > 0);
  assert(0 <= example and example < 20);
  assert(0 <= distrib and distrib <= 2);
  assert(center == 0 or center == 1);

  switch (dimension) {
    case 2: run<2>(example, n_source, n_target, distrib, seed, hfactor, center); break;
    case 3: run<3>(example, n_source, n_target, distrib, seed, hfactor, center); break;
    default: throw std::runtime_error("invalid dimension");
  }

#ifdef PORTAGE_ENABLE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
