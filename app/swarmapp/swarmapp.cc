/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <cassert>
#include <cmath>

#ifdef PORTAGE_ENABLE_MPI
#include <mpi.h>
#else
#define PORTAGE_SERIAL_ONLY
#endif

#include "portage/support/portage.h"
#include "portage/driver/driver_swarm.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/search/search_simple_points.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#ifdef HAVE_NANOFLANN
#include "portage/search/search_kdtree_nanoflann.h"
#include "wonton/support/Point.h"
#endif

using std::shared_ptr;
using std::make_shared;
using Portage::Meshfree::Swarm;
using Portage::Meshfree::SwarmState;
using Portage::Meshfree::SwarmDriver;
using Portage::Meshfree::Accumulate;
using Portage::Meshfree::Estimate;
using Portage::SearchSimplePoints;
#ifdef HAVE_NANOFLANN
using Portage::Search_KDTree_Nanoflann;
#endif
using Portage::SearchPointsByCells;


//////////////////////////////////////////////////////////////////////
// Helper routines and data structures
//////////////////////////////////////////////////////////////////////

template<unsigned int D>
double field_func(int field_order, Wonton::Point<D> coord) {
  double value = 0.0;
  switch (field_order) {
    case -1: {
      double rsqr = 0.0;
      for (int i = 0; i < D; i++)
        rsqr += coord[i]*coord[i];
      value = 1.0;
      for (int i = 0; i < D; i++)
        value *= sin(0.9*2*M_PI*coord[i]);
      value *= exp(-1.5*sqrt(rsqr));
      break;
    }
    case 0:
      value = 25.3;
      break;
    case 1:
      for (int i = 0; i < D; i++) value += coord[i];
      break;
    case 2:
      for (int i = 0; i < D; i++) value += coord[i]*coord[i];
      break;
    case 3:
      for (int i = 0; i < D; i++) value += coord[i]*coord[i]*coord[i];
      break;
    default:
      throw std::runtime_error("Unknown field_order!");
  }

  return value;
}

struct example_properties {
  example_properties(const int field_order0, const int estimation_order0, const int dimension0)
    : field_order(field_order0), estimation_order(estimation_order0), dimension(dimension0)
  { }
  int estimation_order;   // estimation order in example
  int field_order;        // order of the field to map
  int dimension;          // spatial dimension of problem
};

// Use this to add new problems.  If needed, we can extend the
// example_properties struct to contain more information.
std::vector<example_properties> setup_examples() {
  std::vector<example_properties> examples;

  // 2D

  // const func, unitary basis
  examples.emplace_back(0, 0, 2);

  // linear func, unitary basis
  examples.emplace_back(1, 0, 2);

  // linear func, linear basis
  examples.emplace_back(1, 1, 2);

  // quadratic func, unitary basis
  examples.emplace_back(2, 0, 2);

  // quadratic func, linear basis
  examples.emplace_back(2, 1, 2);

  // quadratic func, quadratic basis
  examples.emplace_back(2, 2, 2);

  // cubic func, unitary basis
  examples.emplace_back(3, 0, 2);

  // cubic func, linear basis
  examples.emplace_back(3, 1, 2);

  // cubic func, quadratic basis
  examples.emplace_back(3, 2, 2);

  // exp(-1.5*r)*sin(0.9*2*pi*x)*sin(0.9*2*pi*y)
  examples.emplace_back(-1, 2, 2);

  // 3D

  // const func, unitary basis
  examples.emplace_back(0, 0, 3);

  // linear func, unitary basis
  examples.emplace_back(1, 0, 3);

  // linear func, linear basis
  examples.emplace_back(1, 1, 3);

  // quadratic func, unitary basis
  examples.emplace_back(2, 0, 3);

  // quadratic func, linear basis
  examples.emplace_back(2, 1, 3);

  // quadratic func, quadratic basis
  examples.emplace_back(2, 2, 3);

  // cubic func, unitary basis
  examples.emplace_back(3, 0, 3);

  // cubic func, linear basis
  examples.emplace_back(3, 1, 3);

  // cubic func, quadratic basis
  examples.emplace_back(3, 2, 3);

  // exp(-1.5*r)*sin(0.9*2*pi*x)*sin(0.9*2*pi*y)*sin(0.9*2*pi*z)
  examples.emplace_back(-1, 2, 3);

  return examples;
}

void usage() {
  auto examples = setup_examples();
  std::cout << "Usage: swarmapp example-number nsourcepts ntargetpts distribution seed hfactor center"
            << std::endl;
  std::cout << "example-number must be between 0-9 for 2D and 10-19 for 3D"
            << std::endl;
  std::cout << "nsourcepts is the number of source points**(1/dimension)"
            << std::endl;
  std::cout << "ntargetpts is the number of target points**(1/dimension)"
            << std::endl;
  std::cout << "distribution = 0 for random, 1 for regular grid, 2 for perturbed regular grid"
            << std::endl;
  std::cout << "seed sets the random number seed to get reproducible results" 
            << std::endl;  
  std::cout << "hfactor scales the smoothing length in multiples of the average inter-particle spacing" 
            << std::endl;
  std::cout << "center governs where weights are centered: 0 for target, 1 for source " 
            << std::endl;

  std::cout << "List of example numbers:" << std::endl;
  int i = 0;
  for (const auto &example : examples) {
    std::printf("  %d: order %d remap of order %d function\n",
                i, example.estimation_order,
                example.field_order);
    ++i;
  }
}

// gather form tests
void run_test_2d(int example_num, int n_source, int n_target, 
                 int distribution, int dimension, unsigned int seed, 
                 double hfactor, int center) 
{
  std::cout << "starting swarm app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  example_properties example = setup_examples()[example_num];

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

  for (int p = 0; p < num_source_particles; ++p) {
    auto coord = source_swarm.get_particle_coordinates(p);
    source_data[p] = field_func<2>(example.field_order, coord);
  }
  
  source_state.add_field("remapdata", source_data);
  target_state.add_field("remapdata", target_data);
  remap_fields.emplace_back("remapdata");

  // Smoothing lengths (or in other words "size" of the support
  // function) in each dimension
  double h = 0.;
  int nsmooth = 0;
  Portage::Meshfree::WeightCenter center_type;

  if (center == 0) {
    center_type = Portage::Meshfree::Gather;
    h = 2.0 * hfactor / n_target;
    nsmooth = num_target_particles;
  } else if (center == 1) {
    center_type = Portage::Meshfree::Scatter;
    h = 2.2 * hfactor / n_source;
    nsmooth = num_source_particles;
  }

  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(2, h));
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths(nsmooth, default_lengths);

#ifdef HAVE_NANOFLANN  // Search by kdtree
  SwarmDriver<
    Search_KDTree_Nanoflann,
    Accumulate,
    Estimate,
    2,
    Swarm<2>,
    SwarmState<2>>
      d(source_swarm, source_state, target_swarm, target_state,
        smoothing_lengths, Portage::Meshfree::Weight::B4, 
        Portage::Meshfree::Weight::ELLIPTIC, center_type);
#else
  SwarmDriver<
    SearchPointsByCells,
    Accumulate,
    Estimate,
    2,
    Swarm<2>,
    SwarmState<2>>
      d(source_swarm, source_state, target_swarm, target_state,
        smoothing_lengths, Portage::Meshfree::Weight::B4, 
        Portage::Meshfree::Weight::ELLIPTIC, center_type);
#endif

  if (example.estimation_order == 0)
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Unitary);
  else if (example.estimation_order == 1)
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Linear);
  else if (example.estimation_order == 2)
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Quadratic);

#ifdef PORTAGE_ENABLE_MPI
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  d.run(&mpiexecutor, true);
#else
  d.run();
#endif

  std::vector<double> expected_value(num_target_particles, 0.0);

  std::vector<double> toterr(3,0.);
  for (int p(0); p < num_target_particles; ++p) {
    auto coord = target_swarm.get_particle_coordinates(p);

    expected_value[p] = field_func<2>(example.field_order, coord);

    double error = fabs(expected_value[p] - target_data[p]);
    toterr[0] = fmax(error, toterr[0]);
    toterr[1] = toterr[1]+error;
    toterr[2] = toterr[2]+error*error;
  }
  toterr[1] /= num_source_particles;
  toterr[2] /= num_source_particles;
  toterr[2] = sqrt(toterr[2]);

  std::cout << std::endl;
  std::cout << std::scientific;
  std::cout.precision(17);
  std::cout << "weight center = " << center << std::endl;
  std::cout << "distribution = " << distribution << std::endl;
  std::cout << "dimension = " << dimension << std::endl;
  std::cout << "seed = " << seed << std::endl;
  std::cout << std::endl;
  std::cout << "smoothing length = " << h << std::endl;
  std::cout << "Linf NORM OF ERROR = " << toterr[0] << std::endl;
  std::cout << "L1 NORM OF ERROR = " << toterr[1] << std::endl;
  std::cout << "L2 NORM OF ERROR = " << toterr[2] << std::endl;
  
  // Create the input file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream finp_csv("infield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".csv");
  finp_csv << std::scientific;
  finp_csv.precision(17);
  finp_csv << " X, Y, Value\n";
  for (int p(0); p < num_source_particles; ++p) {
    Wonton::Point<2> coord = source_swarm.get_particle_coordinates(p);
    finp_csv << coord[0] << ", " << coord[1] << ", " << source_data[p] << std::endl;
  }

  // Create the output file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream fout_csv("outfield-"
                     + std::to_string(static_cast<long long>(example_num)) + "-"
                     + std::to_string(static_cast<long long>(center))
                     + ".csv");
  fout_csv << std::scientific;
  fout_csv.precision(17);
  fout_csv << " X, Y, Value\n";
  for (int p(0); p < num_target_particles; ++p) {
    Wonton::Point<2> coord = target_swarm.get_particle_coordinates(p);
    fout_csv << coord[0] << ", " << coord[1] << ", " << target_data[p] << std::endl;
  }
  std::cout << "finishing swarm app..." << std::endl;
}

// gather form tests
void run_test_3d(int example_num, int n_source, int n_target, 
                 int distribution, int dimension, unsigned int seed, 
                 double hfactor, int center) 
{
  std::cout << "starting swarm app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  example_properties example = setup_examples()[example_num];

  // Regularly ordered input swarm; randomly ordered output swarm
  auto source_swarm = SwarmFactory(-1.1, -1.1, -1.1, 1.1, 1.1, 1.1, n_source*n_source*n_source,
                                 distribution, seed);
  auto target_swarm = SwarmFactory(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n_target*n_target*n_target, 
                                  distribution, seed);
  
  auto source_state = make_shared<SwarmState<3>>(*source_swarm);
  auto target_state = make_shared<SwarmState<3>>(*target_swarm);

  int nsrcpts = source_swarm->num_particles();
  auto inputData = make_shared<typename SwarmState<2>::DblVec>(nsrcpts, 0.0);

  for (int p(0); p < nsrcpts; ++p) {
    Wonton::Point<3> coord = source_swarm->get_particle_coordinates(p);
    (*inputData)[p] = field_func<3>(example.field_order, coord);
  }
  
  source_state->add_field("remapdata", inputData);

  int ntarpts = target_swarm->num_particles();
  auto targetData = make_shared<typename SwarmState<2>::DblVec>(ntarpts, 0.0);
  target_state->add_field("remapdata", targetData);

  std::vector<std::string> remap_fields;
  remap_fields.push_back("remapdata");

  // Smoothing lengths (or in other words "size" of the support
  // function) in each dimension
  double h;
  Portage::Meshfree::WeightCenter center_type;
  int nsmooth;
  if(center==0) {
    center_type = Portage::Meshfree::Gather;
    h = 2.0*hfactor/n_target;
    nsmooth = ntarpts;
  } else if (center==1) {
    center_type = Portage::Meshfree::Scatter;
    h = 2.2*hfactor/n_source;
    nsmooth = nsrcpts;
  }
  auto smoothing_lengths =
    Portage::vector<std::vector<std::vector<double>>>
      (nsmooth, std::vector<std::vector<double>>(1, std::vector<double>(3, h)));
                                                                      
#ifdef HAVE_NANOFLANN  // Search by kdtree
  SwarmDriver<
    Search_KDTree_Nanoflann,
    Accumulate,
    Estimate,
    3,
    Swarm<3>,
    SwarmState<3>>
      d(*source_swarm, *source_state, *target_swarm, *target_state,
        smoothing_lengths, Portage::Meshfree::Weight::B4, 
        Portage::Meshfree::Weight::ELLIPTIC, center_type);
#else
  SwarmDriver<
    SearchPointsByCells,
    Accumulate,
    Estimate,
    3,
    Swarm<3>,
    SwarmState<3>>
      d(*source_swarm, *source_state, *target_swarm, *target_state,
        smoothing_lengths, Portage::Meshfree::Weight::B4, 
        Portage::Meshfree::Weight::ELLIPTIC, center_type);
#endif

  if (example.estimation_order == 0)
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Unitary);
  else if (example.estimation_order == 1)
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Linear);
  else if (example.estimation_order == 2)
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Quadratic);

#ifdef PORTAGE_ENABLE_MPI
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  d.run(&mpiexecutor, true);
#else
  d.run();
#endif

  std::vector<double> expected_value(ntarpts, 0.0);

  std::vector<double> toterr(3,0.0);
  for (int p(0); p < ntarpts; ++p) {
  Wonton::Point<3> coord = target_swarm->get_particle_coordinates(p);

    expected_value[p] = field_func<3>(example.field_order, coord);
    double error = fabs(expected_value[p] - (*targetData)[p]);
    toterr[0] = fmax(error, toterr[0]);
    toterr[1] = toterr[1]+error;
    toterr[2] = toterr[2]+error*error;
  }
  toterr[1] /= nsrcpts;
  toterr[2] /= nsrcpts;
  toterr[2] = sqrt(toterr[2]);

  std::cout << std::endl;
  std::cout << std::scientific;
  std::cout.precision(17);
  std::cout << "weight center = " << center << std::endl;
  std::cout << "distribution = " << distribution << std::endl;
  std::cout << "dimension = " << dimension << std::endl;
  std::cout << "seed = " << seed << std::endl;
  std::cout << std::endl;
  std::cout << "smoothing length = " << h << std::endl;
  std::cout << "Linf NORM OF ERROR = " << toterr[0] << std::endl;
  std::cout << "L1 NORM OF ERROR = " << toterr[1] << std::endl;
  std::cout << "L2 NORM OF ERROR = " << toterr[2] << std::endl;
  
  // Create the input file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream finp_csv("infield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".csv");
  finp_csv << std::scientific;
  finp_csv.precision(17);
  finp_csv << " X, Y, Z, Value\n";
  for (int p(0); p < nsrcpts; ++p) {
    Wonton::Point<3> coord = source_swarm->get_particle_coordinates(p);
    finp_csv << coord[0] << ", " << coord[1] << ", " << coord[2] << ", " << (*inputData)[p] << std::endl;
  }

  // Create the output file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream fout_csv("outfield-"
                     + std::to_string(static_cast<long long>(example_num)) + "-"
                     + std::to_string(static_cast<long long>(center))
                     + ".csv");
  fout_csv << std::scientific;
  fout_csv.precision(17);
  fout_csv << " X, Y, Z, Value\n";
  for (int p(0); p < ntarpts; ++p) {
    Wonton::Point<3> coord = target_swarm->get_particle_coordinates(p);
    fout_csv << coord[0] << ", " << coord[1] << ", " << coord[2] << ", " << (*targetData)[p] << std::endl;
  }
  std::cout << "finishing swarm app..." << std::endl;
}

int main(int argc, char** argv) {
  int example_num, n_source, n_target, distribution, dimension, center;
  double hfactor;
  unsigned int seed;
  if (argc < 7) {
    usage();
    return 0;
  }

  example_num = atoi(argv[1]);
  n_source = atoi(argv[2]);
  n_target = atoi(argv[3]);
  distribution = atoi(argv[4]);
  seed = atoi(argv[5]);
  if (example_num<10) dimension = 2;
  else                dimension = 3;
  hfactor = std::stod(argv[6]);
  center = std::atoi(argv[7]);

  // check inputs
  assert(example_num>=0 and example_num<20);
  assert(n_source>0 and n_target>0);
  assert(distribution>=0 and distribution<=2);
  assert(dimension==2 or dimension==3);
  assert(center==0 or center==1);

#ifdef PORTAGE_ENABLE_MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
#endif

  if (dimension==2) run_test_2d(example_num, n_source, n_target, distribution, dimension, seed, hfactor, center);
  if (dimension==3) run_test_3d(example_num, n_source, n_target, distribution, dimension, seed, hfactor, center);

#ifdef PORTAGE_ENABLE_MPI
  MPI_Finalize();
#endif
}
