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

#ifdef ENABLE_MPI
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
#endif

using std::shared_ptr;
using std::make_shared;
using Portage::Meshfree::Swarm;
using Portage::Meshfree::SwarmState;
using Portage::Meshfree::SwarmDriver;
using Portage::Meshfree::SwarmFactory;
using Portage::Meshfree::Accumulate;
using Portage::Meshfree::Estimate;
using Portage::SearchSimplePoints;
#ifdef HAVE_NANOFLANN
using Portage::Search_KDTree_Nanoflann;
#endif
using Portage::SearchPointsByCells;


//////////////////////////////////////////////////////////////////////
// Helper routines and data structures

template<unsigned int D>
double field_func(int field_order, Portage::Point<D> coord) {
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
  example_properties(const int field_order, const int estimation_order)
      : field_order(field_order), estimation_order(estimation_order)
  { }
  int estimation_order;   // estimation order in example
  int field_order;        // order of the field to map
};

// Use this to add new problems.  If needed, we can extend the
// example_properties struct to contain more information.
std::vector<example_properties> setup_examples() {
  std::vector<example_properties> examples;

  // const func, unitary basis
  examples.emplace_back(0, 0);

  // linear func, unitary basis
  examples.emplace_back(1, 0);

  // linear func, linear basis
  examples.emplace_back(1, 1);

  // quadratic func, unitary basis
  examples.emplace_back(2, 0);

  // quadratic func, linear basis
  examples.emplace_back(2, 1);

  // quadratic func, quadratic basis
  examples.emplace_back(2, 2);

  // cubic func, unitary basis
  examples.emplace_back(3, 0);

  // cubic func, linear basis
  examples.emplace_back(3, 1);

  // cubic func, quadratic basis
  examples.emplace_back(3, 2);

  // exp(-1.5*r)*sin(0.9*2*pi*x)*sin(0.9*2*pi*y)
  examples.emplace_back(-1, 2);

  return examples;
}

void usage() {
  auto examples = setup_examples();
  std::cout << "Usage: swarmapp example-number nsourcepts ntargetpts distribution"
            << std::endl;
  std::cout << "distribution = 0 for random, 1 for regular grid, 2 for perturbed regular grid"
            << std::endl;
  std::cout << "List of example numbers:" << std::endl;
  int i = 0;
  bool separated = false;
  for (const auto &example : examples) {
    std::printf("  %d: order %d remap of order %d func\n",
                i, example.estimation_order,
                example.field_order);
    ++i;
  }
}

int main(int argc, char** argv) {
  int example_num, n_source, n_target, distribution;
  if (argc <= 4) {
    usage();
    return 0;
  }

  example_num = atoi(argv[1]);
  n_source = atoi(argv[2]);
  n_target = atoi(argv[3]);
  distribution = atoi(argv[4]);

#ifdef ENABLE_MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
#endif

  std::cout << "starting swarm app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  example_properties example = setup_examples()[example_num];

  // Regularly ordered input swarm; randomly ordered output swarm
  auto inputSwarm = SwarmFactory(-1.1, -1.1, 1.1, 1.1, n_source*n_source,
                                 distribution);
  auto targetSwarm = SwarmFactory(-1.0, -1.0, 1.0, 1.0, n_target*n_target, distribution);
  
  auto inputState = make_shared<SwarmState<2>>(*inputSwarm);
  auto targetState = make_shared<SwarmState<2>>(*targetSwarm);

  int ninpts = inputSwarm->num_particles();
  auto inputData = make_shared<typename SwarmState<2>::DblVec>(ninpts, 0.0);

  for (int p(0); p < ninpts; ++p) {
    Portage::Point<2> coord = inputSwarm->get_particle_coordinates(p);
    (*inputData)[p] = field_func<2>(example.field_order, coord);
  }
  
  inputState->add_field("remapdata", inputData);

  int ntarpts = targetSwarm->num_particles();
  auto targetData = make_shared<typename SwarmState<2>::DblVec>(ntarpts, 0.0);
  targetState->add_field("remapdata", targetData);

  std::vector<std::string> remap_fields;
  remap_fields.push_back("remapdata");

  // Smoothing lengths (or in other words "size" of the support
  // function) in each dimension
  double h = 2.0/n_target;
  auto smoothing_lengths =
    Portage::vector<std::vector<std::vector<double>>>
      (ntarpts, std::vector<std::vector<double>>(1, std::vector<double>(2, 2.05*h)));
                                                                      
#ifdef HAVE_NANOFLANN  // Search by kdtree
  SwarmDriver<
    Search_KDTree_Nanoflann,
    Accumulate,
    Estimate,
    2,
    Swarm<2>,
    SwarmState<2>>
      d(*inputSwarm, *inputState, *targetSwarm, *targetState,
        smoothing_lengths);
#else
  SwarmDriver<
    SearchPointsByCells,
    Accumulate,
    Estimate,
    2,
    Swarm<2>,
    SwarmState<2>>
      d(*inputSwarm, *inputState, *targetSwarm, *targetState,
        smoothing_lengths);
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
  d.run(true, true);


  std::vector<double> expected_value(ntarpts, 0.0);

  double toterr = 0.0;
  for (int p(0); p < ntarpts; ++p) {
  Portage::Point<2> coord = targetSwarm->get_particle_coordinates(p);

    expected_value[p] = field_func<2>(example.field_order, coord);
    double error = expected_value[p] - (*targetData)[p];

    if (ntarpts < 10) {
      std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf)", p,
                  coord[0], coord[1]);
      double dummy=(*targetData)[p];
      std::printf("  Value = % 10.6lf  Err = % lf\n", dummy, error);
    }
    toterr += error*error;
  }

  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
  
  // Create the input file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream finp("infield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".txt");
  finp << std::scientific;
  finp.precision(17);
  std::ofstream finp_csv("infield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".csv");
  finp_csv << std::scientific;
  finp_csv.precision(17);
  finp_csv << " X, Y, Value\n";
  for (int p(0); p < ninpts; ++p) {
    Portage::Point<2> coord = inputSwarm->get_particle_coordinates(p);
    
    finp << coord[0] << " " << coord[1] << " " << (*inputData)[p] << std::endl;
    finp_csv << coord[0] << ", " << coord[1] << ", " << (*inputData)[p] << std::endl;
  }

  // Create the output file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream fout("outfield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".txt");
  fout << std::scientific;
  fout.precision(17);
  std::ofstream fout_csv("outfield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".csv");
  fout_csv << std::scientific;
  fout_csv.precision(17);
  fout_csv << " X, Y, Value\n";
  for (int p(0); p < ntarpts; ++p) {
    Portage::Point<2> coord = targetSwarm->get_particle_coordinates(p);
    
    fout << coord[0] << " " << coord[1] << " " << (*targetData)[p] << std::endl;
    fout_csv << coord[0] << ", " << coord[1] << ", " << (*targetData)[p] << std::endl;
  }
  std::cout << "finishing swarm app..." << std::endl;

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
}
