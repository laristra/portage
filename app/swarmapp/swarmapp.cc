/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
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
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

using std::vector;
using std::shared_ptr;
using std::make_shared;
using Portage::Meshfree::Swarm;
using Portage::Meshfree::SwarmState;
using Portage::Meshfree::SwarmDriver;
using Portage::Meshfree::SwarmFactory;
using Portage::Meshfree::Accumulate;
using Portage::Meshfree::Estimate;
using Portage::SearchSimplePoints;


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
      //      value = exp(4*rsqr);
      value = 1.0;
      for (int i = 0; i < D; i++)
        value *= sin(2*M_PI*coord[i]);
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
  int estimation_order;           // estimation order in example
  int field_order;     // order of the field to map
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

  // exp(4*r^2)*sin(2*pi*x)*sin(2*pi*y)
  examples.emplace_back(-1, 2);

  return examples;
}

void usage() {
  auto examples = setup_examples();
  std::cout << "Usage: swarmapp example-number nsourcepts ntargetpts"
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
  int example_num, n_source, n_target;
  if (argc <= 3) {
    usage();
    return 0;
  }

  example_num = atoi(argv[1]);
  n_source = atoi(argv[2]);
  n_target = atoi(argv[3]);

  // Even though Simple_Mesh is serial only, we still need to initialize MPI
  // for other Portage code.

#ifdef ENABLE_MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  if (numpe > 1) {
    std::cerr << "Particle Remap is only designed for a single process!"
              << std::endl;
    return 1;
  }
#endif

  std::cout << "starting swarm app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  example_properties example = setup_examples()[example_num];

  // Regularly ordered input swarm; randomly ordered output swarm
  auto inputSwarm = SwarmFactory(-1.25, -1.25, 1.25, 1.25, n_source*n_source,
                                 1);
  auto targetSwarm = SwarmFactory(-1.0, -1.0, 1.0, 1.0, n_target*n_target, 2);
  
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
      vector<vector<vector<double>>>(ntarpts, vector<vector<double>>(1, vector<double>(2, 2.05*h)));
                                                                      

  SwarmDriver<
    SearchSimplePoints,
    Accumulate,
    Estimate,
    2,
    Swarm<2>,
    SwarmState<2>>
      d(*inputSwarm, *inputState, *targetSwarm, *targetState,
        smoothing_lengths);
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
  d.run(false);

  std::vector<double> expected_value(ntarpts, 0.0);

  double toterr = 0.0;
  for (int p(0); p < ntarpts; ++p) {
  Portage::Point<2> coord = targetSwarm->get_particle_coordinates(p);

    expected_value[p] = field_func<2>(example.field_order, coord);
    double error = expected_value[p] - (*targetData)[p];

    std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf)", p,
                coord[0], coord[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n", (*targetData)[p], error);
    
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
  for (int p(0); p < ninpts; ++p) {
    Portage::Point<2> coord = inputSwarm->get_particle_coordinates(p);
    
    finp << coord[0] << " " << coord[1] << " " << (*inputData)[p] << std::endl;
  }

  // Create the output file for comparison.
  // The `static_cast` is a workaround for Intel compiler, which is missing a
  // `to_string` function for ints.
  std::ofstream fout("outfield"
                     + std::to_string(static_cast<long long>(example_num))
                     + ".txt");
  fout << std::scientific;
  fout.precision(17);
  for (int p(0); p < ntarpts; ++p) {
    Portage::Point<2> coord = targetSwarm->get_particle_coordinates(p);
    
    fout << coord[0] << " " << coord[1] << " " << (*targetData)[p] << std::endl;
  }
  std::cout << "finishing swarm app..." << std::endl;

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
}
