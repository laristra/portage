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




#include <sys/time.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <utility>
#include <cmath>

#include <mpi.h>

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/mpi_collate.h"
#include "portage/driver/driver.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/jali/jali_state_wrapper.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

using Portage::Jali_Mesh_Wrapper;
using Portage::argsort;
using Portage::reorder;

/*!
  @file portageapp_jali.cc
  @brief A simple application that drives our remap routines.

  This program is used to showcase our capabilities with various types
  of remap operations (e.g. interpolation order) on various types of
  meshes (2d or 3d; node-centered or zone-centered) of some simple
  linear or quadratic data.  For the cases of remapping linear data
  with a second-order interpolator, the L2 norm output at the end
  should be identically zero.
 */

//////////////////////////////////////////////////////////////////////


int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp \n";
  std::cout << "--source_file=srcfilename: file name of source mesh (Exodus II format only)\n";
  std::cout << "--target_file=trgfilename: file name of target mesh (Exodus II format only)\n";

  std::cout << "--remap_order=1|2 (default = 1): order of accuracy of interpolation\n\n";
  std::cout << "--limit=0|1 (default = 1): limiting to preserve bounds (1) or no (0)";
  std::cout << "--entity_kind=cell|node (default = cell): kind of entity remap fields live on\n\n";

  std::cout << "  Also, the target mesh is written out with the attached field values are output to a file 'output.exo'\n";
  return 0;
}


int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif


  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int dim = 0;
  int interp_order = 1;
  int poly_order = 0;
  std::string srcfile, trgfile;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  Portage::LimiterType limiter = Portage::NOLIMITER;

  if (argc < 3) return print_usage();
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of("=");
    std::string keyword = arg.substr(keyword_beg, keyword_end-keyword_beg);
    std::string valueword = (keyword_end == std::string::npos) ? "" :
        arg.substr(keyword_end+1, len-(keyword_end+1));

    if (keyword == "entity_kind") {
      if (valueword == "cell" || valueword == "CELL")
        entityKind = Jali::Entity_kind::CELL;
      else if (valueword == "node" || valueword == "NODE")
        entityKind = Jali::Entity_kind::NODE;
      else {
        std::cerr <<
            "Only node and cell based remapping supported at this time" <<
            std::endl;
        exit(-1);
      }
    }
    else if (keyword == "source_file")
      srcfile = valueword;
    else if (keyword == "target_file")
      trgfile = valueword;
    else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
    }
    else if (keyword == "limit") {
      if (valueword == "1" || valueword == "")
        limiter = Portage::BARTH_JESPERSEN;
    }
    else
      std::cerr << "Unrecognized option " << keyword << std::endl;
  }

  if ((srcfile.length() == 0 && trgfile.length() != 0) ||
      (srcfile.length() != 0 && trgfile.length() == 0)) {
    if (rank == 0)
      std::cerr <<
          "portageapp_jali_file - ERROR - both source_file and target_file must be specified\n" <<
          std::endl;
  }

  if (rank == 0) {
    std::cout << "starting portageapp_jali_file...\n";
    std::cout << "Source mesh file is " << srcfile << "\n";
    std::cout << "Target mesh file is " << trgfile << "\n";
    std::cout << "   Field lives on entity kind " << entityKind << "\n";
    std::cout << "   Interpolation order is " << interp_order << "\n";
    if (interp_order == 2)
      std::cout << "   Limiter type is " << limiter << "\n";
  }

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE,
          Jali::Entity_kind::EDGE,
          Jali::Entity_kind::WEDGE,
          Jali::Entity_kind::CORNER});
  sourceMesh = mf(srcfile);
  targetMesh = mf(trgfile);
  dim = sourceMesh->space_dimension();
  if (dim != targetMesh->space_dimension()) {
    if (rank == 0)
      std::cerr << "Source and target mesh dimensions are different\n";
    MPI_Finalize();
    return -1;
  }

  // Wrappers for interfacing with the underlying mesh data structures.
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  const int nsrccells = sourceMeshWrapper.num_owned_cells() +
      sourceMeshWrapper.num_ghost_cells();
  const int ntarcells = targetMeshWrapper.num_owned_cells();

  const int nsrcnodes = sourceMeshWrapper.num_owned_nodes() +
      sourceMeshWrapper.num_ghost_nodes();
  const int ntarnodes = targetMeshWrapper.num_owned_nodes();

  // Native jali state managers for source and target
  Jali::State sourceState(sourceMesh);
  sourceState.init_from_mesh();
  Jali::State targetState(targetMesh);

  std::vector<double> sourceData;
  std::vector<std::string> remap_fields;

  Jali::State::iterator it = sourceState.begin();
  if (entityKind == Jali::Entity_kind::CELL) {  // Cell-centered remaps
    while (it != sourceState.end()) {
      std::shared_ptr<Jali::BaseStateVector> bv = *it;
      if (bv->entity_kind() == Jali::Entity_kind::CELL &&
          bv->get_type() == typeid(double)) {
        remap_fields.push_back(bv->name());

        targetState.add(bv->name(), targetMesh, Jali::Entity_kind::CELL,
                        Jali::Entity_type::ALL, 0.0);
      }
      ++it;
    }
  } else {  // Node-centered remaps
    while (it != sourceState.end()) {
      std::shared_ptr<Jali::BaseStateVector> bv = *it;
      if (bv->entity_kind() == Jali::Entity_kind::NODE) {
        remap_fields.push_back(bv->name());

        targetState.add(bv->name(), targetMesh, Jali::Entity_kind::NODE,
                        Jali::Entity_type::ALL, 0.0);
      }
      ++it;
    }
  }

  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  const float seconds_init = diff.tv_sec + 1.0E-6*diff.tv_usec;
  if (rank == 0) std::cout << "Mesh Initialization Time: " << seconds_init <<
                     std::endl;
  
  gettimeofday(&begin, 0);


  // Portage wrappers for source and target fields

  Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);
  Portage::Jali_State_Wrapper targetStateWrapper(targetState);

  if (dim == 2) {
    if (interp_order == 1) {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR2D,
        Portage::Interpolate_1stOrder,
        2,
        Portage::Jali_Mesh_Wrapper,
        Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    } else if (interp_order == 2) {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR2D,
        Portage::Interpolate_2ndOrder,
        2,
        Portage::Jali_Mesh_Wrapper,
        Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields, limiter);
      d.run(numpe > 1);
    }
  } else {
    if (interp_order == 1) {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_1stOrder,
        3,
        Portage::Jali_Mesh_Wrapper,
        Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    } else {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_2ndOrder,
        3,
        Portage::Jali_Mesh_Wrapper,
        Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields, limiter);
      d.run(numpe > 1);
    }
  }



  // Dump some timing information
  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  const float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  if (rank == 0) std::cout << "Time: " << seconds << std::endl;



  // Dump output
  targetState.export_to_mesh();
  targetMesh->write_to_exodus_file("output.exo");
  if (rank == 0)
    std::cout << "...done." << std::endl;

  if (rank == 0)
    std::printf("finishing portageapp_jali_file...\n");

  MPI_Finalize();
}
