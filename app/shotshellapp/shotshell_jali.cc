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





#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "mpi.h"

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

#include "portage/support/mpi_collate.h"
#include "portage/driver/driver.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_r2d.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#define MSTK_HAVE_MPI
#include "Mesh_MSTK.hh"

using Portage::Jali_Mesh_Wrapper;
using Portage::Jali_State_Wrapper;
using Portage::argsort;
using Portage::reorder;

double const_val = 90.36;

int print_usage() {
  std::cout << "Usage: shotshellapp --entity_kind=cell|node " <<
      "--remap_order=1|2 --field_order=0|1|2 --output_results=y|n " <<
      " --source=input-mesh --target=output-mesh" << std::endl;
  std::cout << "      remap order - default = 1" << std::endl << std::endl;
  std::cout << "      field order - default = 0" << std::endl << std::endl;
  std::cout << "      ---- 0th order source function = " << const_val <<
      std::endl;
  std::cout << "      ---- 1st order source function = x+y+z where (x,y,z) " <<
      std::endl;
  std::cout << "           is cell centroid or node coordinate" << std::endl;
  std::cout << "      ---- 2nd order source function = x*x+y*y+z*z where " <<
      std::endl;
  std::cout << "           (x,y,z) is cell centroid or node coordinate" <<
      std::endl << std::endl;
  std::cout << "      output_results - default = y" << std::endl;
  std::cout << "      ---- if output is requested, the target mesh is written"
            << std::endl;
  std::cout << "           with the remapped field attached. Also, a text file"
            << std::endl;
  std::cout << "           called field_cell_rA_fB.txt or field_node_rA_rB.txt"
            << std::endl;
  std::cout << "           containing target field values is written where " <<
      std::endl;
  std::cout << "           'A' is the remap order and B is the field order" <<
      std::endl;
  return 0;
}

// This is a 2-D test!  Find the example data in
// test_data/shotshell.exo, shotshell-v.exo.

int main(int argc, char** argv) {
  const double TOL = 1e-4;

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


  int interp_order = 1;
  int poly_order = 0;
  bool dumpMesh = true;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;
  std::string infilename;
  std::string outfilename;

  // Get the example to run from command-line parameter
  if (argc < 3) return print_usage();

  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of("=");
    std::string keyword = arg.substr(keyword_beg, keyword_end-keyword_beg);
    std::string valueword = arg.substr(keyword_end+1, len-(keyword_end+1));
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
  } else if (keyword == "remap_order")
      interp_order = stoi(valueword);
    else if (keyword == "field_order")
      poly_order = stoi(valueword);
    else if (keyword == "output_results")
      dumpMesh == (valueword == "y");
    else if (keyword == "source")
      infilename = valueword;
    else if (keyword == "target")
      outfilename = valueword;
    else
      std::cerr << "Unrecognized option " << keyword << std::endl;
  }

  if (numpe > 1 && entityKind == Jali::Entity_kind::NODE) {
    if (rank == 0)
      std::cerr << std::endl <<
          "shotshellapp_jali - ERROR - NODE-CENTERED REMAPPING NOT IMPLEMENTED FOR DISTRIBUTED MESHES" << std::endl;
    MPI_Finalize();
    return -1;
  }
    
      
  std::cout << "starting shotshellapp...\n";
  std::cout << "   Field is of polynomial order" << poly_order << "\n";
  std::cout << "   Field lives on entity kind " << entityKind << "\n";
  std::cout << "   Interpolation order is " << interp_order << "\n";

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE,
                        Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  const std::shared_ptr<Jali::Mesh> sourceMesh = mf(infilename);
  const Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  const int inputDim = sourceMesh->space_dimension();

  const std::shared_ptr<Jali::Mesh> targetMesh = mf(outfilename);
  const Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  const int targetDim = targetMesh->space_dimension();

  assert(inputDim == targetDim);

  std::cout << "Target mesh stats: " <<
      targetMeshWrapper.num_owned_cells() << " " <<
      targetMeshWrapper.num_owned_nodes() << std::endl;

  Jali::State sourceState(sourceMesh);
  Jali::State targetState(targetMesh);

  std::vector<double> sourceData;
  std::vector<std::string> remap_fields;

  int ncall_src = sourceMeshWrapper.num_owned_cells() +
      sourceMeshWrapper.num_ghost_cells();
  int nnall_src = sourceMeshWrapper.num_owned_nodes() +
        sourceMeshWrapper.num_ghost_nodes();
    

  if (entityKind == Jali::Entity_kind::CELL) {
    sourceData.resize(ncall_src, 0.0);

    for (int i = 0; i < ncall_src; i++) {
      JaliGeometry::Point coord = sourceMesh->cell_centroid(i);
      double x = coord[0];
      double y = coord[1];
      double z = (inputDim > 2) ? coord[2] : 0.0;
      if (poly_order == 0)
        sourceData[i] = const_val;
      else if (poly_order == 1)
        sourceData[i] = x + y + z;
      else if (poly_order == 2)
        sourceData[i] = x*x + y*y + z*z;
    }

    sourceState.add("celldata", sourceMesh, entityKind,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    targetState.add("celldata", targetMesh, entityKind,
                    Jali::Entity_type::ALL, 0.0);

    remap_fields.push_back("celldata");

  } else {
    sourceData.resize(nnall_src, 0.0);

    for (int i = 0; i < nnall_src; i++) {
      JaliGeometry::Point nxyz;
      sourceMesh->node_get_coordinates(i, &nxyz);
      double x = nxyz[0];
      double y = nxyz[1];
      double z = (inputDim > 2) ? nxyz[2] : 0.0;
      if (poly_order == 0)
        sourceData[i] = const_val;
      else if (poly_order == 1)
        sourceData[i] = x+y+z;
      else if (poly_order == 2)
        sourceData[i] = x*x+y*y+z*z;
    }

    sourceState.add("nodedata", sourceMesh, entityKind,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    targetState.add("nodedata", targetMesh, entityKind,
                    Jali::Entity_type::ALL, 0.0);

    remap_fields.push_back("nodedata");
  }


  const Jali_State_Wrapper sourceStateWrapper(sourceState);
  Jali_State_Wrapper targetStateWrapper(targetState);


  if (inputDim == 2) {
    if (interp_order == 1) {
      Portage::Driver<Portage::SearchKDTree,
                      Portage::IntersectR2D,
                      Portage::Interpolate_1stOrder,
                      2,
                      Portage::Jali_Mesh_Wrapper,
                      Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    } else if (interp_order == 2) {
      Portage::Driver<Portage::SearchKDTree,
                      Portage::IntersectR2D,
                      Portage::Interpolate_2ndOrder,
                      2,
                      Portage::Jali_Mesh_Wrapper,
                      Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    }
  } else if (inputDim == 3) {
    if (interp_order == 1) {
      Portage::Driver<Portage::SearchKDTree,
                      Portage::IntersectR3D,
                      Portage::Interpolate_1stOrder,
                      3,
                      Portage::Jali_Mesh_Wrapper,
                      Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    } else if (interp_order == 2) {
      Portage::Driver<Portage::SearchKDTree,
                      Portage::IntersectR3D,
                      Portage::Interpolate_2ndOrder,
                      3,
                      Portage::Jali_Mesh_Wrapper,
                      Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    }
  }


  double l2norm = 0.0;
  double const * cellvecout;
  double const * nodevecout;
  int ncown_trg = targetMeshWrapper.num_owned_cells();
  int nnown_trg = targetMeshWrapper.num_owned_nodes();

  if (entityKind == Jali::Entity_kind::CELL) {

    targetStateWrapper.get_data<double>(Portage::CELL, "celldata", &cellvecout);
    std::cerr << "Last result: " << cellvecout[ncown_trg] << std::endl;
    
    for (int c = 0; c < ncown_trg; c++) {
      JaliGeometry::Point ccen;
      double x, y, z = 0.0;
      double diff;
      double expected_val = 0.0;
      ccen = targetMesh->cell_centroid(c);
      x = ccen[0];
      y = ccen[1];
      z = (inputDim == 3) ? ccen[2] : 0.0;
      if (poly_order == 0)
        expected_val = const_val;
      else if (poly_order == 1)
        expected_val = x + y + z;
      else if (poly_order == 2)
        expected_val = x*x + y*y + z*z;

      diff = expected_val - cellvecout[c];
      l2norm += diff*diff;
    }
  } else {

    targetStateWrapper.get_data<double>(Portage::NODE, "nodedata", &nodevecout);

    std::cerr << "Last result: " << nodevecout[nnown_trg-1] <<
        std::endl;

    for (int n = 0; n < nnown_trg; n++) {
      JaliGeometry::Point nxyz;
      double x, y, z = 0.0;
      double diff;
      double expected_val = 0.0;
      targetMesh->node_get_coordinates(n, &nxyz);
      x = nxyz[0];
      y = nxyz[1];
      z = (inputDim == 3) ? nxyz[2] : 0.0;
      if (poly_order == 0)
        expected_val = const_val;
      else if (poly_order == 1)
        expected_val = x + y + z;
      else if (poly_order == 2)
        expected_val = x*x + y*y + z*z;

      diff = expected_val - nodevecout[n];
      l2norm += diff*diff;
    }
  }

  if (numpe == 1) {
    l2norm = sqrt(l2norm);
    std::cout << "L2 norm of error " << l2norm << "\n";
  } else {
    std::cout << std::flush << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    double globalerr;
    MPI_Reduce(&l2norm, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) std::cerr << "L2 norm of error " << sqrt(globalerr) << "\n";
  }
    

  if (dumpMesh) {
    std::string entstr = (entityKind == Jali::Entity_kind::CELL) ? "cell" :
        "node";

    // current version of MSTK (2.27rc2) has a bug in exporting exodus
    // files with node based variables in parallel. Can remove the
    // following 'if' condition when Jali TPLs are upgraded

    if (numpe == 1) {
      std::cerr << "Saving the source mesh" << std::endl;
      sourceState.export_to_mesh();
      
      std::string inbasename =
          infilename.substr(0, infilename.find(".exo")) + "_" + entstr +
          "_f" + std::to_string(static_cast<long long>(poly_order)) + "_r" +
          std::to_string(static_cast<long long>(interp_order));
      
      dynamic_cast<Jali::Mesh_MSTK*>(sourceMesh.get())->
          write_to_exodus_file(inbasename + ".exo");
      
      std::cerr << "Saving the target mesh" << std::endl;
      targetState.export_to_mesh();
      
      std::string outbasename =
          outfilename.substr(0, outfilename.find(".exo")) + "_" + entstr +
          "_f" + std::to_string(static_cast<long long>(poly_order)) + "_r" +
          std::to_string(static_cast<long long>(interp_order));
      
      dynamic_cast<Jali::Mesh_MSTK*>(targetMesh.get())->
          write_to_exodus_file(outbasename + ".exo");
    }
      
    // Collect and sort field values by global IDs of entities
    
    std::vector<int> lgid;
    std::vector<double> lvalues;

    if (entityKind == Jali::Entity_kind::CELL) {
      lgid.resize(ncown_trg);
      lvalues.resize(ncown_trg);
      for (int i = 0; i < ncown_trg; i++) {
        lgid[i] = targetMesh->GID(i, entityKind);
        lvalues[i] = cellvecout[i];
      }
    } else {
      lgid.resize(nnown_trg);
      lvalues.resize(nnown_trg);
      for (int i = 0; i < nnown_trg; i++) {
        lgid[i] = targetMesh->GID(i, entityKind);
        lvalues[i] = nodevecout[i];
      }
    }

    std::vector<int> idx;
    argsort(lgid, idx);
    reorder(lgid, idx);
    reorder(lvalues, idx);
    
    // Build filename and open the file
    std::string fieldfilename = "field_" +
        std::to_string(static_cast<long long>(inputDim)) + "d_" +
        entstr + "_f" + std::to_string(static_cast<long long>(poly_order)) + "_r" +
        std::to_string(static_cast<long long>(interp_order)) + ".txt";
    if (numpe > 1) {
      int maxwidth = static_cast<long long>(std::ceil(std::log10(numpe)));
      char rankstr[10];
      std::snprintf(rankstr, sizeof(rankstr), "%0*d", maxwidth, rank);
      fieldfilename = fieldfilename + "." + std::string(rankstr);
    }
    std::ofstream fout(fieldfilename);
    fout << std::scientific;
    fout.precision(17);
    for (int i=0; i < lgid.size(); i++)
      fout << lgid[i] << " " << lvalues[i] << std::endl;
  }

  std::printf("finishing shotshellapp...\n");

  MPI_Finalize();

  return 0;
}
