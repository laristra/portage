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
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

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


double const const_val = 73.98;   // some random number

int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp " <<
      "--dim=2|3 --nsourcecells=N --ntargetcells=M --conformal=y|n \n" << 
      "--reverse_ranks=y|n --entity_kind=cell|node --field_order=0|1|2 \n" <<
      "--remap_order=1|2 --output_results=y|n \n\n";

  std::cout << "--dim (default = 2): spatial dimension of mesh\n\n";
  std::cout << "--nsourcecells (NO DEFAULT): Num cells per in each " <<
      "coord dir in source mesh\n\n";
  std::cout << "--ntargetcells (NO DEFAULT): Num of cells in each " <<
      "coord dir in target mesh\n\n";

  std::cout << "--conformal (default = y): 'y' means mesh boundaries match\n\n";

  std::cout << "--reverse_ranks (default = n): " <<
      "Applicable only for distributed runs.\n";
  std::cout << "  if 'y', then a greater mismatch is created " <<
      "between source and target meshes \n";
  std::cout << "  by reversing the assignment order of source mesh partitions to MPI ranks\n\n";

  std::cout << "--entity_kind (default = cell): entities on which remapping is to be done\n\n";

  std::cout << "--field order (default = 0): " <<
      "polynomial order of source field\n";
  std::cout << "  0th order function = " << const_val << "\n";
  std::cout << "  1st order function = x+y+z \n";
  std::cout << "  2nd order function = x*x+y*y+z*z \n";
  std::cout << "  Here, (x,y,z) is cell centroid or node coordinate\n\n";

  std::cout << "--remap order (default = 1): " <<
      "order of accuracy of interpolation\n\n";

  std::cout << "--output_results (default = y)" << std::endl;
  std::cout << "  If 'y', the two meshes are output with the " <<
      "remapped field attached. \n";
  std::cout << "  Also, the target field values are output to a text file called\n";
  std::cout << "  'field_cell_rA_fB.txt' or 'field_node_rA_rB.txt' where 'A' is\n";
  std::cout << "  the field polynomial order and B is the remap/interpolation order\n\n";
  return 0;
}
//////////////////////////////////////////////////////////////////////


int create_meshes(int const dim, int const n_source, int const n_target,
                  bool const conformal_meshes,
                  std::shared_ptr<Jali::Mesh> *sourceMesh,
                  std::shared_ptr<Jali::Mesh> *targetMesh,
                  MPI_Comm mpicomm) {

  int numpe, rank;
  MPI_Comm_size(mpicomm, &numpe);
  MPI_Comm_rank(mpicomm, &rank);

  // The mesh factory and mesh setup
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});

  if (dim == 2) {
    // 2d quad input mesh from (0,0) to (1,1) with n_source x n_source zones
    *sourceMesh = mf(0.0, 0.0, 1.0, 1.0, n_source, n_source);

    if (conformal_meshes) {
      // 2d quad output mesh from (0,0) to (1,1) with n_target x n_target zones
      *targetMesh = mf(0.0, 0.0, 1.0, 1.0, n_target, n_target);
    } else {
      // 2d quad output mesh from (0,0) to (1+1.5dx,1) with (n+1)x(n+1)
      // zones and dx equal to the sourceMesh grid spacing
      double dx = 1.0/static_cast<double>(n_target);
      *targetMesh = mf(0.0, 0.0, 1.0+1.5*dx, 1.0+1.5*dx, n_target, n_target);
    }
  } else if (dim == 3) {
    if (numpe == 1) {
      // 3d hex input mesh from (0,0,0) to (1,1,1) with n_source x
      // n_source x n_source zones

      *sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_source, n_source,
                       n_source);

      if (conformal_meshes) {
        // 3d hex output mesh from (0,0,0) to (1,1,1) with n_target
        // x n_target x n_target zones

        *targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_target, n_target,
                         n_target);
      } else {
        // 3d hex output mesh from (0,0,0) to
        // (1+1.5dx,1+1.5dx,1+1.5dx) with
        // (n_target)x(n_target)x(n_target) zones and dx equal to
        // the sourceMesh grid spacing

        double dx = 1.0/static_cast<double>(n_target);
        *targetMesh = mf(0.0, 0.0, 0.0, 1.0+1.5*dx, 1.0+1.5*dx, 1.0+1.5*dx,
                         n_target, n_target, n_target);
      }
    } else {  // generate the input and target meshes for the distributed case

      int source_dim = cbrt(1.0f*numpe) + 0.01f;

#ifdef MANUAL_SOURCE_DECOMPOSITION
      // Set up a local communicator so that we can define mesh partitions
      // explicitly on each rank without Jali distributing it for us
      MPI_Group world_group, local_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      int ranks[1];  ranks[0] = rank;
      MPI_Group_incl(world_group, 1, ranks, &local_group);
      MPI_Comm local_comm;
      MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm);
      Jali::MeshFactory mf_local(local_comm);

      mf_local.included_entities({Jali::Entity_kind::ALL_KIND});
      mf_local.boundary_ghosts_requested(false);
      mf_local.num_ghost_layers_distmesh(1);

      // compute the local partition of the source mesh based on the rank;
      // n_source is the number of cells in each dimension in each partition;
      // the number of ranks must be a perfect cube (1, 8, 27, etc.)
      double source_step = 1.0f / source_dim;
      int rrank = reverse_source_ranks ? numpe - rank - 1 : rank;
      int source_x = rrank % source_dim;
      int source_y = (rrank / source_dim) % source_dim;
      int source_z = rrank / (source_dim*source_dim);

      *sourceMesh = mf_local(source_step*source_x, source_step*source_y,
                             source_step*source_z, source_step*(source_x+1),
                             source_step*(source_y+1), source_step*(source_z+1),
                             n_source, n_source, n_source);

#else

      mf.included_entities({Jali::Entity_kind::FACE, Jali::Entity_kind::EDGE,
              Jali::Entity_kind::WEDGE});
      mf.boundary_ghosts_requested(false);
      mf.num_ghost_layers_distmesh(1);
      mf.partitioner(Jali::Partitioner_type::BLOCK);
      *sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_source*source_dim,
                       n_source*source_dim, n_source*source_dim);
#endif


      // n_target is the number of cells in each dimension in each partition;
      // the number of ranks must be a perfect cube (1, 8, 27, etc.)

      int target_dim = cbrt(1.0f*numpe) + 0.01f;
      if (conformal_meshes) {
        *targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_target*target_dim,
                         n_target*target_dim, n_target*target_dim);
      } else {
        double dx = 1.0/static_cast<double>(n_target*target_dim);
        *targetMesh = mf(0.0, 0.0, 0.0, 1.0+1.5*dx, 1.0+1.5*dx, 1.0+1.5*dx,
                         n_target*target_dim, n_target*target_dim,
                         n_target*target_dim);
      }

    }  // distributed case
  }

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


  // Get the example to run from command-line parameter
  int n_source = 0, n_target = 0;
  int dim = 2;
  bool conformal = true;

  int interp_order = 1;
  int poly_order = 0;
  bool dump_output = true;
  bool reverse_source_ranks = false;
  Jali::Entity_kind entityKind = Jali::Entity_kind::CELL;

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
    } else if (keyword == "dim") {
      dim = stoi(valueword);
      assert(dim == 2 || dim == 3);
    } else if (keyword == "nsourcecells")
      n_source = stoi(valueword);
    else if (keyword == "ntargetcells")
      n_target = stoi(valueword);
    else if (keyword == "conformal")
      conformal = (valueword == "y");
    else if (keyword == "reverse_ranks")
      reverse_source_ranks = (numpe > 1 && valueword == "y");
    else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
    } else if (keyword == "field_order") {
      poly_order = stoi(valueword);
      assert(poly_order >=0 && poly_order < 3);
    } else if (keyword == "output_results")
      dump_output = (valueword == "y");
    else
      std::cerr << "Unrecognized option " << keyword << std::endl;
  }

  if (numpe > 1 && entityKind == Jali::Entity_kind::NODE) {
    if (rank == 0)
      std::cerr << std::endl <<
          "portageapp_jali - ERROR - NODE CENTERED REMAP NOT IMPLEMENTED FOR DISTRIBUTED MESHES\n" << std::endl;
    MPI_Finalize();
    return -1;
  }

  if (rank == 0)
  {
    std::cout << "starting portageapp...\n";
    std::cout << "   Problem is in " << dim << " dimensions\n";
    std::cout << "   Source mesh has " << n_source << " cells in each dir\n";
    std::cout << "   Target mesh has " << n_target << " cells in each dir\n";
    std::cout << "   Field is of polynomial order " << poly_order << "\n";
    std::cout << "   Field lives on entity kind " << entityKind << "\n";
    std::cout << "   Interpolation order is " << interp_order << "\n";
  }

  struct timeval begin, end, diff;
  gettimeofday(&begin, 0);

  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;

  create_meshes(dim, n_source, n_target, conformal, &sourceMesh, &targetMesh,
                MPI_COMM_WORLD);

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
  Jali::State targetState(targetMesh);

  std::vector<double> sourceData;
  std::vector<std::string> remap_fields;

  // Cell-centered remaps
  if (entityKind == Jali::Entity_kind::CELL) {
    sourceData.resize(nsrccells);

    if (poly_order == 0) {
      for (unsigned int c = 0; c < nsrccells; ++c)
        sourceData[c] = const_val;
    } else if (poly_order == 1) {
      for (unsigned int c = 0; c < nsrccells; ++c) {
        JaliGeometry::Point cen = sourceMesh->cell_centroid(c);
        sourceData[c] = cen[0]+cen[1];
        if (dim == 3)
          sourceData[c] += cen[2];
      }
    } else {  // quadratic function
      for (unsigned int c = 0; c < nsrccells; ++c) {
        JaliGeometry::Point cen = sourceMesh->cell_centroid(c);
        sourceData[c] = cen[0]*cen[0]+cen[1]*cen[1];
        if (dim == 3)
          sourceData[c] += cen[2]*cen[2];
      }
    }

    sourceState.add("celldata", sourceMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    targetState.add("celldata", targetMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, 0.0);

    // Register the variable name and interpolation order with the driver
    remap_fields.push_back("celldata");

  } else {
    sourceData.resize(nsrcnodes);

    /*!
      @todo make node_get_coordinates be consistent in data type with
      cell_centroid?
    */
    if (poly_order == 0) {
      for (int i = 0; i < nsrcnodes; ++i)
        sourceData[i] = const_val;
    } else {
      if (dim == 2) {
        Portage::Point<2> nodexy;
        for (int i = 0; i < nsrcnodes; ++i) {
          sourceMeshWrapper.node_get_coordinates(i, &nodexy);
          if (poly_order == 1)
            sourceData[i] = nodexy[0] + nodexy[1];
          else
            sourceData[i] = nodexy[0]*nodexy[0] + nodexy[1]*nodexy[1];
        }
      } else {  // 3d
        Portage::Point<3> nodexyz;
        for (int i = 0; i < nsrcnodes; ++i) {
          sourceMeshWrapper.node_get_coordinates(i, &nodexyz);
          if (poly_order == 1)
            sourceData[i] = nodexyz[0] + nodexyz[1] + nodexyz[2];
          else
            sourceData[i] = nodexyz[0]*nodexyz[0] + nodexyz[1]*nodexyz[1]
                + nodexyz[2]*nodexyz[2];
        }
      }
    }

    sourceState.add("nodedata", sourceMesh, Jali::Entity_kind::NODE,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    targetState.add("nodedata", targetMesh, Jali::Entity_kind::NODE,
                    Jali::Entity_type::ALL, 0.0);


    // Register the variable name and remap order with the driver
    remap_fields.push_back("nodedata");

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
      d.set_remap_var_names(remap_fields);
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
      d.set_remap_var_names(remap_fields);
      d.run(numpe > 1);
    }
  }



  // Dump some timing information
  if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  const float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  if (rank == 0) std::cout << "Time: " << seconds << std::endl;



  // Output results for small test cases
  double error, toterr = 0.0;
  double const * cellvecout;
  double const * nodevecout;

  if (entityKind == Jali::Entity_kind::CELL) {
    targetStateWrapper.get_data<double>(Portage::CELL, "celldata",
                                        &cellvecout);

    if (numpe == 1 && n_target < 10)
      std::cout << "celldata vector on target mesh after remapping is:"
                << std::endl;

    if (dim == 2) {
      for (int c = 0; c < ntarcells; ++c) {
        Portage::Point<2> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
        
        if (poly_order == 0)
          error = const_val - cellvecout[c];
        else if (poly_order == 1)
          error = ccen[0] + ccen[1] - cellvecout[c];
        else  // quadratic
          error = ccen[0]*ccen[0] + ccen[1]*ccen[1] + ccen[2]*ccen[2] -
              cellvecout[c];

        if (numpe == 1 && n_target < 10) {
          std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                      ccen[0], ccen[1]);
          std::printf("  Value = % 10.6lf  Err = % lf\n",
                      cellvecout[c], error);
        }
        toterr += error*error;
      }
    } else {  // dim == 3
      for (int c = 0; c < ntarcells; ++c) {
        Portage::Point<3> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
        
        if (poly_order == 0)
          error = const_val - cellvecout[c];
        else if (poly_order == 1)
          error = ccen[0] + ccen[1] + ccen[2] - cellvecout[c];
        else  // quadratic
          error = ccen[0]*ccen[0] + ccen[1]*ccen[1] + ccen[2]*ccen[2] -
              cellvecout[c];
        
        if (numpe == 1 && n_target < 10) {
          std::printf("%d Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)",
                      rank, c, ccen[0], ccen[1], ccen[2]);
          std::printf("  Value = % 10.6lf  Err = % lf\n",
                      cellvecout[c], error);
        }
        toterr += error*error;
      } 
    }
    
  } else {

    targetStateWrapper.get_data<double>(Portage::NODE, "nodedata",
                                        &nodevecout);

    if (numpe == 1 && n_target < 10)
      std::cout << "nodedata vector on target mesh after remapping is:"
                << std::endl;

    if (dim == 2) {
      Portage::Point<2> nodexy;
      for (int i = 0; i < ntarnodes; ++i) {
        targetMeshWrapper.node_get_coordinates(i, &nodexy);

        if (poly_order == 0)
          error = const_val - nodevecout[i];
        else if (poly_order == 1)
          error = nodexy[0] + nodexy[1] - nodevecout[i];
        else
          error = nodexy[0]*nodexy[0] + nodexy[1]*nodexy[1] - nodevecout[i];

        if (n_target < 10) {
          std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf) ", i,
                      nodexy[0], nodexy[1]);
          std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], error);
        }

        toterr += error*error;
      }
    } else {  // dim == 3
      Portage::Point<3> nodexyz;
      for (int i = 0; i < ntarnodes; ++i) {
        targetMeshWrapper.node_get_coordinates(i, &nodexyz);
        
        if (poly_order == 0)
          error = const_val - nodevecout[i];
        else if (poly_order == 1)
          error = nodexyz[0] + nodexyz[1] + nodexyz[2] - nodevecout[i];
        else
          error = nodexyz[0]*nodexyz[0] + nodexyz[1]*nodexyz[1]
              + nodexyz[2]*nodexyz[2] - nodevecout[i];

        if (numpe == 1 && n_target < 10) {
          std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf,% 5.3lf) ", i,
                      nodexyz[0], nodexyz[1], nodexyz[2]);
          std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], error);
        }

        toterr += error*error;
      }
    }
  }

  if (numpe == 1) {
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
  } else {
    std::cout << std::flush << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    double globalerr;
    MPI_Reduce(&toterr, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(globalerr));
  }


  // Dump output, if requested
  if (dump_output) {

    // The current version of MSTK (2.27rc2) has a bug in writing out 
    // exodus files with node variables in parallel and so we will avoid
    // the exodus export in this situation. The 'if' statement can be
    // removed once we upgrade to the next version of MSTK

    if (numpe == 1) {
      if (rank == 0)
        std::cout << "Dumping data to Exodus files..." << std::endl;
      sourceState.export_to_mesh();
      targetState.export_to_mesh();
      sourceMesh->write_to_exodus_file("input.exo");
      targetMesh->write_to_exodus_file("output.exo");
      if (rank == 0)
        std::cout << "...done." << std::endl;
    }

    std::vector<int> lgid;
    std::vector<double> lvalues;
    std::string entstr;
    if (entityKind == Jali::Entity_kind::CELL) {
      entstr = "cell";
      lgid.resize(ntarcells);
      lvalues.resize(ntarcells);
      for (int i=0; i < ntarcells; i++) {
        lgid[i] = targetMesh->GID(i, Jali::Entity_kind::CELL);
        lvalues[i] = cellvecout[i];
      }
    } else {
      entstr = "node";
      lgid.resize(ntarnodes);
      lvalues.resize(ntarnodes);
      for (int i=0; i < ntarnodes; i++) {
        lgid[i] = targetMesh->GID(i, Jali::Entity_kind::NODE);
        lvalues[i] = nodevecout[i];
      }
    }

    // sort the field values by global ID

    std::vector<int> idx;
    argsort(lgid, idx);   // find sorting indices based on global IDS
    reorder(lgid, idx);   // sort the global ids
    reorder(lvalues, idx);  // sort the values

    // construct the field file name and open the file

    std::string fieldfilename = "field_" +
        std::to_string(static_cast<long long>(dim)) + "d_" +
        entstr + "_f" + std::to_string(static_cast<long long>(poly_order)) + "_r" +
        std::to_string(static_cast<long long>(interp_order));
    if (!conformal) fieldfilename = fieldfilename + "_nc";
    if (reverse_source_ranks) fieldfilename = fieldfilename + "_rev";
    fieldfilename = fieldfilename + ".txt";
    if (numpe > 1) {
      int maxwidth = static_cast<long long>(std::ceil(std::log10(numpe)));
      char rankstr[10];
      std::snprintf(rankstr, sizeof(rankstr), "%0*d", maxwidth, rank);
      fieldfilename = fieldfilename + "." + std::string(rankstr);
    }
    std::ofstream fout(fieldfilename);
    fout << std::scientific;
    fout.precision(17);

    // write out the values

    for (int i=0; i < lgid.size(); i++)
      fout << lgid[i] << " " << lvalues[i] << std::endl;
  }  // if (dump_output)

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}
