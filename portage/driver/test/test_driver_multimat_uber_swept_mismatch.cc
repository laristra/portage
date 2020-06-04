/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include "portage/support/portage.h"
#ifdef PORTAGE_HAS_TANGRAM

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#ifdef WONTON_ENABLE_MPI
#include "mpi.h"
#endif

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/driver/uberdriver.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/search/search_swept_face.h"


#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/search/search_kdtree.h"

double TOL = 1e-6;


// Tests for multi-material swept-face remap

// The conceptual material layout that includes a T-junction and
// interfaces that are aligned with the coordinate axes/planes. The
// three materials have very disparate but constant densities. Given
// this setup we know the volume, mass and centroid of each material.

// In each of the four tests, the rightmost boundary nodes are moved
// by some delta amount outside the initial domain of the source mesh.
// The source domain along each direction is from 0.0 to 1.0.
// The target domain along each direction is from 0.0 to 1.0+DX. 
//
// The current behavior of the swept-face based intersection skips the 
// facet-connected neighbor of a source cell if that facet is on the 
// domain boundary. As a result, the swept-polygon covering the mismatched 
// area is not accounted for in the list of weights for interpolation. 
// Therefore, the remapped volume fraction is partial and the exact values that
// are used for comparison accounts for this partial volume fraction. 
// These tests confirm that the swept face remap, in its current form, 
// still need a mismatch fix. This needs to be addressed later !!  

TEST(UberDriverSwept, ThreeMat2D_1stOrder_MisMatch) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 3, 3);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 3, 3);

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  // move a single node of the target mesh in x direction by DX
  double DX = 0.01;
  std::array<double, 2> pnt {};
  targetMesh->node_get_coordinates(5, &pnt);
  pnt[0] += DX;
  targetMesh->node_set_coordinates(5, pnt.data());

  //move the right boundary nodes for the target, so 
  //that the last layer of target cells are partially
  //covered by the source mesh
  for (int i = 0; i < 4; i++) {
   targetMesh->node_get_coordinates(12+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(12+i, pnt.data());
  }

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  // The material geometry in the overall domain will look like this
  // and we will put down a rectangular mesh that has multiple cells
  // in each direction on this domain so that we get some pure and
  // some mixed cells
  //
  // Note that only MOF type algorithms or VOF algorithms with the
  // material ordering 0,1,2 will get the T-junction geometry right
  //
  //    0,1           0.5,1         1,1
  //     *-------------:------------*
  //     |             :            |
  //     |             :        2   |
  //     |             :     mat2   |
  //     |             :            |
  //     |             :            |
  //     |     0       +............|1,0.5
  //     |   mat0      :            |
  //     |             :            |
  //     |             :     mat1   |
  //     |             :        1   |
  //     |             :            |
  //     *-------------:------------*
  //    0,0           0.5,0         1,0

  constexpr int nmats = 3;
  std::string matnames[nmats] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain

  Wonton::Point<2> matlo[nmats], mathi[nmats];
  matlo[0] = Wonton::Point<2>(0.0, 0.0);
  mathi[0] = Wonton::Point<2>(0.5, 1.0);
  matlo[1] = Wonton::Point<2>(0.5, 0.0);
  mathi[1] = Wonton::Point<2>(1.0, 0.5);
  matlo[2] = Wonton::Point<2>(0.5, 0.5);
  mathi[2] = Wonton::Point<2>(1.0, 1.0);

  double matrho[nmats] = {0.1, 10.0, 100.0};  // material density
  double matvol[nmats] = {0.5, 0.25, 0.25};
  double matmass[nmats] = {0.05, 2.5, 25.0};

  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Wonton::Point<2>> matcen_src[nmats];


  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  //
  // Based on the material geometry, make a list of SOURCE cells that
  // are in each material and collect their material volume fractions
  // and material centroids

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < nsrccells; c++) {
    std::vector<Wonton::Point<2>> ccoords;
    sourceMeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = sourceMeshWrapper.cell_volume(c);

    Wonton::Point<2> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<2>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<2>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_src[m].push_back(c);
          matvf_src[m].push_back(xmoments[0]/cellvol);

          Wonton::Point<2> mcen(xmoments[1]/xmoments[0],
                                 xmoments[2]/xmoments[0]);
          matcen_src[m].push_back(mcen);
        }
      }
    }
  }


  //-------------------------------------------------------------------
  // Now add the material and material cells to the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.add_material(matnames[m], matcells_src[m]);

  // Create multi-material variables to store the volume fractions and
  // centroids for each material in the cells
  for (int m = 0; m < nmats; m++) {
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &(matvf_src[m][0]));
    sourceStateWrapper.mat_add_celldata("mat_centroids", m, &(matcen_src[m][0]));
  }

  // Also assign a different constant density value for each material
  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.mat_add_celldata("density", m, matrho[m]);
  
  //-------------------------------------------------------------------
  // Sanity check - do we get the right volumes and masses for
  // materials from the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    sourceStateWrapper.mat_get_cells(m, &matcells);
    double const *vf;
    double const *rho;
    sourceStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    sourceStateWrapper.mat_get_celldata("density", m, &rho);

    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*sourceMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
    }

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);
  }


  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  std::vector<std::string> remap_fields = {"density"};

  //-------------------------------------------------------------------
  // Add the materials and fields to the target mesh.
  //-------------------------------------------------------------------

  std::vector<int> dummymatcells;
  targetStateWrapper.add_material("mat0", dummymatcells);
  targetStateWrapper.add_material("mat1", dummymatcells);
  targetStateWrapper.add_material("mat2", dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Wonton::Point<2>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  //  Run the remap driver using MOF as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2
  //-------------------------------------------------------------------

  Portage::UberDriver<2,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  // Set all_convex to false, this is a requirement for MM swept-face remap
  d.set_interface_reconstructor_options(false); 

  d.compute_interpolation_weights<Portage::SearchSweptFace, Portage::IntersectSweptFace2D>();

  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_1stOrder>("density", "density", 0.0, dblmax);


  //-------------------------------------------------------------------
  // Based on the material geometry, make a list of the TARGET cells
  // that are in each material
  //-------------------------------------------------------------------

  std::vector<int> matcells_trg[nmats];
  std::vector<double> matvf_trg[nmats];

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < ntrgcells; c++) {
    std::vector<Wonton::Point<2>> ccoords;
    targetMeshWrapper.cell_get_coordinates(c, &ccoords);

    Wonton::Point<2> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<2>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<2>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_trg[m].push_back(c);
        }
      }
    }
  }

  //------------------------------------------------------------------
  // Exact analytic volume fractions on the target mesh
  //------------------------------------------------------------------
  //  ___________
  // |   |0: |  2|
  // |2__|5:_|8__|
  // |   | :.|...|
  // |1__\4:_|7_1|
  // |   / : |   |
  // |0__|3:_|6__|

  // Difference between the source volume fraction of 1/2 and a target volume
  // fraction in a multi-material cell for this particular case. For target,
  // we subtract the area of the triangle with base=DX and height DY (cell
  // height) from the volume of one material in the numerator and from
  // the cell volume in the denominator. 
  double DY = 1./3.;
  double CELL_VOL = DY*DY;
  double TRIANGLE_VOL = DX*DY/2.;
  double vf_diff = 1./2.-(CELL_VOL/2.-TRIANGLE_VOL)/(CELL_VOL-TRIANGLE_VOL);
  double partial_vf1 = (DY*DY)/((DY+DX)*DY); 
  double partial_vf2 = partial_vf1/2.; 
 
  matvf_trg[0].push_back(1.0);              // cell 0
  matvf_trg[0].push_back(1.0);              // cell 1
  matvf_trg[0].push_back(1.0);              // cell 2
  matvf_trg[0].push_back(0.5 - vf_diff);    // cell 3
  matvf_trg[0].push_back(0.5 - vf_diff);    // cell 4
  matvf_trg[0].push_back(0.5);              // cell 5

  matvf_trg[1].push_back(0.5 + vf_diff);    // cell 3   
  matvf_trg[1].push_back((0.5+vf_diff)/2.); // cell 4
  matvf_trg[1].push_back(partial_vf1);      // cell 6
  matvf_trg[1].push_back(partial_vf2);      // cell 7

  matvf_trg[2].push_back((0.5+vf_diff)/2.); // cell 4
  matvf_trg[2].push_back(0.5);              // cell 5
  matvf_trg[2].push_back(partial_vf2);      // cell 7
  matvf_trg[2].push_back(partial_vf1);      // cell 8

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  // First we check that we got 'nmats' materials in the target state
  //-------------------------------------------------------------------

  ASSERT_EQ(nmats, targetStateWrapper.num_materials());
  for (int m = 0; m < nmats; m++)
  ASSERT_EQ(matnames[m], targetStateWrapper.material_name(m));

  std::cerr << "\n\n\n";
  std::cerr << " Number of materials in target mesh: " << nmats << "\n";
  std::cerr << " Material names: ";
  for (int m = 0; m < nmats; m++) std::cerr << " " << matnames[m];
  std::cerr << "\n\n";

  // We compare the material sets calculated by the remapper to
  // the ones calculated independently above

  std::vector<int> matcells_remap[nmats];
  for (int m = 0; m < nmats; m++) {
    targetStateWrapper.mat_get_cells(m, &matcells_remap[m]);
    int nmatcells = matcells_remap[m].size();

    ASSERT_EQ(matcells_trg[m].size(), unsigned(nmatcells));

    std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    std::sort(matcells_trg[m].begin(), matcells_trg[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_trg[m][ic]);
  }
  

  // Then check volume fracs with independently calculated vals
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matvf_trg[m][ic], matvf_remap[ic], 1.0e-10);
   
    for (int ic = 0; ic < nmatcells; ic++)
      std::cout<<" m "<<m<<" cell "<<matcells_remap[m][ic]
               <<" vf: "<<matvf_trg[m][ic]<<" "<<matvf_remap[ic]<<std::endl;

    double const *density_remap;
    targetStateWrapper.mat_get_celldata("density", m, &density_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matrho[m], density_remap[ic], 1.0e-10);
   
    for (int ic = 0; ic < nmatcells; ic++)
      std::cout<<" m "<<m<<" cell "<<matcells_remap[m][ic]
               <<" density: "<<matrho[m]<<" "<<density_remap[ic]<<std::endl;

#ifdef DEBUG
    std::cerr << "Number of cells in material " << m << " is " << nmatcells << "\n";
    std::cerr << "Material " << m << " cell ids, volume fractions:"
              << "\n";
    for (int ic = 0; ic < nmatcells; ic++)
      std::cerr <<
          "  ID = " << std::setw(2) << matcells_remap[m][ic] <<
          "  Vol.Frac. = " << std::setw(6) << matvf_remap[ic] <<
          "  Density = " << std::setw(4) << density_remap[ic] << "\n";
    std::cerr << "\n";
#endif
  }

  // Also check total material volume and mass on the target side
  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);
    double *vf, *rho;
    Wonton::Point<2> *cen;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &cen);
    targetStateWrapper.mat_get_celldata("density", m, &rho);

    Wonton::Point<2> totcen;
    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*targetMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
      totcen += rho[ic]*cen[ic]*cellvol;
    }
    totcen /= mass;

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);

#ifdef DEBUG
    std::cerr << "\nMaterial " << m << "\n";
    std::cerr << " Expected volume " << std::setw(3) <<
        matvol[m] << "    Computed volume " << std::setw(3) << volume << "\n";
    std::cerr << " Expected mass " << std::setw(3) <<
        matmass[m] << "   Computed mass " << std::setw(3) << mass << "\n";
    std::cerr << "\n\n";
#endif
  }

} //ThreeMat2D_1stOrder


TEST(UberDriverSwept, ThreeMat2D_2ndOrder_MisMatch) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 3, 3);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 3, 3);

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  // move few(4) nodes (along x=1,y=0:3 w.r.t grid) within material 0 of the 
  // target mesh in x direction by DX
  double DX = 0.05;
  std::array<double, 2> pnt;
  
  for (int i = 0; i < 4; i++) {
   targetMesh->node_get_coordinates(4+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(4+i, pnt.data());
  }

  //move the right boundary nodes for the target, so 
  //that the last layer of target cells are partially
  //covered by the source mesh
  for (int i = 0; i < 4; i++) {
   targetMesh->node_get_coordinates(12+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(12+i, pnt.data());
  }

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  // The material geometry in the overall domain will look like this
  // and we will put down a rectangular mesh that has multiple cells
  // in each direction on this domain so that we get some pure and
  // some mixed cells
  //
  // Note that only MOF type algorithms or VOF algorithms with the
  // material ordering 0,1,2 will get the T-junction geometry right
  //
  //    0,1           0.5,1         1,1
  //     *-------------:------------*
  //     |             :            |
  //     |             :        2   |
  //     |             :     mat2   |
  //     |             :            |
  //     |             :            |
  //     |     0       +............|1,0.5
  //     |   mat0      :            |
  //     |             :            |
  //     |             :     mat1   |
  //     |             :        1   |
  //     |             :            |
  //     *-------------:------------*
  //    0,0           0.5,0         1,0

  constexpr int nmats = 3;
  std::string matnames[nmats] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  Wonton::Point<2> matlo[nmats], mathi[nmats];
  matlo[0] = Wonton::Point<2>(0.0, 0.0);
  mathi[0] = Wonton::Point<2>(0.5, 1.0);
  matlo[1] = Wonton::Point<2>(0.5, 0.0);
  mathi[1] = Wonton::Point<2>(1.0, 0.5);
  matlo[2] = Wonton::Point<2>(0.5, 0.5);
  mathi[2] = Wonton::Point<2>(1.0, 1.0);

  double matvol[nmats] = {0.5, 0.25, 0.25};
  double matmass[nmats] = {0.375, 1.0, 3.3750};

  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Wonton::Point<2>> matcen_src[nmats];
  std::vector<double> matrho_src[nmats];

  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  //
  // Based on the material geometry, make a list of SOURCE cells that
  // are in each material and collect their material volume fractions
  // and material centroids

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < nsrccells; c++) {
    std::vector<Wonton::Point<2>> ccoords;
    sourceMeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = sourceMeshWrapper.cell_volume(c);

    Wonton::Point<2> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<2>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<2>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_src[m].push_back(c);
          matvf_src[m].push_back(xmoments[0]/cellvol);

          Wonton::Point<2> mcen(xmoments[1]/xmoments[0],
                                xmoments[2]/xmoments[0]);
          matcen_src[m].push_back(mcen);
          matrho_src[m].push_back((m+1)*(m+1)*(mcen[0]+mcen[1]));
        }
      }
    }
  }

  //-------------------------------------------------------------------
  // Now add the material and material cells to the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.add_material(matnames[m], matcells_src[m]);

  // Create multi-material variables to store the volume fractions and
  // centroids for each material in the cells
  for (int m = 0; m < nmats; m++) {
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &(matvf_src[m][0]));
    sourceStateWrapper.mat_add_celldata("mat_centroids", m, &(matcen_src[m][0]));
  }

  // Also assign a different constant density value for each material
  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.mat_add_celldata("density", m, &(matrho_src[m][0]));
  
  //-------------------------------------------------------------------
  // Sanity check - do we get the right volumes and masses for
  // materials from the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    sourceStateWrapper.mat_get_cells(m, &matcells);
    double const *vf;
    double const *rho;
    sourceStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    sourceStateWrapper.mat_get_celldata("density", m, &rho);

    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*sourceMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
    }

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);
  }


  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  std::vector<std::string> remap_fields = {"density"};

  //-------------------------------------------------------------------
  // Add the materials and fields to the target mesh.
  //-------------------------------------------------------------------

  std::vector<int> dummymatcells;
  targetStateWrapper.add_material("mat0", dummymatcells);
  targetStateWrapper.add_material("mat1", dummymatcells);
  targetStateWrapper.add_material("mat2", dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Wonton::Point<2>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  //  Run the remap driver using MOF as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2
  //-------------------------------------------------------------------

  Portage::UberDriver<2,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  // Set all_convex to false, this is a requirement for MM swept-face remap
  d.set_interface_reconstructor_options(false); 

  d.compute_interpolation_weights<Portage::SearchSweptFace, Portage::IntersectSweptFace2D>();

  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_2ndOrder>("density", "density", 0.0, dblmax);

  //-------------------------------------------------------------------
  // Based on the material geometry, make a list of the TARGET cells
  // that are in each material and collect their material volume
  // fractions and material centroids
  //-------------------------------------------------------------------

  std::vector<int> matcells_trg[nmats];
  std::vector<double> matvf_trg[nmats];
  std::vector<Wonton::Point<2>> matcen_trg[nmats];
  std::vector<double> matrho_trg[nmats];

  int ntrgcells = targetMeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                                 Portage::Entity_type::ALL);
  for (int c = 0; c < ntrgcells; c++) {
    std::vector<Wonton::Point<2>> ccoords;
    targetMeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = targetMeshWrapper.cell_volume(c);

    Wonton::Point<2> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<2>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<2>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_trg[m].push_back(c);
          matvf_trg[m].push_back(xmoments[0]/cellvol);

          Wonton::Point<2> mcen(xmoments[1]/xmoments[0],
                                 xmoments[2]/xmoments[0]);
          matcen_trg[m].push_back(mcen);
          matrho_trg[m].push_back((m+1)*(m+1)*(mcen[0]+mcen[1]));
        }
      }
    }
  }

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  // First we check that we got 'nmats' materials in the target state
  //-------------------------------------------------------------------

  ASSERT_EQ(nmats, targetStateWrapper.num_materials());
  for (int m = 0; m < nmats; m++)
  ASSERT_EQ(matnames[m], targetStateWrapper.material_name(m));

  std::cerr << "\n\n\n";
  std::cerr << " Number of materials in target mesh: " << nmats << "\n";
  std::cerr << " Material names: ";
  for (int m = 0; m < nmats; m++) std::cerr << " " << matnames[m];
  std::cerr << "\n\n";

  // We compare the material sets calculated by the remapper to
  // the ones calculated independently above

  std::vector<int> matcells_remap[nmats];
  for (int m = 0; m < nmats; m++) {
    targetStateWrapper.mat_get_cells(m, &matcells_remap[m]);
    int nmatcells = matcells_remap[m].size();

    ASSERT_EQ(matcells_trg[m].size(), unsigned(nmatcells));

    std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    std::sort(matcells_trg[m].begin(), matcells_trg[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_trg[m][ic]);
  }
  

  // Then check volume fracs with independently calculated vals
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    for (int ic = 0; ic < nmatcells; ic++)
     ASSERT_NEAR(matvf_trg[m][ic], matvf_remap[ic], 1.0e-10);

    double const *density_remap;
    targetStateWrapper.mat_get_celldata("density", m, &density_remap);

    //The second order interpolation is not linearity preserving for the
    //swept face based remap due to significant deviation in material
    //centroids after remap. Hence, a very loose tolerance is used for 
    // checks on cell-wise density values. 
    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matrho_trg[m][ic], density_remap[ic], 1.0e-2);

    for (int ic = 0; ic < nmatcells; ic++)
      std::cout<<" m "<<m<<" cell "<<matcells_remap[m][ic]
               <<" density: "<<matrho_trg[m][ic]<<" "<<density_remap[ic]<<std::endl;

#ifdef DEBUG
    std::cerr << "Number of cells in material " << m << " is " << nmatcells << "\n";
    std::cerr << "Material " << m << " cell ids, volume fractions:"
              << "\n";
    for (int ic = 0; ic < nmatcells; ic++)
      std::cerr <<
          "  ID = " << std::setw(2) << matcells_remap[m][ic] <<
          "  Vol.Frac. = " << std::setw(6) << matvf_remap[ic] <<
          "  Density = " << std::setw(4) << density_remap[ic] << "\n";
    std::cerr << "\n";
#endif
  }

  // Also check total material volume and mass on the target side
  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);
    double *vf, *rho;
    Wonton::Point<2> *cen;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &cen);
    targetStateWrapper.mat_get_celldata("density", m, &rho);

    Wonton::Point<2> totcen;
    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*targetMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
      totcen += rho[ic]*cen[ic]*cellvol;
    }
    totcen /= mass;

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);

#ifdef DEBUG
    std::cerr << "\nMaterial " << m << "\n";
    std::cerr << " Expected volume " << std::setw(3) <<
        matvol[m] << "    Computed volume " << std::setw(3) << volume << "\n";
    std::cerr << " Expected mass " << std::setw(3) <<
        matmass[m] << "   Computed mass " << std::setw(3) << mass << "\n";
    std::cerr << "\n\n";
#endif
  }

} //ThreeMat2D_2ndOrder


TEST(UberDriverSwept, ThreeMat3D_1stOrder_MisMatch) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  // move few(8) nodes (along x=1, y=1,2, z=0:3 w.r.t grid) within 
  // material 0 of the target mesh in x direction by DX
  double DX = 0.01;
  std::array<double, 3> pnt;

  for (int i = 0; i < 8; i++) {
   targetMesh->node_get_coordinates(20+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(20+i, pnt.data());
  }

  //move the right boundary nodes for the target, so 
  //that the last layer of target cells are partially
  //covered by the source mesh
  for (int i = 0; i < 16; i++) {
   targetMesh->node_get_coordinates(48+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(48+i, pnt.data());
  }

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  // The material geometry in the overall domain will look like this
  // and we will put down a rectangular mesh that has multiple cells
  // in each direction on this domain so that we get some pure and
  // some mixed cells
  //
  // Note that only MOF type algorithms or VOF algorithms with the
  // material ordering 0,1,2 will get the T-junction geometry right
  //
  //    0,1           0.5,1         1,1
  //     *-------------:------------*
  //     |             :            |
  //     |             :        2   |
  //     |             :     mat2   |
  //     |             :            |
  //     |             :            |
  //     |     0       +............|1,0.5
  //     |   mat0      :            |
  //     |             :            |
  //     |             :     mat1   |
  //     |             :        1   |
  //     |             :            |
  //     *-------------:------------*
  //    0,0           0.5,0         1,0

  constexpr int nmats = 3;
  std::string matnames[nmats] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  Wonton::Point<3> matlo[nmats], mathi[nmats];
  matlo[0] = Wonton::Point<3>(0.0, 0.0, 0.0);
  mathi[0] = Wonton::Point<3>(0.5, 1.0, 1.0);
  matlo[1] = Wonton::Point<3>(0.5, 0.0, 0.0);
  mathi[1] = Wonton::Point<3>(1.0, 0.5, 1.0);
  matlo[2] = Wonton::Point<3>(0.5, 0.5, 0.0);
  mathi[2] = Wonton::Point<3>(1.0, 1.0, 1.0);

  double matrho[nmats] = {0.1, 10.0, 100.0};  // material density
  double matvol[nmats] = {0.5, 0.25, 0.25};
  double matmass[nmats] = {0.05, 2.5, 25.0};

  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Wonton::Point<3>> matcen_src[nmats];


  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  //
  // Based on the material geometry, make a list of SOURCE cells that
  // are in each material and collect their material volume fractions
  // and material centroids

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < nsrccells; c++) {
    std::vector<Wonton::Point<3>> ccoords;
    sourceMeshWrapper.cell_get_coordinates(c, &ccoords);
    
    double cellvol = sourceMeshWrapper.cell_volume(c);

    Wonton::Point<3> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<3>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<3>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_src[m].push_back(c);
          matvf_src[m].push_back(xmoments[0]/cellvol);

          Wonton::Point<3> mcen(xmoments[1]/xmoments[0],
                                 xmoments[2]/xmoments[0],
                                 xmoments[3]/xmoments[0]);
          matcen_src[m].push_back(mcen);
        }
      }
    }
  }


  //-------------------------------------------------------------------
  // Now add the material and material cells to the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.add_material(matnames[m], matcells_src[m]);

  // Create multi-material variables to store the volume fractions and
  // centroids for each material in the cells
  for (int m = 0; m < nmats; m++) {
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &(matvf_src[m][0]));
    sourceStateWrapper.mat_add_celldata("mat_centroids", m, &(matcen_src[m][0]));
  }

  // Also assign a different constant density value for each material
  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.mat_add_celldata("density", m, matrho[m]);
  
  //-------------------------------------------------------------------
  // Sanity check - do we get the right volumes and masses for
  // materials from the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    sourceStateWrapper.mat_get_cells(m, &matcells);
    double const *vf;
    double const *rho;
    sourceStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    sourceStateWrapper.mat_get_celldata("density", m, &rho);

    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*sourceMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
    }

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);
  }


  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  std::vector<std::string> remap_fields = {"density"};

  //-------------------------------------------------------------------
  // Add the materials and fields to the target mesh.
  //-------------------------------------------------------------------

  std::vector<int> dummymatcells;
  targetStateWrapper.add_material("mat0", dummymatcells);
  targetStateWrapper.add_material("mat1", dummymatcells);
  targetStateWrapper.add_material("mat2", dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Wonton::Point<3>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  //  Run the remap driver using MOF as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2
  //-------------------------------------------------------------------

  Portage::UberDriver<3,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  // Set all_convex to false, this is a requirement for MM swept-face remap
  d.set_interface_reconstructor_options(false); 

  d.compute_interpolation_weights<Portage::SearchSweptFace, Portage::IntersectSweptFace3D>();

  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_1stOrder>("density", "density", 0.0, dblmax);


  //-------------------------------------------------------------------
  // Based on the material geometry, make a list of the TARGET cells
  // that are in each material
  //-------------------------------------------------------------------

  std::vector<int> matcells_trg[nmats];
  std::vector<double> matvf_trg[nmats];

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < ntrgcells; c++) {
    std::vector<Wonton::Point<3>> ccoords;
    targetMeshWrapper.cell_get_coordinates(c, &ccoords);

    Wonton::Point<3> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<3>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<3>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_trg[m].push_back(c);
        }
      }
    }
  }

  //------------------------------------------------------------------
  // Exact analytic volume fractions on the target mesh
  //------------------------------------------------------------------
  // material layout
  // ---------------
  // |   \  :  | 2  |
  // ---------------
  // |  0 | : .|....|
  // ---------------
  // |   /  :  | 1  |
  // ---------------
  //
  // cells: index z = 0 
  /// --------------
  // | 6  \  15 | 24 |
  // ---------------
  // | 3   | 12 | 21 |
  // ---------------
  // | 0  /  9  | 18 |
  // ----------------
  // index z = 1
  /// ----------------
  // | 7 \  16 | 25 |
  // -----------------
  // | 4  | 13 | 22 |
  // -----------------
  // | 1 /  10 | 19 |
  // -----------------
  // index z = 2
  /// ----------------
  // | 8 \  17 | 26 |
  // -----------------
  // | 5  | 14 | 23 |
  // -----------------
  // | 2 /  11 | 20 |
  // -----------------

  // Difference between the source volume fraction of 1/2 and a target volume
  // fraction in a multi-material cell for this particular case. For target,
  // we subtract the area of the triangle with base=DX and height DY (cell
  // height) from the volume of one material in the numerator and from
  // the cell volume in the denominator. 
  double DH = 1./3.;
  double CVOL = DH*DH*DH;
  double CVOL1 = DX*DH*DH/2.;
  double CVOL2 = 2*CVOL1;  
  double vf_diff1 = 1./2.-(CVOL/2.-CVOL1)/(CVOL-CVOL1);
  double vf_diff2 = 1./2.-(CVOL/2.-CVOL2)/(CVOL-CVOL2);
  double partial_vf1 = (DH*DH*DH)/((DH+DX)*DH*DH); 
  double partial_vf2 = partial_vf1/2.; 

  //mat 0: cells 0, 1, 2, 3, 4, 5, 6, 7, 8 
  //             9, 10, 11, 12, 13, 14, 15, 16, 17
  matvf_trg[0].push_back(1.0);             
  matvf_trg[0].push_back(1.0);             
  matvf_trg[0].push_back(1.0);            
  matvf_trg[0].push_back(1.0);             
  matvf_trg[0].push_back(1.0);             
  matvf_trg[0].push_back(1.0);            
  matvf_trg[0].push_back(1.0);             
  matvf_trg[0].push_back(1.0);             
  matvf_trg[0].push_back(1.0);            
  matvf_trg[0].push_back(0.5 - vf_diff1);   
  matvf_trg[0].push_back(0.5 - vf_diff1);   
  matvf_trg[0].push_back(0.5 - vf_diff1);   
  matvf_trg[0].push_back(0.5 - vf_diff2);   
  matvf_trg[0].push_back(0.5 - vf_diff2);   
  matvf_trg[0].push_back(0.5 - vf_diff2);   
  matvf_trg[0].push_back(0.5 - vf_diff1);   
  matvf_trg[0].push_back(0.5 - vf_diff1);   
  matvf_trg[0].push_back(0.5 - vf_diff1);   

  //mat 0: cells 9, 10, 11, 12, 13, 14, 
  //             18, 19, 20, 21, 22, 23
  matvf_trg[1].push_back(0.5 + vf_diff1);      
  matvf_trg[1].push_back(0.5 + vf_diff1);      
  matvf_trg[1].push_back(0.5 + vf_diff1);      
  matvf_trg[1].push_back((0.5+vf_diff2)/2.);
  matvf_trg[1].push_back((0.5+vf_diff2)/2.);
  matvf_trg[1].push_back((0.5+vf_diff2)/2.);
  matvf_trg[1].push_back(partial_vf1);             
  matvf_trg[1].push_back(partial_vf1);             
  matvf_trg[1].push_back(partial_vf1);             
  matvf_trg[1].push_back(partial_vf2);             
  matvf_trg[1].push_back(partial_vf2);             
  matvf_trg[1].push_back(partial_vf2);             

  //mat 0: cells 12, 13, 14, 15, 16, 17, 
  //             21, 22, 23, 24, 25, 26
  matvf_trg[2].push_back((0.5+vf_diff2)/2.);
  matvf_trg[2].push_back((0.5+vf_diff2)/2.);
  matvf_trg[2].push_back((0.5+vf_diff2)/2.);
  matvf_trg[2].push_back(0.5+vf_diff1);             
  matvf_trg[2].push_back(0.5+vf_diff1);             
  matvf_trg[2].push_back(0.5+vf_diff1);             
  matvf_trg[2].push_back(partial_vf2);             
  matvf_trg[2].push_back(partial_vf2);             
  matvf_trg[2].push_back(partial_vf2);             
  matvf_trg[2].push_back(partial_vf1);             
  matvf_trg[2].push_back(partial_vf1);             
  matvf_trg[2].push_back(partial_vf1);             

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  // First we check that we got 'nmats' materials in the target state
  //-------------------------------------------------------------------

  ASSERT_EQ(nmats, targetStateWrapper.num_materials());
  for (int m = 0; m < nmats; m++)
  ASSERT_EQ(matnames[m], targetStateWrapper.material_name(m));

  std::cerr << "\n\n\n";
  std::cerr << " Number of materials in target mesh: " << nmats << "\n";
  std::cerr << " Material names: ";
  for (int m = 0; m < nmats; m++) std::cerr << " " << matnames[m];
  std::cerr << "\n\n";

  // We compare the material sets calculated by the remapper to
  // the ones calculated independently above

  std::vector<int> matcells_remap[nmats];
  for (int m = 0; m < nmats; m++) {
    targetStateWrapper.mat_get_cells(m, &matcells_remap[m]);
    int nmatcells = matcells_remap[m].size();

    ASSERT_EQ(matcells_trg[m].size(), unsigned(nmatcells));

    std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    std::sort(matcells_trg[m].begin(), matcells_trg[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_trg[m][ic]);
  }
  

  // Then check volume fracs with independently calculated vals
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    for (int ic = 0; ic < nmatcells; ic++)
     ASSERT_NEAR(matvf_trg[m][ic], matvf_remap[ic], 1.0e-10);

    double const *density_remap;
    targetStateWrapper.mat_get_celldata("density", m, &density_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matrho[m], density_remap[ic], 1.0e-10);

#ifdef DEBUG
    std::cerr << "Number of cells in material " << m << " is " << nmatcells << "\n";
    std::cerr << "Material " << m << " cell ids, volume fractions:"
              << "\n";
    for (int ic = 0; ic < nmatcells; ic++)
      std::cerr <<
          "  ID = " << std::setw(2) << matcells_remap[m][ic] <<
          "  Vol.Frac. = " << std::setw(6) << matvf_remap[ic] <<
          "  Density = " << std::setw(4) << density_remap[ic] << "\n";
    std::cerr << "\n";
#endif
  }

  // Also check total material volume and mass on the target side
  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);
    double *vf, *rho;
    Wonton::Point<3> *cen;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &cen);
    targetStateWrapper.mat_get_celldata("density", m, &rho);

    Wonton::Point<3> totcen;
    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*targetMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
      totcen += rho[ic]*cen[ic]*cellvol;
    }
    totcen /= mass;

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);

#ifdef DEBUG
    std::cerr << "\nMaterial " << m << "\n";
    std::cerr << " Expected volume " << std::setw(3) <<
        matvol[m] << "    Computed volume " << std::setw(3) << volume << "\n";
    std::cerr << " Expected mass " << std::setw(3) <<
        matmass[m] << "   Computed mass " << std::setw(3) << mass << "\n";
    std::cerr << "\n\n";
#endif
  }

} //ThreeMat3D_1stOrder


TEST(UberDriverSwept, ThreeMat3D_2ndOrder_MisMatch) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  // move few(16) nodes (along x=1, y=0:3, z=0:3 w.r.t grid) within 
  // material 0 of the target mesh in x direction by DX
  double DX = 0.01;
  std::array<double, 3> pnt;

  for (int i = 0; i < 16; i++) {
   targetMesh->node_get_coordinates(16+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(16+i, pnt.data());
  }

  //move the right boundary nodes for the target, so 
  //that the last layer of target cells are partially
  //covered by the source mesh
  for (int i = 0; i < 16; i++) {
   targetMesh->node_get_coordinates(48+i, &pnt);
   pnt[0] += DX;
   targetMesh->node_set_coordinates(48+i, pnt.data());
  }

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  // The material geometry in the overall domain will look like this
  // and we will put down a rectangular mesh that has multiple cells
  // in each direction on this domain so that we get some pure and
  // some mixed cells
  //
  // Note that only MOF type algorithms or VOF algorithms with the
  // material ordering 0,1,2 will get the T-junction geometry right
  //
  //    0,1           0.5,1         1,1
  //     *-------------:------------*
  //     |             :            |
  //     |             :        2   |
  //     |             :     mat2   |
  //     |             :            |
  //     |             :            |
  //     |     0       +............|1,0.5
  //     |   mat0      :            |
  //     |             :            |
  //     |             :     mat1   |
  //     |             :        1   |
  //     |             :            |
  //     *-------------:------------*
  //    0,0           0.5,0         1,0

  constexpr int nmats = 3;
  std::string matnames[nmats] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  Wonton::Point<3> matlo[nmats], mathi[nmats];
  matlo[0] = Wonton::Point<3>(0.0, 0.0, 0.0);
  mathi[0] = Wonton::Point<3>(0.5, 1.0, 1.0);
  matlo[1] = Wonton::Point<3>(0.5, 0.0, 0.0);
  mathi[1] = Wonton::Point<3>(1.0, 0.5, 1.0);
  matlo[2] = Wonton::Point<3>(0.5, 0.5, 0.0);
  mathi[2] = Wonton::Point<3>(1.0, 1.0, 1.0);

  double matvol[nmats] = {0.5, 0.25, 0.25};
  double matmass[nmats] = {0.625, 1.5, 4.5};

  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Wonton::Point<3>> matcen_src[nmats];
  std::vector<double> matrho_src[nmats];

  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  //
  // Based on the material geometry, make a list of SOURCE cells that
  // are in each material and collect their material volume fractions
  // and material centroids

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < nsrccells; c++) {
    std::vector<Wonton::Point<3>> ccoords;
    sourceMeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = sourceMeshWrapper.cell_volume(c);

    Wonton::Point<3> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<3>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<3>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_src[m].push_back(c);
          matvf_src[m].push_back(xmoments[0]/cellvol);

          Wonton::Point<3> mcen(xmoments[1]/xmoments[0],
                                xmoments[2]/xmoments[0],
                                xmoments[3]/xmoments[0]);
          matcen_src[m].push_back(mcen);
          matrho_src[m].push_back((m+1)*(m+1)*(mcen[0]+mcen[1]+mcen[2]));
        }
      }
    }
  }

  //-------------------------------------------------------------------
  // Now add the material and material cells to the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.add_material(matnames[m], matcells_src[m]);

  // Create multi-material variables to store the volume fractions and
  // centroids for each material in the cells
  for (int m = 0; m < nmats; m++) {
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &(matvf_src[m][0]));
    sourceStateWrapper.mat_add_celldata("mat_centroids", m, &(matcen_src[m][0]));
  }

  // Also assign a different constant density value for each material
  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.mat_add_celldata("density", m, &(matrho_src[m][0]));
  
  //-------------------------------------------------------------------
  // Sanity check - do we get the right volumes and masses for
  // materials from the source state
  //-------------------------------------------------------------------

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    sourceStateWrapper.mat_get_cells(m, &matcells);
    double const *vf;
    double const *rho;
    sourceStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    sourceStateWrapper.mat_get_celldata("density", m, &rho);

    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*sourceMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
    }

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);
  }


  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  std::vector<std::string> remap_fields = {"density"};

  //-------------------------------------------------------------------
  // Add the materials and fields to the target mesh.
  //-------------------------------------------------------------------

  std::vector<int> dummymatcells;
  targetStateWrapper.add_material("mat0", dummymatcells);
  targetStateWrapper.add_material("mat1", dummymatcells);
  targetStateWrapper.add_material("mat2", dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Wonton::Point<3>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  //  Run the remap driver using MOF as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2
  //-------------------------------------------------------------------

  Portage::UberDriver<3,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  // Set all_convex to false, this is a requirement for MM swept-face remap
  d.set_interface_reconstructor_options(false); 

  d.compute_interpolation_weights<Portage::SearchSweptFace, Portage::IntersectSweptFace3D>();

  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_2ndOrder>("density", "density", 0.0, dblmax);

  //-------------------------------------------------------------------
  // Based on the material geometry, make a list of the TARGET cells
  // that are in each material and collect their material volume
  // fractions and material centroids
  //-------------------------------------------------------------------

  std::vector<int> matcells_trg[nmats];
  std::vector<double> matvf_trg[nmats];
  std::vector<Wonton::Point<3>> matcen_trg[nmats];
  std::vector<double> matrho_trg[nmats];

  int ntrgcells = targetMeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                                 Portage::Entity_type::ALL);
  for (int c = 0; c < ntrgcells; c++) {
    std::vector<Wonton::Point<3>> ccoords;
    targetMeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = targetMeshWrapper.cell_volume(c);

    Wonton::Point<3> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<3>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < nmats; m++) {
      if (BOX_INTERSECT::intersect_boxes<3>(matlo[m], mathi[m],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          matcells_trg[m].push_back(c);
          matvf_trg[m].push_back(xmoments[0]/cellvol);

          Wonton::Point<3> mcen(xmoments[1]/xmoments[0],
                                 xmoments[2]/xmoments[0],
                                 xmoments[3]/xmoments[0]);
          matcen_trg[m].push_back(mcen);
          matrho_trg[m].push_back((m+1)*(m+1)*(mcen[0]+mcen[1]+mcen[2]));
        }
      }
    }
  }

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  // First we check that we got 'nmats' materials in the target state
  //-------------------------------------------------------------------

  ASSERT_EQ(nmats, targetStateWrapper.num_materials());
  for (int m = 0; m < nmats; m++)
  ASSERT_EQ(matnames[m], targetStateWrapper.material_name(m));

  std::cerr << "\n\n\n";
  std::cerr << " Number of materials in target mesh: " << nmats << "\n";
  std::cerr << " Material names: ";
  for (int m = 0; m < nmats; m++) std::cerr << " " << matnames[m];
  std::cerr << "\n\n";

  // We compare the material sets calculated by the remapper to
  // the ones calculated independently above

  std::vector<int> matcells_remap[nmats];
  for (int m = 0; m < nmats; m++) {
    targetStateWrapper.mat_get_cells(m, &matcells_remap[m]);
    int nmatcells = matcells_remap[m].size();

    ASSERT_EQ(matcells_trg[m].size(), unsigned(nmatcells));

    std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    std::sort(matcells_trg[m].begin(), matcells_trg[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_trg[m][ic]);
  }
  

  // Then check volume fracs with independently calculated vals
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    for (int ic = 0; ic < nmatcells; ic++)
     ASSERT_NEAR(matvf_trg[m][ic], matvf_remap[ic], 1.0e-10);

    double const *density_remap;
    targetStateWrapper.mat_get_celldata("density", m, &density_remap);

    //The second order interpolation is not linearity preserving for the
    //swept face based remap due to significant deviation in material
    //centroids after remap. Hence, a very loose tolerance is used for 
    // checks on cell-wise density values. 
    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matrho_trg[m][ic], density_remap[ic], 1.0e-2);

#ifdef DEBUG
    std::cerr << "Number of cells in material " << m << " is " << nmatcells << "\n";
    std::cerr << "Material " << m << " cell ids, volume fractions:"
              << "\n";
    for (int ic = 0; ic < nmatcells; ic++)
      std::cerr <<
          "  ID = " << std::setw(2) << matcells_remap[m][ic] <<
          "  Vol.Frac. = " << std::setw(6) << matvf_remap[ic] <<
          "  Density = " << std::setw(4) << density_remap[ic] << "\n";
    std::cerr << "\n";
#endif
  }

  // Also check total material volume and mass on the target side
  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    targetStateWrapper.mat_get_cells(m, &matcells);
    double *vf, *rho;
    Wonton::Point<3> *cen;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &vf);
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &cen);
    targetStateWrapper.mat_get_celldata("density", m, &rho);

    Wonton::Point<3> totcen;
    double volume = 0.0, mass = 0.0;
    int const num_matcells = matcells.size();
    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*targetMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
      totcen += rho[ic]*cen[ic]*cellvol;
    }
    totcen /= mass;

    ASSERT_NEAR(matvol[m], volume, 1.0e-10);
    ASSERT_NEAR(matmass[m], mass, 1.0e-10);

#ifdef DEBUG
    std::cerr << "\nMaterial " << m << "\n";
    std::cerr << " Expected volume " << std::setw(3) <<
        matvol[m] << "    Computed volume " << std::setw(3) << volume << "\n";
    std::cerr << " Expected mass " << std::setw(3) <<
        matmass[m] << "   Computed mass " << std::setw(3) << mass << "\n";
    std::cerr << "\n\n";
#endif
  }

} //ThreeMat3D_2ndOrder

#endif  // ifdef PORTAGE_HAS_TANGRAM
