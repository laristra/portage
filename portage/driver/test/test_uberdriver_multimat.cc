/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifdef HAVE_TANGRAM

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "portage/driver/uberdriver.h"

double TOL = 1e-6;


// Tests for multi-material remap with 1st Order Accurate Remap

// The conceptual material layout that includes a T-junction and
// interfaces that are aligned with the coordinate axes/planes. The
// three materials have very disparate but constant densities. Given
// this setup we know the volume, mass and centroid of each material.

// Knowing the simple material geometry, we can do box intersections
// to compute exact values of material volume fractions and centroids
// in each cell on the source and target meshes. We can then compare
// the geometrically computed values of the volume fractions and
// centroids to the ones we get out of the code. Additionally, we can
// compute the mass of the materials based on the REMAPPED densities
// on the target mesh and compare it to the analytical values


TEST(UberDriver, ThreeMat2D_1stOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 7, 6);

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

  double meshtemp = 55;  // scalar temperature on mesh
  double matrho[nmats] = {0.1, 10.0, 100.0};  // material density
  double matvol[nmats] = {0.5, 0.25, 0.25};
  double matmass[nmats] = {0.05, 2.5, 25.0};
  Wonton::Point<2> matcen[nmats] = {Wonton::Point<2>(0.25,0.5),
                                    Wonton::Point<2>(0.75,0.25),
                                    Wonton::Point<2>(0.75,0.75)};

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
  // Now add temperature field to the mesh
  //-------------------------------------------------------------------

  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "temperature", meshtemp);

  
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
    for (int ic = 0; ic < matcells.size(); ic++) {
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

  std::vector<std::string> remap_fields = {"density", "temperature"};


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

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", 0.0);

  //-------------------------------------------------------------------
  //  Run the remap driver using XMOF-2D as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2
  //-------------------------------------------------------------------

  Portage::UberDriver<2,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Tangram::XMOF2D_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  d.set_num_tols<Portage::Entity_kind::CELL>(default_num_tols);

  d.compute_interpolation_weights<Portage::SearchKDTree, Portage::IntersectR2D>();

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_1stOrder>("density", "density",
  					       0.0, dblmax);
  d.interpolate<double, Portage::Entity_kind::CELL,
  		Portage::Interpolate_1stOrder>("temperature", dblmin, dblmax);


  //-------------------------------------------------------------------
  // Based on the material geometry, make a list of the TARGET cells
  // that are in each material and collect their material volume
  // fractions and material centroids
  //-------------------------------------------------------------------

  std::vector<int> matcells_trg[nmats];
  std::vector<double> matvf_trg[nmats];
  std::vector<Wonton::Point<2>> matcen_trg[nmats];

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
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

    ASSERT_EQ(matcells_trg[m].size(), nmatcells);

    std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    std::sort(matcells_trg[m].begin(), matcells_trg[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_trg[m][ic]);
  }



  // Then check volume fracs and centroids with independently calculated vals
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matvf_trg[m][ic], matvf_remap[ic], 1.0e-10);

    Wonton::Point<2> const *matcen_remap;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      for (int d = 0; d < 2; d++)
        ASSERT_NEAR(matcen_trg[m][ic][d], matcen_remap[ic][d], 1.0e-10);

    double const *density_remap;
    targetStateWrapper.mat_get_celldata("density", m, &density_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matrho[m], density_remap[ic], 1.0e-10);

#ifdef DEBUG
    std::cerr << "Number of cells in material " << m << " is " << nmatcells << "\n";
    std::cerr << "Material " << m << " cell ids, volume fractions, centroids:"
              << "\n";
    for (int ic = 0; ic < nmatcells; ic++)
      std::cerr <<
          "  ID = " << std::setw(2) << matcells_remap[m][ic] <<
          "  Vol.Frac. = " << std::setw(6) << matvf_remap[ic] <<
          "  Centroid = " << std::setw(6) << matcen_remap[ic][0] << std::setw(6) << matcen_remap[ic][1] <<
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
    for (int ic = 0; ic < matcells.size(); ic++) {
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
    std::cerr << " Expected centroid " << std::setw(6) <<
        matcen[m] << "  Computed centroid " << totcen << "\n";
    std::cerr << "\n\n";
#endif
  }

  // Finally check that we got the right target temperature values
  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "temperature",
                                   &targettemp);
  for (int i = 0; i < ntrgcells; i++)
    ASSERT_NEAR(targettemp[i], meshtemp, 1.0e-10);

}  // ThreeMat2D_1stOrder





TEST(UberDriver, ThreeMat3D_1stOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

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

  double meshtemp = 55;  // scalar temperature on mesh
  double matrho[nmats] = {0.1, 10.0, 100.0};  // material density
  double matvol[nmats] = {0.5, 0.25, 0.25};
  double matmass[nmats] = {0.05, 2.5, 25.0};
  Wonton::Point<3> matcen[nmats] = {Wonton::Point<3>(0.25,0.5,0.5),
                                    Wonton::Point<3>(0.75,0.25,0.5),
                                    Wonton::Point<3>(0.75,0.75,0.5)};

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
  // Now add temperature field to the mesh
  //-------------------------------------------------------------------

  sourceStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", meshtemp);
  
  
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
    for (int ic = 0; ic < matcells.size(); ic++) {
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

  std::vector<std::string> remap_fields = {"density", "temperature"};


  //-------------------------------------------------------------------
  // Add the materials and fields to the target mesh.
  //-------------------------------------------------------------------

  std::vector<int> dummymatcells;
  targetStateWrapper.add_material("mat0", dummymatcells);
  targetStateWrapper.add_material("mat1", dummymatcells);
  targetStateWrapper.add_material("mat2", dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Wonton::Point<3>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("rho", 0.0);

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "TEMP", 0.0);


  //-------------------------------------------------------------------
  //  Run the remap driver using MOF-3D as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2 and may not get the right
  //  orientation of interface on boundary cells
  //-------------------------------------------------------------------

  std::vector<Portage::Entity_kind> entity_kinds({Wonton::Entity_kind::CELL});
  std::vector<Portage::Field_type> field_types({Portage::Field_type::MESH_FIELD, Portage::Field_type::MULTIMATERIAL_FIELD});
  
  Wonton::SerialExecutor_type executor;

  Portage::UberDriver<3,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper, &executor);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  d.set_num_tols<Portage::Entity_kind::CELL>(default_num_tols);

  d.compute_interpolation_weights<Portage::SearchKDTree, Portage::IntersectR3D>();

  
  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_1stOrder>("density", "rho", 0.0, dblmax);
  d.interpolate<double, Portage::Entity_kind::CELL,
                  Portage::Interpolate_1stOrder>("temperature", "TEMP", dblmin,
                                                 dblmax);


  //-------------------------------------------------------------------
  // Based on the material geometry, make a list of the TARGET cells
  // that are in each material and collect their material volume
  // fractions and material centroids
  //-------------------------------------------------------------------

  std::vector<int> matcells_trg[nmats];
  std::vector<double> matvf_trg[nmats];
  std::vector<Wonton::Point<3>> matcen_trg[nmats];

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

    ASSERT_EQ(matcells_trg[m].size(), nmatcells);

    std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    std::sort(matcells_trg[m].begin(), matcells_trg[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_trg[m][ic]);
  }



  // Then check volume fracs and centroids with independently calculated vals
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matvf_trg[m][ic], matvf_remap[ic], 1.0e-9);

    Wonton::Point<3> const *matcen_remap;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);

    // MOF cannot match moments and centroids as well as it can volume
    // fractions - so use looser tolerances
    for (int ic = 0; ic < nmatcells; ic++)
      for (int d = 0; d < 3; d++)
        ASSERT_NEAR(matcen_trg[m][ic][d], matcen_remap[ic][d], 1.0e-9);

    double const *density_remap;
    targetStateWrapper.mat_get_celldata("rho", m, &density_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR(matrho[m], density_remap[ic], 1.0e-12);

#ifdef DEBUG
    std::cerr << "Number of cells in material " << m << " is " << nmatcells << "\n";
    std::cerr << "Material " << m << " cell ids, volume fractions, centroids:"
              << "\n";
    for (int ic = 0; ic < nmatcells; ic++)
      std::cerr <<
          "  ID = " << std::setw(2) << matcells_remap[m][ic] <<
          "  Vol.Frac. = " << std::setw(6) << matvf_remap[ic] <<
          "  Centroid = " << std::setw(6) << matcen_remap[ic][0] << " " <<
          std::setw(6) << matcen_remap[ic][1] << " " << std::setw(6) << matcen_remap[ic][2] <<
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
    targetStateWrapper.mat_get_celldata("rho", m, &rho);

    Wonton::Point<3> totcen;
    double volume = 0.0, mass = 0.0;
    for (int ic = 0; ic < matcells.size(); ic++) {
      double cellvol = vf[ic]*targetMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
      totcen += rho[ic]*cen[ic]*cellvol;
    }
    totcen /= mass;

    ASSERT_NEAR(matvol[m], volume, 1.0e-12);

#ifdef DEBUG
    std::cerr << "\nmaterial " << m << "\n";
    std::cerr << " expected volume " << std::setw(3) <<
        matvol[m] << "    computed volume " << std::setw(3) << volume << "\n";
    std::cerr << " expected mass " << std::setw(3) <<
        matmass[m] << "   computed mass " << std::setw(3) << mass << "\n";
    std::cerr << " expected centroid " << std::setw(6) <<
        matcen[m] << "  computed centroid " << totcen << "\n";
    std::cerr << "\n\n";
#endif
  }


  // Finally check that we got the right target temperature values
  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "TEMP",
                                   &targettemp);
  for (int i = 0; i < ntrgcells; i++)
    ASSERT_NEAR(targettemp[i], meshtemp, 1.0e-10);

}  // ThreeMat3D_1stOrder


#endif  // ifdef HAVE_TANGRAM
