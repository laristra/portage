/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "portage/support/portage.h"
#ifdef PORTAGE_HAS_TANGRAM

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector> 

//Wonton includes
#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

//This test is distributed
#include "mpi.h"
#include "gtest/gtest.h"

//Portage includes
#include "portage/driver/uberdriver.h"
#include "portage/driver/mmdriver.h"
#include "portage/intersect/simple_intersect_for_tests.h"

//Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

//Tangram includes
#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/reconstruct/VOF.h"

// Tests for distributed multi-material remap with 1st and 
// 2nd Order Accurate Remap on 4 ranks. 

// This distributed test uses two material configurations for 
// testing purposes. 
//
// MATERIAL CONFIGURATION: NESTED_BOX
// The conceptual material layout is two nested boxes inside an
// outer box domain with unit length. The material interfaces
// align with the coordinate axes/planes. With this layout, all
// ranks in the distributed case will have all three materials. 
//
//
// 2D 
//
//    0,1            0.5,1.          1,1
//     *------------------------------*
//     |            0, mat0           |
//     |                              |    
//     |    *....................*    |
//     |    :       1, mat1      :    |
//     |    :     *.........*    :    |
//     |    :     :         :    :    |
//     |    :     : 2, mat2 :    :    |
//     |    :     :         :    :    | 1,0.5
//     |    :     :         :    :    |
//     |    :     *.........*    :    |
//     |    :                    :    |
//     |    *....................*    |
//     |                              |                     
//     |                              |
//     *------------------------------*
//    0,0            0.5,0           1,0    

// MATERIAL CONFIGURATION: LAYER
// The conceptual material layout is three layered materials
// within the box.  The material interfaces
// align with the coordinate axes/planes. 
//
// 2D 
//
//    0,1      0.3,1     0.65,1      1,1
//     *------------------------------*
//     |        |          |          |
//     |        |          |          |
//     |        |          |          |
//     |        |          |          |
//     | mat 0  |  mat1    |  mat 2   |
//     |        |          |          |
//     |        |          |          |
//     |        |          |          |
//     |        |          |          |
//     |        |          |          |
//     |        |          |          |
//     *------------------------------*
//    0,0      0.3,0     0.65,0      1,0

// The three materials have densities initialized with either constant
// or linear functions. Given these setups we know the volume, mass and
// centroid of each material.

// Knowing the material geometry, we can do box intersections
// to compute exact values of material volume fractions and centroids
// in each cell on the source and target meshes. We can then compare
// the geometrically computed values of the volume fractions and
// centroids to the ones we get out of the code. Additionally, we can
// compute the mass of the materials based on the REMAPPED densities
// on the target mesh and compare it to the analytical values.
// Finally, we make sure that the constant/linear fields are recovered
// on the target mesh side

double TOL = 1e-6;

enum DENSITY_FUNCTION
{ 
  CONSTANT,
  LINEAR
};

enum MAT_CONFIG
{
  NESTED_BOX,
  LAYER
};

template<int dim, MAT_CONFIG mconfig>
struct material_metadata {};

template <>
struct material_metadata<2, MAT_CONFIG::LAYER>
{
  const int nmats = 3; 
  std::string matnames[3] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  std::vector<Portage::Point<2>> matlo[3], mathi[3];

  material_metadata(){
    matlo[0] = {Portage::Point<2>(0.0, 0.0)}; 
    mathi[0] = {Portage::Point<2>(0.3, 1.0)}; 
    matlo[1] = {Portage::Point<2>(0.3, 0.0)};
    mathi[1] = {Portage::Point<2>(0.65, 1.0)};
    matlo[2] = {Portage::Point<2>(0.65, 0.0)}; 
    mathi[2] = {Portage::Point<2>(1.0, 1.0)}; 
  }
  
  int mat_nboxes[3] = {1,1,1}; 

  // Reference Volume 
  std::vector<double> matvol[4] = {{ 0.1714285714286, 0.1551020408163, 0.0}, 
				   { 0.2142857142857, 0.1938775510204, 0.0},
				   { 0.008163265306122, 0.2, 0.2}, 
				   { 0.01020408163265, 0.25, 0.25}}; 

  // Reference Mass
  std::vector<double> matmass_const[4] = {{ 0.01714285714286, 1.551020408163, 0.0}, 
					  { 0.02142857142857, 1.938775510204, 0.0},
				          { 0.0008163265306122, 2.0, 20.0}, 
					  { 0.001020408163265, 2.5, 25.0}}; 


  std::vector<double> matmass_linear[4] = {{ 0.07469387755102, 0.1118950437318, 0.0}, 
					   { 0.1698979591837, 0.209110787172, 0.0},
				           { 0.004723032069971, 0.1521428571429, 0.2221428571429}, 
					   { 0.009548104956268, 0.2794642857143, 0.3669642857143}}; 

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0}; 

  // Compute density
  double compute_density(Portage::Point<2> &coords, int matid, DENSITY_FUNCTION dtype)
  {
    double value = 0;
    switch(dtype)
    {
      case CONSTANT : value = matdensity_const[matid]; break; 
      case LINEAR : value = coords[0]+coords[1]; break; 
    }
    return value; 
  }   
}; //material_metadata

template <>
struct material_metadata<3, MAT_CONFIG::LAYER>
{
  const int nmats = 3; 
  std::string matnames[3] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  std::vector<Portage::Point<3>> matlo[3], mathi[3];

  material_metadata(){
    matlo[0] = {Portage::Point<3>(0.0, 0.0, 0.0)}; 
    mathi[0] = {Portage::Point<3>(0.3, 1.0, 1.0)}; 
    matlo[1] = {Portage::Point<3>(0.3, 0.0, 0.0)};
    mathi[1] = {Portage::Point<3>(0.65, 1.0, 1.0)};
    matlo[2] = {Portage::Point<3>(0.65, 0.0, 0.0)}; 
    mathi[2] = {Portage::Point<3>(1.0, 1.0, 1.0)}; 
  }
  
  int mat_nboxes[3] = {1,1,1}; 

  // Reference Volume 
  std::vector<double> matvol[4] = {{ 0.1714285714286, 0.1551020408163, 0.0}, 
				   { 0.2142857142857, 0.1938775510204, 0.0},
				   { 0.008163265306122, 0.2, 0.2}, 
				   { 0.01020408163265, 0.25, 0.25}}; 

  // Reference Mass
  std::vector<double> matmass_const[4] = {{ 0.01714285714286, 1.551020408163, 0.0}, 
					  { 0.02142857142857, 1.938775510204, 0.0},
				          { 0.0008163265306122, 2.0, 20.0}, 
					  { 0.001020408163265, 2.5, 25.0}}; 


  std::vector<double> matmass_linear[4] = {{ 0.1604081632653, 0.1894460641399, 0.0}, 
					   { 0.2770408163265, 0.3060495626822 , 0.0},
				           { 0.008804664723032, 0.2521428571429, 0.3221428571429}, 
					   { 0.01465014577259, 0.4044642857143, 0.4919642857143}}; 

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0}; 

  // Compute density
  double compute_density(Portage::Point<3> &coords, int matid, DENSITY_FUNCTION dtype)
  {
    double value = 0;
    switch(dtype)
    {
      case CONSTANT : value = matdensity_const[matid]; break; 
      case LINEAR : value = coords[0]+coords[1]+coords[2]; break; 
    }
    return value; 
  }   
}; //material_metadata

template <>
struct material_metadata<2, MAT_CONFIG::NESTED_BOX>
{
  const int nmats = 3;
  std::string matnames[3] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  std::vector<Portage::Point<2>> matlo[3], mathi[3];

  material_metadata(){
    matlo[0] = {Portage::Point<2>(0.0, 0.0), Portage::Point<2>(0.0, 0.2),
                Portage::Point<2>(0.8, 0.2), Portage::Point<2>(0.0, 0.8)};
    mathi[0] = {Portage::Point<2>(1.0, 0.2), Portage::Point<2>(0.2, 0.8),
                Portage::Point<2>(1.0, 0.8), Portage::Point<2>(1.0, 1.0)};

    matlo[1] = {Portage::Point<2>(0.2, 0.2), Portage::Point<2>(0.2, 0.4),
                Portage::Point<2>(0.6, 0.4), Portage::Point<2>(0.2, 0.6)};
    mathi[1] = {Portage::Point<2>(0.8, 0.4), Portage::Point<2>(0.4, 0.6),
                Portage::Point<2>(0.8, 0.6), Portage::Point<2>(0.8, 0.8)};

    matlo[2] = {Portage::Point<2>(0.4, 0.4)};
    mathi[2] = {Portage::Point<2>(0.6, 0.6)};
  }

  int mat_nboxes[3] = {4,4,1};

  // Reference Volume
  std::vector<double> matvol[4] = {{0.2, 0.12, 0.04},
				   {0.2, 0.12, 0.04},
				   {0.2, 0.12, 0.04},
				   {0.2, 0.12, 0.04}}; 

  // Reference Mass
  std::vector<double> matmass_const[4] = {{0.02, 1.2, 4.0},
				   {0.02, 1.2, 4.0},
				   {0.02, 1.2, 4.0},
				   {0.02, 1.2, 4.0}}; 

  std::vector<double> matmass_linear[4] = {{0.088, 0.088, 0.04},
				   {0.2, 0.12, 0.04},
				   {0.2, 0.12, 0.04},
				   {0.312, 0.152, 0.04}}; 

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0};

  // Compute density
  double compute_density(Portage::Point<2> &coords, int matid, DENSITY_FUNCTION dtype)
  {
    double value = 0;
    switch(dtype)
    {
      case CONSTANT : value = matdensity_const[matid]; break;
      case LINEAR : value = coords[0]+coords[1]; break;
    }
    return value;
  }
}; //material_metadata

template<>
struct material_metadata<3, MAT_CONFIG::NESTED_BOX>
{
  const int nmats = 3;
  std::string matnames[3] = {"mat0", "mat1", "mat2"};

  // Extents of the materials in the overall domain
  std::vector<Portage::Point<3>> matlo[3], mathi[3];

  material_metadata(){
    double x[2][4] = {{0, 0.2, 0.8, 1.0}, {0.2, 0.4, 0.6, 0.8}};
    double y[2][4] = {{0, 0.2, 0.8, 1.0}, {0.2, 0.4, 0.6, 0.8}};
    double z[2][4] = {{0, 0.2, 0.8, 1.0}, {0.2, 0.4, 0.6, 0.8}};

    matlo[0] = {Portage::Point<3>(x[0][0], y[0][0], z[0][0]),
                Portage::Point<3>(x[0][0], y[0][2], z[0][0]),
                Portage::Point<3>(x[0][0], y[0][1], z[0][0]),
                Portage::Point<3>(x[0][2], y[0][1], z[0][0]),
                Portage::Point<3>(x[0][1], y[0][1], z[0][0]),
                Portage::Point<3>(x[0][1], y[0][1], z[0][2])};

    mathi[0] = {Portage::Point<3>(x[0][3], y[0][1], z[0][3]),
                Portage::Point<3>(x[0][3], y[0][3], z[0][3]),
                Portage::Point<3>(x[0][1], y[0][2], z[0][3]),
                Portage::Point<3>(x[0][3], y[0][2], z[0][3]),
                Portage::Point<3>(x[0][2], y[0][2], z[0][1]),
                Portage::Point<3>(x[0][2], y[0][2], z[0][3])};

    matlo[1] = {Portage::Point<3>(x[1][0], y[1][0], z[1][0]),
                Portage::Point<3>(x[1][0], y[1][2], z[1][0]),
                Portage::Point<3>(x[1][0], y[1][1], z[1][0]),
                Portage::Point<3>(x[1][2], y[1][1], z[1][0]),
                Portage::Point<3>(x[1][1], y[1][1], z[1][0]),
                Portage::Point<3>(x[1][1], y[1][1], z[1][2])};

    mathi[1] = {Portage::Point<3>(x[1][3], y[1][1], z[1][3]),
                Portage::Point<3>(x[1][3], y[1][3], z[1][3]),
                Portage::Point<3>(x[1][1], y[1][2], z[1][3]),
                Portage::Point<3>(x[1][3], y[1][2], z[1][3]),
                Portage::Point<3>(x[1][2], y[1][2], z[1][1]),
                Portage::Point<3>(x[1][2], y[1][2], z[1][3])};

    matlo[2] = {Portage::Point<3>(0.4, 0.4, 0.4)};
    mathi[2] = {Portage::Point<3>(0.6, 0.6, 0.6)};
  }

  int mat_nboxes[3] = {6,6,1};

  // Reference Volume
  std::vector<double> matvol[4] = {{0.264, 0.088, 0.008},
				   {0.264, 0.088, 0.008},
				   {0.264, 0.088, 0.008},
				   {0.264, 0.088, 0.008}}; 

  // Reference Mass
  std::vector<double> matmass_const[4] = {{0.0264, 0.88, 0.8},
				   {0.0264, 0.88, 0.8},
				   {0.0264, 0.88, 0.8},
				   {0.0264, 0.88, 0.8}}; 

  std::vector<double> matmass_linear[4] = {{0.2712, 0.1128, 0.012},
				   {0.396, 0.132, 0.012},
				   {0.396, 0.132, 0.012},
				   {0.5208, 0.1512, 0.012}}; 

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0};

  // Compute density
  double compute_density(Portage::Point<3> coords, int matid, DENSITY_FUNCTION dtype)
  {
    double value = 0;
    switch(dtype)
    {
      case CONSTANT : value = matdensity_const[matid]; break;
      case LINEAR : value = coords[0]+coords[1]+coords[2]; break;
    }
    return value;
  }
}; // material_metadata

template<int dim, MAT_CONFIG mconfig> 
void compute_material_fields_on_mesh(
  Wonton::Jali_Mesh_Wrapper &MeshWrapper, 
  Portage::Entity_type etype,
  material_metadata<dim, mconfig> &mdata,  
  std::vector<int> (&matcells)[3],
  std::vector<double> (&matvf)[3], 
  std::vector<Portage::Point<dim>> (&matcen)[3], 
  std::vector<double> (&matrho)[3], 
  DENSITY_FUNCTION dtype)
{
  int ncells = MeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                        etype);
  
  for (int c = 0; c < ncells; c++) {
    std::vector<Portage::Point<dim>> ccoords;
    MeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = MeshWrapper.cell_volume(c);

    Portage::Point<dim> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<dim>(ccoords, &cell_lo, &cell_hi);

    std::vector<double> xmoments;
    for (int m = 0; m < mdata.nmats; m++) {
     double vol = 0.0;    
     Portage::Point<dim> mcen; 
     for (int d = 0; d < dim; d++) mcen[d] = 0.0; 

     for (int nb = 0; nb < mdata.mat_nboxes[m]; nb++){

      if (BOX_INTERSECT::intersect_boxes<dim>(mdata.matlo[m][nb], mdata.mathi[m][nb],
                                            cell_lo, cell_hi, &xmoments)) {
        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          vol += xmoments[0];
          for (int d = 0; d < dim; d++)
            mcen[d] += xmoments[d+1];
        }
       }
      }//num_boxes

      if (vol > 1.0e-06) {  // non-trivial intersection
          matcells[m].push_back(c);
          matvf[m].push_back(vol/cellvol);
          for (int d = 0; d < dim; d++)
            mcen[d] = mcen[d]/vol;
          matcen[m].push_back(mcen);
          double val = mdata.compute_density(mcen, m, dtype);
          matrho[m].push_back(val);
        }
     } //material_loop
  }//loop over cells
} //compute_material_fields_on_mesh

template<int dim, MAT_CONFIG mconfig>
void check_fields_after_remap(
  Wonton::Jali_Mesh_Wrapper &targetMeshWrapper, 
  Wonton::Jali_State_Wrapper &targetStateWrapper,
  material_metadata<dim, mconfig> &mdata, 
  DENSITY_FUNCTION dtype, bool dbg_print=false)
{
    int commRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    //Compute reference values by using exact intersections of 
    //target mesh with the given material bounds
    constexpr int nmats = 3; 
    std::vector<int> matcells_ref[nmats];
    std::vector<double> matvf_ref[nmats];
    std::vector<Portage::Point<dim>> matcen_ref[nmats];
    std::vector<double> matrho_ref[nmats];

    compute_material_fields_on_mesh<dim, mconfig>(targetMeshWrapper, 
    Portage::Entity_type::PARALLEL_OWNED, mdata, 
    matcells_ref, matvf_ref, matcen_ref, matrho_ref, dtype); 

    if (dbg_print) {
      for (int m = 0; m < nmats; m++) {

        // Get target material cells and their vf, centroids, density values.
        std::vector<int> matcells_remap;
        double const *matvf_remap;
        double const *matrho_remap;
        Portage::Point<2> const *matcen_remap;

        targetStateWrapper.mat_get_cells(m, &matcells_remap);
        targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);
        targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);
        targetStateWrapper.mat_get_celldata("density", m, &matrho_remap);

        // Print details
        std::cerr<<"Mat ID "<<m<<std::endl;
        std::cerr<<"matcells_ref.size() = "<<matcells_ref[m].size()<<
                   ", matcells_remap.size() = "<< matcells_remap.size()<<std::endl;

        std::cerr<<" Reference cell id :: { vf,  centroid, density}"<<std::endl;

        int const num_current_matcells_ref = matcells_ref[m].size();
        int const num_matcells_remap = matcells_remap.size();

        for (int c = 0; c < num_current_matcells_ref; c++) {
          int gid = targetMeshWrapper.get_global_id(matcells_ref[m][c], Wonton::CELL);
          std::cerr << gid << ":: { " << matvf_ref[m][c] <<", {"
                    << matcen_ref[m][c][0] <<", " << matcen_ref[m][c][1] <<"}, "
                    << matrho_ref[m][c]<<" }" << std::endl;
        }
        std::cerr<<std::endl;

        std::cerr<<" Remapped cell id :: { vf,  centroid, density}"<<std::endl;
        for (int c = 0; c < num_matcells_remap; c++) {
          int gid = targetMeshWrapper.get_global_id(matcells_remap[c], Wonton::CELL);
          std::cerr << gid << ":: { " << matvf_remap[c] << ", {"
                    << matcen_remap[c][0] <<", " << matcen_remap[c][1] <<"}, "
                    << matrho_remap[c] << " }" << std::endl;
        }
    } 
  } else {

    //------------------------------------------------------------------------
    // CHECK 1: Number of materials on target 
    //-----------------------------------------------------------------------

    ASSERT_EQ(mdata.nmats, targetStateWrapper.num_materials());
    for (int m = 0; m < nmats; m++)
      ASSERT_EQ(mdata.matnames[m], targetStateWrapper.material_name(m));

    //-----------------------------------------------------------------------
    // CHECK 2: Cells in each material on target 
    //-----------------------------------------------------------------------

    std::vector<int> matcells_remap[nmats];
    for (int m = 0; m < nmats; m++) {
      targetStateWrapper.mat_get_cells(m, &matcells_remap[m]);
      int nmatcells = matcells_remap[m].size();

      for (int ic = 0; ic < nmatcells; ic++)
        ASSERT_EQ(matcells_remap[m][ic], matcells_ref[m][ic]);
    }

    //------------------------------------------------------------------------
    // CHECK 3: Material volume fractions and density
    // in each cell in target 
    //------------------------------------------------------------------------
    for (int m = 0; m < nmats; m++) {
      int nmatcells = matcells_remap[m].size();

      double const *matvf_remap;
      targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

      for (int ic = 0; ic < nmatcells; ic++)
        ASSERT_NEAR( matvf_remap[ic], matvf_ref[m][ic], 1.0e-09);

      // MOF cannot match moments and centroids as well as it can volume
      // fractions - so use looser tolerances
      // In this test, we don't check centroids for the nested box configuration 
      // as there will be cells where two materials interface is a box corner, 
      // and MOF will approximate it with a linear interface within that cell, 
      // and the centroid of the materials in the splitted cell would not match
      // the exact one. 
      
      if (mconfig == MAT_CONFIG::LAYER) {

        Portage::Point<dim> const *matcen_remap;
        targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);

        for (int ic = 0; ic < nmatcells; ic++)
          for (int d = 0; d < dim; d++)
            ASSERT_NEAR(matcen_remap[ic][d], matcen_ref[m][ic][d], 1.0e-9);
      }

      double const *matrho_remap;
      targetStateWrapper.mat_get_celldata("density", m, &matrho_remap);

      for (int ic = 0; ic < nmatcells; ic++)
        ASSERT_NEAR( matrho_remap[ic], matrho_ref[m][ic], 1.0e-10);
    }
 } //dbg_print     
} //check_fields_after_remap

template<int dim, MAT_CONFIG mconfig>
void run(std::shared_ptr<Jali::Mesh> &sourceMesh, 
    std::shared_ptr<Jali::Mesh> &targetMesh,
    std::shared_ptr<Jali::State> &sourceState,
    std::shared_ptr<Jali::State> &targetState,
    DENSITY_FUNCTION dtype)
{
  int nranks = 0; 
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  int commRank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  std::cout<<"On rank "<<commRank<<std::endl;

  //------------------------------------------------------------------------
  // Create source/target mesh/state wrappers
  //------------------------------------------------------------------------
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  //------------------------------------------------------------------------
  // Create material metadata and obtain field values for source
  //------------------------------------------------------------------------
  material_metadata<dim, mconfig> mdata; 
  constexpr int nmats = 3; 
  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Portage::Point<dim>> matcen_src[nmats];
  std::vector<double> matrho_src[nmats];

  compute_material_fields_on_mesh<dim, mconfig>(sourceMeshWrapper, Portage::Entity_type::ALL, 
  mdata, matcells_src, matvf_src, matcen_src, matrho_src, dtype);

  //------------------------------------------------------------------------
  // Add fields to source state 
  //------------------------------------------------------------------------
   for (int m = 0; m < nmats; m++)
    sourceStateWrapper.add_material(mdata.matnames[m], matcells_src[m]);

  // Create multi-material variables to store the volume fractions and
  // centroids for each material in the cells
  for (int m = 0; m < nmats; m++) {
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &(matvf_src[m][0]));
    sourceStateWrapper.mat_add_celldata("mat_centroids", m, &(matcen_src[m][0]));
  }

  // Add density profile for each material
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

    double volume = 0.0;
    double mass = 0.0;
    int const num_matcells = matcells.size();

    for (int ic = 0; ic < num_matcells; ic++) {
      double cellvol = vf[ic]*sourceMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
    }
  
#ifndef NDEBUG
    std::cerr<<std::setprecision(13)<<"On rank "<<commRank<<" Mat "<<m
    <<": volume = "<<volume<<", mass = "<<mass<<std::endl;  
 
#endif 
   
    ASSERT_NEAR(mdata.matvol[commRank][m], volume, 1.0e-12);

    if (dtype == CONSTANT)
      ASSERT_NEAR(mdata.matmass_const[commRank][m], mass, 1.0e-12);
    else 
      ASSERT_NEAR(mdata.matmass_linear[commRank][m], mass, 1.0e-12);
   
  }

  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------
  std::vector<std::string> remap_fields = {"density"};
  std::vector<int> dummymatcells;
  for (int m = 0; m < nmats; m++)
   targetStateWrapper.add_material(mdata.matnames[m], dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Portage::Point<dim>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  // Run the remap driver using MOF as the interface
  // reconstructor which is guaranteed to recover the correct
  // geometry of the linear interfaces. 
  //-------------------------------------------------------------------
  int numpe; 
  MPI_Comm_size(MPI_COMM_WORLD, &numpe); 

  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  Wonton::Executor_type *executor = (numpe > 1) ? &mpiexecutor : nullptr;

  if ((dim == 2) && (dtype == CONSTANT)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR2D,
		      Portage::Interpolate_1stOrder, 2,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.set_bnd_limiter(Portage::Boundary_Limiter_type::BND_NOLIMITER);
    d.run(executor);  // run in parallel
  } else if ((dim == 2) && (dtype == LINEAR)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR2D,
		      Portage::Interpolate_2ndOrder, 2,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.set_bnd_limiter(Portage::Boundary_Limiter_type::BND_NOLIMITER);
    d.run(executor);  // run in parallel
  } else if ((dim == 3) && (dtype == CONSTANT)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR3D,
		      Portage::Interpolate_1stOrder, 3,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitRnD<3>, Tangram::ClipRnD<3>>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.set_bnd_limiter(Portage::Boundary_Limiter_type::BND_NOLIMITER);
    d.run(executor);  // run in parallel
  } else if ((dim == 3) && (dtype == LINEAR)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR3D,
		      Portage::Interpolate_2ndOrder, 3,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitRnD<3>, Tangram::ClipRnD<3>>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.set_bnd_limiter(Portage::Boundary_Limiter_type::BND_NOLIMITER);
    d.run(executor);  // run in parallel
  } else
   std::cerr<<"Remapping requested for dim != 2 or 3"<<std::endl;

  //-------------------------------------------------------------------
  // Check remapping results on target mesh 
  //-------------------------------------------------------------------
   check_fields_after_remap<dim, mconfig>(targetMeshWrapper, targetStateWrapper, mdata, 
   dtype);
} //run

TEST(MMDriver2D, Layer_Const1stOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 7, 7);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap 
  run<2, LAYER>(sourceMesh, targetMesh, sourceState, targetState, CONSTANT);
}

TEST(MMDriver2D, Layer_Linear2ndOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 7, 7);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<2, LAYER>(sourceMesh, targetMesh, sourceState, targetState, LINEAR);
}

TEST(MMDriver3D, Layer_Const1stOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 7, 7);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<3, LAYER>(sourceMesh, targetMesh, sourceState, targetState, CONSTANT);
}


TEST(MMDriver3D, Layer_Linear2ndOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 7, 7);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<3, LAYER>(sourceMesh, targetMesh, sourceState, targetState, LINEAR);
}


TEST(MMDriver2D, NestedBox_Const1stOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 10, 10);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<2, NESTED_BOX>(sourceMesh, targetMesh, sourceState, targetState, CONSTANT);
}

TEST(MMDriver2D, NestedBox_Linear2ndOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 10, 10);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<2, NESTED_BOX>(sourceMesh, targetMesh, sourceState, targetState, LINEAR);
}

TEST(MMDriver3D, NestedBox_Const1stOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<3, NESTED_BOX>(sourceMesh, targetMesh, sourceState, targetState, CONSTANT);
}

TEST(MMDriver3D, NestedBox_Linear2ndOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<3, NESTED_BOX>(sourceMesh, targetMesh, sourceState, targetState, LINEAR);
}

// 2-material problem with VOF

TEST(MMDriver, TwoMat2D_VOF_1stOrderRemap) {

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
  //     |             :            |
  //     |             :            |
  //     |             :            |
  //     |             :            |
  //     |     0       +     1      |     
  //     |   mat0      :   mat1     |
  //     |             :            |
  //     |             :            |
  //     |             :            |
  //     |             :            |
  //     *-------------:------------*
  //    0,0           0.5,0         1,0

  constexpr int nmats = 2;
  std::string matnames[nmats] = {"mat0", "mat1"};

  // Extents of the materials in the overall domain

  Wonton::Point<2> matlo[nmats], mathi[nmats];
  matlo[0] = Wonton::Point<2>(0.0, 0.0);
  mathi[0] = Wonton::Point<2>(0.5, 1.0);
  matlo[1] = Wonton::Point<2>(0.5, 0.0);
  mathi[1] = Wonton::Point<2>(1.0, 1.0);


  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<double> matrho_src[nmats];


  //-------------------------------------------------------------------
  // COMPUTE MATERIAL DATA ON SOURCE SIDE - VOLUME FRACTIONS, CENTROID
  // CELL LISTS
  //-------------------------------------------------------------------
  //
  // Based on the material geometry, make a list of SOURCE cells that
  // are in each material and collect their material volume fractions

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
          matrho_src[m].push_back(m);
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
  // for each material in the cells
  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.mat_add_celldata("mat_volfracs", m, &(matvf_src[m][0]));

  // Also assign a different constant density value for each material
  for (int m = 0; m < nmats; m++)
    sourceStateWrapper.mat_add_celldata("density", m, &(matrho_src[m][0]));


  //-------------------------------------------------------------------
  // Add the materials and fields to the target mesh.
  //-------------------------------------------------------------------

  std::vector<int> dummymatcells;
  targetStateWrapper.add_material("mat0", dummymatcells);
  targetStateWrapper.add_material("mat1", dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);


  //-------------------------------------------------------------------
  // The driver is no longer responsible for distributing the mesh.
  // The calling program is responsible, so we create the flat mesh/state
  // prior to instantiating the driver distribute here.
  //-------------------------------------------------------------------

  // create the flat mesh
  Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
  source_mesh_flat.initialize(sourceMeshWrapper);

  // create the flat state
  Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);

  // mat_volfracs and mat_centroids are always imported from the state wrapper
  source_state_flat.initialize(sourceStateWrapper, {"density"});

  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  
  // Use a bounding box distributor to send the source cells to the target
  // partitions where they are needed
  Portage::MPI_Bounding_Boxes distributor(&mpiexecutor);
  distributor.distribute(source_mesh_flat, source_state_flat, targetMeshWrapper,
                         targetStateWrapper);


  //-------------------------------------------------------------------
  //  Run the remap driver using XMOF-2D as the interface
  //  reconstructor which is guaranteed to recover the correct
  //  geometry of the T-junction. A VOF-based interface reconstructor
  //  will get the right geometry for the T-junction only if the
  //  material ordering is 0, 1, 2
  //-------------------------------------------------------------------

  Portage::UberDriver<2,
                      Wonton::Flat_Mesh_Wrapper<>, 
                      Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>>,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Tangram::VOF, Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>
    d(source_mesh_flat, source_state_flat,
        targetMeshWrapper, targetStateWrapper, &mpiexecutor);


  d.compute_interpolation_weights<Portage::SearchKDTree, Portage::IntersectR2D>();

  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate<double, Portage::Entity_kind::CELL,
                Portage::Interpolate_2ndOrder>("density", "density",
  					       0.0, dblmax,
                                               Portage::Limiter_type::NOLIMITER,
                                               Portage::Boundary_Limiter_type::BND_NOLIMITER);

}  // TwoMat2D_VOF_1stOrderRemap

#endif  // ifdef PORTAGE_HAS_TANGRAM
