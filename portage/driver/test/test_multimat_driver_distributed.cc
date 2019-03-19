#ifdef HAVE_TANGRAM

#include <cmath>
#include <iostream>
#include <memory>
#include <vector> 

//This test is distributed
#include "mpi.h"
#include "gtest/gtest.h"

//Tangram includes
#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"

//Portage includes
#include "portage/driver/mmdriver.h"
#include "portage/intersect/simple_intersect_for_tests.h"

//Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

//Wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// Tests for distributed multi-material remap with 1st and 
// 2nd Order Accurate Remap on 4 ranks. 

// The conceptual material layout is two nested boxes inside an
// outer box domain with unit length. The material interfaces
// align with the coordinate axes/planes. With this layout, all
// four ranks in the distributed case will have all three materials. 
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

// The three materials have densities initialized with either constant
// or linear functions. Given this setup we know the volume, mass and
// centroid of each material.

// Knowing the simple material geometry, we can do box intersections
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

template<int dim>
struct material_metadata {};

template <>
struct material_metadata<2>
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

  //Volume 
  double matvol[3] = {0.64, 0.32, 0.04};

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0};

  // Mass
  double matmass_const[3] = {0.064, 3.2, 4.0};
  double matmass_linear[3] = {0.64, 0.32, 0.04};

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
}; 

template<>
struct material_metadata<3>
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

  //Volume 
  double matvol[3] = {0.784, 0.208, 0.008};

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0};

  // Mass
  double matmass_const[3] = {0.0784, 2.08, 0.8};
  double matmass_linear[3] = {0.676, 0.312, 0.012};

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
};

template<int dim> 
void compute_material_fields_on_mesh(
  Wonton::Jali_Mesh_Wrapper &MeshWrapper, 
  material_metadata<dim> &mdata,  
  std::vector<int> (&matcells)[3],
  std::vector<double> (&matvf)[3], 
  std::vector<Point<dim>> (&matcen)[3], 
  std::vector<double> (&matrho)[3], 
  DENSITY_FUNCTION dtype)
{
  int ncells = MeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                        Portage::Entity_type::ALL);
  
  for (int c = 0; c < ncells; c++) {
    std::vector<Portage::Point<dim>> ccoords;
    MeshWrapper.cell_get_coordinates(c, &ccoords);

    double cellvol = MeshWrapper.cell_volume(c);

    Portage::Point<dim> cell_lo, cell_hi;
    BOX_INTERSECT::bounding_box<dim>(ccoords, &cell_lo, &cell_hi);

    //DBG
    std::cout<<"Cell "<<c<<" :: Vol = "<<cellvol<<std::endl;
    //DBG

    std::vector<double> xmoments;
    for (int m = 0; m < mdata.nmats; m++) {
     int cnt = 0; 
     double vol = 0.0;    
     Portage::Point<dim> mcen; 
     for (int d = 0; d < dim; d++) mcen[d] = 0.0; 

     for (int nb = 0; nb < mdata.mat_nboxes[m]; nb++){

      if (BOX_INTERSECT::intersect_boxes<dim>(mdata.matlo[m][nb], mdata.mathi[m][nb],
                                            cell_lo, cell_hi, &xmoments)) {

      //DBG
      std::cout<<"----> Mat "<<m<<" :: nb "<<nb<<" :: Vol = "<<xmoments[0]<<std::endl;
      //DBG


        if (xmoments[0] > 1.0e-06) {  // non-trivial intersection
          vol += xmoments[0];
          for (int d = 0; d < dim; d++)
            mcen[d] += xmoments[d+1];
          //cnt += 1; 
          //for (int d = 0; d < dim; d++)
          //  mcen[d] += xmoments[d+1]/xmoments[0];
          //cnt += 1; 
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
}

template<int dim>
void check_fields_after_remap(
  Wonton::Jali_Mesh_Wrapper &targetMeshWrapper, 
  Wonton::Jali_State_Wrapper &targetStateWrapper,
  material_metadata<dim> &mdata, 
  DENSITY_FUNCTION dtype)
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

  compute_material_fields_on_mesh(targetMeshWrapper, mdata, 
    matcells_ref, matvf_ref, matcen_ref, matrho_ref, dtype); 

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

    ASSERT_EQ(nmatcells, matcells_ref[m].size());

    //std::sort(matcells_remap[m].begin(), matcells_remap[m].end());
    //std::sort(matcells_ref[m].begin(), matcells_ref[m].end());

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_EQ(matcells_remap[m][ic], matcells_ref[m][ic]);
  }

  //------------------------------------------------------------------------
  // CHECK 3: Material volume fractions, centroids and density
  // in each cell in target 
  //------------------------------------------------------------------------
  /*for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    //for (int ic = 0; ic < nmatcells; ic++)
    //  ASSERT_NEAR( matvf_remap[ic], matvf_ref[m][ic], 1.0e-12);

    Portage::Point<2> const *matcen_remap;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);

    // MOF cannot match moments and centroids as well as it can volume
    // fractions - so use looser tolerances
    for (int ic = 0; ic < nmatcells; ic++)
      for (int d = 0; d < dim; d++)
        ASSERT_NEAR(matcen_remap[ic][d], matcen_ref[m][ic][d], 1.0e-9);

    double const *matrho_remap;
    targetStateWrapper.mat_get_celldata("density", m, &matrho_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR( matrho_remap[ic], matrho_ref[m][ic], 1.0e-10);
  }
  */
  //DBG
  for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    Portage::Point<2> const *matcen_remap;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);

    double const *matrho_remap;
    targetStateWrapper.mat_get_celldata("density", m, &matrho_remap);

    std::cout<<"Mat ID "<<m<<std::endl;
    for (int ic = 0; ic < nmatcells; ic++)
    {
      std::cout<<"----> Ref for cell id  "<<matcells_ref[m][ic]<<" = {";
      std::cout<<matvf_ref[m][ic]<<", "<<
                 matcen_ref[m][ic][0]<<", "<<matcen_ref[m][ic][1]<<", "<<
                 matrho_ref[m][ic]<<" }"<<std::endl; 
  
      std::cout<<"----> After remap for cell id  "<<matcells_remap[m][ic]<<" = {";
      std::cout<<matvf_remap[ic]<<", "<<
                 matcen_remap[ic][0]<<", "<<matcen_remap[ic][1]<<", "<<
                 matrho_remap[ic]<<" }"<<std::endl; 
    } 


  }
  //DBG
  

  //------------------------------------------------------------------------
  // CHECK 4: Total mass, this step requires global communication
  //------------------------------------------------------------------------

}

TEST(MMDriver2D, Const1stOrder)
{
  //------------------------------------------------------------------------
  // Create source/target meshes and states
  //------------------------------------------------------------------------
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 7, 7);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 7, 7);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  //------------------------------------------------------------------------
  // Create material metadata and obtain field values for source
  //------------------------------------------------------------------------
  material_metadata<2> mdata; 
  constexpr int nmats = 3; 
  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Portage::Point<2>> matcen_src[nmats];
  std::vector<double> matrho_src[nmats];

  compute_material_fields_on_mesh(sourceMeshWrapper, mdata, matcells_src,
    matvf_src, matcen_src, matrho_src, CONSTANT);

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

    double volume = 0.0, mass = 0.0;
    for (int ic = 0; ic < matcells.size(); ic++) {
      double cellvol = vf[ic]*sourceMeshWrapper.cell_volume(matcells[ic]);
      volume += cellvol;
      mass += rho[ic]*cellvol;
    }

    ASSERT_NEAR(mdata.matvol[m], volume, 1.0e-12);
    ASSERT_NEAR(mdata.matmass_const[m], mass, 1.0e-12);
  }





  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------
  std::vector<std::string> remap_fields = {"density"};
  std::vector<int> dummymatcells;
  for (int m = 0; m < nmats; m++)
   targetStateWrapper.add_material(mdata.matnames[m], dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Portage::Point<2>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  // Run the remap driver using XMOF-2D as the interface
  // reconstructor which is guaranteed to recover the correct
  // geometry of the linear interfaces. 
  //-------------------------------------------------------------------
  int numpe; 
  MPI_Comm_size(MPI_COMM_WORLD, &numpe); 

  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  Wonton::Executor_type *executor = (numpe > 1) ? &mpiexecutor : nullptr;

  Portage::MMDriver<Portage::SearchKDTree,
                    Portage::IntersectR2D,
                    Portage::Interpolate_1stOrder,
                    2,
                    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                    Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);
  d.set_remap_var_names(remap_fields);
  d.set_limiter(Portage::Limiter_type::NOLIMITER);
  d.run(executor);  // run in parallel

  //-------------------------------------------------------------------
  // Check remapping results on target mesh 
  //-------------------------------------------------------------------
  check_fields_after_remap(targetMeshWrapper, targetStateWrapper, mdata, 
  DENSITY_FUNCTION::CONSTANT);
}
/*
TEST(MMDriver2D, Linear2ndOrder)
{
  //------------------------------------------------------------------------
  // Create source/target meshes and states
  //------------------------------------------------------------------------
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 10, 10);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 10, 10);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  //------------------------------------------------------------------------
  // Create material metadata and obtain field values for source
  //------------------------------------------------------------------------
  material_metadata<2> mdata; 
  constexpr int nmats = 3; 
  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Portage::Point<2>> matcen_src[nmats];
  std::vector<double> matrho_src[nmats];

  compute_material_fields_on_mesh(sourceMeshWrapper, mdata, matcells_src,
    matvf_src, matcen_src, matrho_src, LINEAR);

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
  // Field(s) we have to remap
  //-------------------------------------------------------------------
  std::vector<std::string> remap_fields = {"density"};
  std::vector<int> dummymatcells;
  for (int m = 0; m < nmats; m++)
   targetStateWrapper.add_material(mdata.matnames[m], dummymatcells);

  targetStateWrapper.mat_add_celldata<double>("mat_volfracs");
  targetStateWrapper.mat_add_celldata<Portage::Point<2>>("mat_centroids");
  targetStateWrapper.mat_add_celldata<double>("density", 0.0);

  //-------------------------------------------------------------------
  // Run the remap driver using XMOF-2D as the interface
  // reconstructor which is guaranteed to recover the correct
  // geometry of the linear interfaces. 
  //-------------------------------------------------------------------
  int numpe; 
  MPI_Comm_size(MPI_COMM_WORLD, &numpe); 

  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  Wonton::Executor_type *executor = (numpe > 1) ? &mpiexecutor : nullptr;

  Portage::MMDriver<Portage::SearchKDTree,
                    Portage::IntersectR2D,
                    Portage::Interpolate_2ndOrder,
                    2,
                    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                    Tangram::XMOF2D_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);
  d.set_remap_var_names(remap_fields);
  d.set_limiter(Portage::Limiter_type::NOLIMITER);
  d.run(executor);  // run in parallel

  //-------------------------------------------------------------------
  // Check remapping results on target mesh 
  //-------------------------------------------------------------------
  check_fields_after_remap(targetMeshWrapper, targetStateWrapper, mdata, 
  DENSITY_FUNCTION::LINEAR);
}

TEST(MMDriver3D, Const1stOrder)
{

}

TEST(MMDriver3D, Linear2ndOrder)
{

}
*/





#endif  // ifdef have_tangram
