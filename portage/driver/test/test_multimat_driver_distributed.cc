#ifdef HAVE_TANGRAM

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
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

  // Reference Volume 
  std::map<int, int> ref_rank = {{1,0}, {2,1}, {4,2}}; 

  std::vector<double> matvol[3][3] = {{{ 0.64}, { 0.3485714285714, 0.4057142857143}, 
                    		       { 0.1885714285714, 0.2171428571429, 
					 0.2171428571429, 0.2457142857143 }}, 
    				      {{ 0.32}, { 0.1885714285714, 0.2685714285714}, 
				       { 0.1085714285714, 0.1567346938776,  
					 0.1567346938776, 0.2244897959184 }},
				      {{ 0.04}, { 0.03428571428571, 0.04}, 
				       { 0.02938775510204, 0.03428571428571, 
					 0.03428571428571, 0.04 }}}; 

  // Reference Mass
  std::vector<double> matmass_const[3][3] = 
				    {{{ 0.064}, { 0.03485714285714, 0.04057142857143}, 
                    		       { 0.01885714285714, 0.02171428571429, 
					 0.02171428571429, 0.02457142857143 }}, 
    				      {{ 3.2}, { 1.885714285714, 2.685714285714}, 
				       { 1.085714285714, 1.567346938776,  
					 1.567346938776, 2.244897959184 }},
				      {{ 4.0}, { 3.428571428571, 4.0}, 
				       { 2.938775510204, 3.428571428571, 
					 3.428571428571, 4.0 }}}; 

  //double matmass_linear[3][3] = {{0.64, 0.32, 0.04}, {}, {}};

  std::vector<double> matmass_linear[3][3] = 
				    {{{ 0.64}, { 0.2515918367347, 0.4945306122449 }, 
                    		       { 0.08016326530612, 0.2016326530612 , 
					 0.2016326530612, 0.3688163265306 }}, 
    				      {{ 0.32}, {  0.1635918367347, 0.2817959183673 }, 
				       { 0.07787755102041, 0.1435801749271 ,  
					 0.1435801749271, 0.2471603498542 }},
				      {{ 0.04}, { 0.03379591836735, 0.04 }, 
				       { 0.02854810495627, 0.03379591836735 , 
					 0.03379591836735, 0.04 }}}; 

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0}; 

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

  // Reference Volume 
  std::map<int, int> ref_rank = {{1,0}, {2,1}, {4,2}}; 

  std::vector<double> matvol[3][3] = {{{ 0.784}, { 0.4377142857143, 0.5291428571429}, 
                    		       { 0.2437551020408, 0.2935510204082 , 
					 0.2935510204082, 0.3515102040816 }}, 
    				      {{ 0.208}, { 0.1268571428571, 0.1771428571429}, 
				       { 0.07689795918367, 0.1077551020408,  
					 0.1077551020408, 0.150693877551}},
				      {{ 0.008}, { 0.006857142857143, 0.008}, 
				       { 0.005877551020408, 0.006857142857143, 
					 0.006857142857143,  0.008}}}; 

  // Reference Mass
  std::vector<double> matmass_const[3][3] = 
				    {{{ 0.0784}, { 0.04377142857143, 0.05291428571429}, 
                    		       { 0.02437551020408, 0.02935510204082, 
					 0.02935510204082, 0.03515102040816}}, 
    				      {{ 2.08}, { 1.268571428571, 1.771428571429}, 
				       { 0.7689795918367, 1.077551020408,  
					 1.077551020408, 1.50693877551}},
				      {{ 0.8}, { 0.6857142857143, 0.8}, 
				       { 0.5877551020408, 0.6857142857143, 
					 0.6857142857143, 0.8}}}; 

  std::vector<double> matmass_linear[3][3] = 
				    {{{ 1.176}, { 0.5494040816327, 0.8878204081633}, 
                    		       { 0.2446110787172, 0.4193586005831, 
					 0.4193586005831, 0.6594355685131}}, 
    				      {{ 0.312}, {  0.1751020408163, 0.2736489795918}, 
				       { 0.09659475218659, 0.1535440233236,  
					 0.1535440233236, 0.239643148688}},
				      {{ 0.012}, { 0.01018775510204, 0.012 }, 
				       { 0.008648396501458, 0.01018775510204, 
					 0.01018775510204, 0.012}}}; 

  // Constant density
  double matdensity_const[3] = {0.1, 10.0, 100.0};
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

    std::vector<double> xmoments;
    for (int m = 0; m < mdata.nmats; m++) {
     int cnt = 0; 
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
      ASSERT_NEAR( matvf_remap[ic], matvf_ref[m][ic], 1.0e-10);

    // MOF cannot match moments and centroids as well as it can volume
    // fractions - so use looser tolerances
    // In this test, we don't check centroids at all as there will be cells
    // where two materials will form a t-junction, and MOF will approximate
    // it with a linear interface within that cell, and the centroid of the 
    // materials in the divided cell would not match the exact one. 
    //
    //Portage::Point<dim> const *matcen_remap;
    //targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);
    //for (int ic = 0; ic < nmatcells; ic++)
    //  for (int d = 0; d < dim; d++)
    //    ASSERT_NEAR(matcen_remap[ic][d], matcen_ref[m][ic][d], 1.0e-9);

    double const *matrho_remap;
    targetStateWrapper.mat_get_celldata("density", m, &matrho_remap);

    for (int ic = 0; ic < nmatcells; ic++)
      ASSERT_NEAR( matrho_remap[ic], matrho_ref[m][ic], 1.0e-10);
  }
  
#ifdef DEBUG
   /*
   for (int m = 0; m < nmats; m++) {
    int nmatcells = matcells_remap[m].size();

    double const *matvf_remap;
    targetStateWrapper.mat_get_celldata("mat_volfracs", m, &matvf_remap);

    Portage::Point<2> const *matcen_remap;
    targetStateWrapper.mat_get_celldata("mat_centroids", m, &matcen_remap);

    double const *matrho_remap;
    targetStateWrapper.mat_get_celldata("density", m, &matrho_remap);

    std::cerr<<"Mat ID "<<m<<std::endl;
    for (int ic = 0; ic < nmatcells; ic++)
    {
      std::cerr<<"----> Ref for cell id  "<<matcells_ref[m][ic]<<" = {";
      std::cerr<<matvf_ref[m][ic]<<", "<<
                 matcen_ref[m][ic][0]<<", "<<matcen_ref[m][ic][1]<<", "<<
                 matrho_ref[m][ic]<<" }"<<std::endl; 
  
      std::cerr<<"----> After remap for cell id  "<<matcells_remap[m][ic]<<" = {";
      std::cerr<<matvf_remap[ic]<<", "<<
                 matcen_remap[ic][0]<<", "<<matcen_remap[ic][1]<<", "<<
                 matrho_remap[ic]<<" }"<<std::endl; 
    } 
  }*/
#endif 
} //check_fields_after_remap

template<int dim>
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
  material_metadata<dim> mdata; 
  constexpr int nmats = 3; 
  std::vector<int> matcells_src[nmats];
  std::vector<double> matvf_src[nmats];
  std::vector<Portage::Point<dim>> matcen_src[nmats];
  std::vector<double> matrho_src[nmats];

  compute_material_fields_on_mesh(sourceMeshWrapper, mdata, matcells_src,
  matvf_src, matcen_src, matrho_src, dtype);

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
  
#ifdef DEBUG
 /*  
    std::cerr<<std::setprecision(13)<<"On rank "<<commRank<<" Mat "<<m
    <<": volume = "<<volume<<", mass = "<<mass<<std::endl;  
 */
#endif 

    int id = mdata.ref_rank.find(nranks)->second; 
    ASSERT_NEAR(mdata.matvol[m][id][commRank], volume, 1.0e-12);

    if (dtype == CONSTANT)
      ASSERT_NEAR(mdata.matmass_const[m][id][commRank], mass, 1.0e-12);
    else 
      ASSERT_NEAR(mdata.matmass_linear[m][id][commRank], mass, 1.0e-12);
  }

#ifdef DEBUG
  /*
  std::cerr<<"****** BEFORE REMAP ******"<<std::endl;
  int nsrccells = sourceMeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                                 Portage::Entity_type::ALL);
  std::cerr<<"SOURCE MESH ---->"<<std::endl;
  for (int c = 0; c < nsrccells; c++) {
    int gid = sourceMeshWrapper.get_global_id(c, Portage::Entity_kind::CELL);
    std::cerr<<" ---- Cell "<<gid<<std::endl;
  } 

  int  ntgtcells = targetMeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                                  Portage::Entity_type::ALL);
  std::cerr<<"TARGET MESH ---->"<<std::endl;
  for (int c = 0; c < ntgtcells; c++) {
    int gid = targetMeshWrapper.get_global_id(c, Portage::Entity_kind::CELL);
    std::cerr<<" ---- Cell "<<gid<<std::endl;
  }
  */
#endif

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
		      Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.run(executor);  // run in parallel
  }
  else if ((dim == 2) && (dtype == LINEAR)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR2D,
		      Portage::Interpolate_2ndOrder, 2,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.run(executor);  // run in parallel
  }
  else if ((dim == 3) && (dtype == CONSTANT)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR3D,
		      Portage::Interpolate_1stOrder, 3,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.run(executor);  // run in parallel
  }
  else if ((dim == 3) && (dtype == LINEAR)){
    Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR3D,
		      Portage::Interpolate_2ndOrder, 3,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
		      Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    d.set_limiter(Portage::Limiter_type::NOLIMITER);
    d.run(executor);  // run in parallel
  }
 else
   std::cerr<<"Remapping requested for dim != 2 or 3"<<std::endl;

#ifdef DEBUG
  /*
  std::cerr<<"****** AFTER REMAP ******"<<std::endl;
  nsrccells = sourceMeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                                 Portage::Entity_type::ALL);
  std::cerr<<"SOURCE MESH ---->"<<std::endl;
  for (int c = 0; c < nsrccells; c++) {
    int gid = sourceMeshWrapper.get_global_id(c, Portage::Entity_kind::CELL);
    std::cerr<<" ---- Cell "<<gid<<std::endl;
  } 

  ntgtcells = targetMeshWrapper.num_entities(Portage::Entity_kind::CELL,
                                                 Portage::Entity_type::ALL);
  std::cerr<<"TARGET MESH ---->"<<std::endl;
  for (int c = 0; c < ntgtcells; c++) {
    int gid = targetMeshWrapper.get_global_id(c, Portage::Entity_kind::CELL);
    std::cerr<<" ---- Cell "<<gid<<std::endl;
  } 
  */
#endif

  //-------------------------------------------------------------------
  // Check remapping results on target mesh 
  //-------------------------------------------------------------------
   check_fields_after_remap(targetMeshWrapper, targetStateWrapper, mdata, 
   dtype);
} //run

TEST(MMDriver2D, Const1stOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 7, 7);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 7, 7);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap 
  run<2>(sourceMesh, targetMesh, sourceState, targetState, CONSTANT);
}

TEST(MMDriver2D, Linear2ndOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 7, 7);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 7, 7);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<2>(sourceMesh, targetMesh, sourceState, targetState, LINEAR);
}

/*
TEST(MMDriver3D, Const1stOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 7, 7);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 7, 7);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<3>(sourceMesh, targetMesh, sourceState, targetState, CONSTANT);
}


TEST(MMDriver3D, Linear2ndOrder)
{
  // Create source/target meshes and states
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;

  std::shared_ptr<Jali::Mesh> targetMesh;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK); 

  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 7, 7);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 7, 7);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  //Remap
  run<3>(sourceMesh, targetMesh, sourceState, targetState, LINEAR);
}
*/


#endif  // ifdef have_tangram
