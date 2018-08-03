/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#include "portage/wonton/state/jali/jali_state_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"


// Vector type for 2d doubles
struct Vec2d
{
  double x;
  double y;

  void set(double xvalue, double yvalue)
  {
    x = xvalue;  y = yvalue;
  }

  friend std::ostream &operator<<(std::ostream &output, const Vec2d &v)
  {
    output << "(" << v.x << ", " << v.y << ")";
    return output;
  }
};


TEST(Jali_State_Wrapper, DataTypes) {

  // Add multiple state vector types
  int const n_cells = 4;
  int const n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  double dtest[] = {62., 78., 43., 22.};
  int itest[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec2d vtest[n_cells];
  for (unsigned int i = 0; i < n_cells; i++) vtest[i].set(1.0*i, 2.0*i);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  Wonton::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  std::shared_ptr<Jali::State> state(Jali::State::create(inputMesh));
  Wonton::Jali_State_Wrapper wrapper(*state);

  state->add("f1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, ftest);
  state->add("i1", inputMesh, Jali::Entity_kind::NODE,
            Jali::Entity_type::ALL, itest);
  state->add("v1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, vtest);

  // Get raw float data using wrapper
  float* fdata;
  wrapper.mesh_get_data(Portage::CELL, "f1", &fdata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(fdata[i], ftest[i]);

  // Get raw int data using wrapper
  int* idata;
  wrapper.mesh_get_data(Portage::NODE, "i1", &idata);
  for (unsigned int i = 0; i < n_nodes; i++) ASSERT_EQ(idata[i], itest[i]);

  // Get raw Vec2d data using wrapper
  Vec2d* vdata;
  wrapper.mesh_get_data(Portage::CELL, "v1", &vdata);
  for (unsigned int i = 0; i < n_cells; i++) 
  {
    ASSERT_EQ(vdata[i].x, vtest[i].x);
    ASSERT_EQ(vdata[i].y, vtest[i].y);
  }

  // add data through wrapper
  wrapper.mesh_add_data<double>(Portage::CELL, "d1", dtest);
  double dval = 123.456;
  wrapper.mesh_add_data<double>(Portage::CELL, "d2", dval);

  double *ddata;
  wrapper.mesh_get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest[i]);
  wrapper.mesh_get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dval);

  // Iterate through a vector of names
  std::vector<std::string> fields;
  fields.push_back("v1");  fields.push_back("f1");  fields.push_back("i1");
  for (auto it = fields.begin(); it != fields.end(); it++)
  {
    Portage::Entity_kind on_what = wrapper.get_entity(*it);
    int nent = inputMeshWrapper.num_entities(on_what);
    if (typeid(float) == wrapper.get_data_type(*it))
    {
      float* fdata;
      wrapper.mesh_get_data(on_what, *it, &fdata);
      for (unsigned int i = 0; i < nent; i++) ASSERT_EQ(fdata[i], ftest[i]);
    }
    else if (typeid(int) == wrapper.get_data_type(*it))
    {
      int* idata;
      wrapper.mesh_get_data(on_what, *it, &idata);
      for (unsigned int i = 0; i < nent; i++) ASSERT_EQ(idata[i], itest[i]);
    }
    else if (typeid(Vec2d) == wrapper.get_data_type(*it))
    {
      Vec2d* vdata;
      wrapper.mesh_get_data(on_what, *it, &vdata);
      for (unsigned int i = 0; i < nent; i++)
      {
        ASSERT_EQ(vdata[i].x, vtest[i].x);
        ASSERT_EQ(vdata[i].y, vtest[i].y);
      }
    }
    else
    {
      ASSERT_EQ(0, 1);    // This else should never be reached in this test
    }
  }

  // Iterate through all fields using the wrapper
  for (auto it = wrapper.names_begin(); it != wrapper.names_end(); it++)
  {
    Portage::Entity_kind on_what = wrapper.get_entity(*it);
    int nent = inputMeshWrapper.num_entities(on_what);
    if (typeid(float) == wrapper.get_data_type(*it))
    {
      float* fdata;
      wrapper.mesh_get_data(on_what, *it, &fdata);
      for (unsigned int i = 0; i < nent; i++) ASSERT_EQ(fdata[i], ftest[i]);
    }
    else if (typeid(double) == wrapper.get_data_type(*it))
    {
      double* ddata;
      wrapper.mesh_get_data(on_what, *it, &ddata);
      if (*it == "d1")
        for (unsigned int i = 0; i < nent; i++) ASSERT_EQ(ddata[i], dtest[i]);
      else
        for (unsigned int i = 0; i < nent; i++) ASSERT_EQ(ddata[i], dval);
    }
    else if (typeid(int) == wrapper.get_data_type(*it))
    {
      int* idata;
      wrapper.mesh_get_data(on_what, *it, &idata);
      for (unsigned int i = 0; i < nent; i++) ASSERT_EQ(idata[i], itest[i]);
    }
    else if (typeid(Vec2d) == wrapper.get_data_type(*it))
    {
      Vec2d* vdata;
      wrapper.mesh_get_data(on_what, *it, &vdata);
      for (unsigned int i = 0; i < nent; i++)
      {
        ASSERT_EQ(vdata[i].x, vtest[i].x);
        ASSERT_EQ(vdata[i].y, vtest[i].y);
      }
    }
    else
    {
      ASSERT_EQ(0, 1);    // This else should never be reached in this test
    }
  }
}


TEST(Jali_State_Wrapper, MMState) {
  
  // Create a 3x3 mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 3.0, 3.0, 3, 3);

  ASSERT_NE(nullptr, mesh);

  // Mesh wrapper

  Wonton::Jali_Mesh_Wrapper meshwrapper(*mesh);

  // Create a state object

  std::shared_ptr<Jali::State> state = Jali::State::create(mesh);

  // State wrapper

  Wonton::Jali_State_Wrapper statewrapper(*state);

  int ncells = meshwrapper.num_owned_cells();

  // Define a single valued state vector on cells with string ids and
  // initialized to 0
  
  statewrapper.mesh_add_data(Portage::Entity_kind::CELL, "cell_density", 0.0);

  // Try to get const and non-const version of the raw data

  Wonton::Jali_State_Wrapper const statewrapper2(*state);

  //  double *rho;
  //  statewrapper2.mesh_get_data(Portage::Entity_kind::CELL, "cell_density", &rho);
  //  double dummyrho = rho[0];

  double const *constrho;
  statewrapper2.mesh_get_data(Portage::Entity_kind::CELL, "cell_density", &constrho);
  double const constdummyrho = constrho[0];

  

  // Define a multi-material state vector on cells to store volume fractions
  // Create 3 material sets in the state corresponding to a T-junction
  // configuration
  //
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    0   | 0  : 2 |    2   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*
  //     |        |    : 2 |    2   |
  //     |        |    +............|
  //     |    0   |  0 : 1 |    1   |
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    0   |  0 : 1 |    1   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*

  std::vector<std::vector<int>> matcells = {{0, 1, 3, 4, 6, 7},
                                            {1, 2, 4, 5},
                                            {4, 5, 7, 8}};
  std::vector<std::vector<int>> cellmats = {{0}, {0, 1}, {1},
                                            {0}, {0, 1, 2}, {1, 2},
                                            {0}, {0, 2}, {2}};

  statewrapper.add_material("steel1", matcells[0]);
  statewrapper.add_material("aluminum1", matcells[1]);

  int nmats = statewrapper.num_materials();
  ASSERT_EQ(2, nmats);

  statewrapper.add_material("aluminum1", matcells[2]);  // Name reused
  ASSERT_EQ(nmats, statewrapper.num_materials());  // no change expected

  statewrapper.add_material("aluminum2", matcells[2]);  // Unique ID and name
  nmats++;
  ASSERT_EQ(nmats, statewrapper.num_materials());  // should have incremented

  for (int i = 0; i < 9; i++) {
    ASSERT_EQ(cellmats[i].size(), statewrapper.cell_get_num_mats(i));
    std::vector<int> cellmats2;
    statewrapper.cell_get_mats(i, &cellmats2);
    for (int j = 0; j < cellmats[i].size(); j++)
      ASSERT_EQ(cellmats[i][j], cellmats2[j]);
  }


  ASSERT_EQ(statewrapper.material_name(0), "steel1");
  ASSERT_EQ(statewrapper.material_name(1), "aluminum1");
  ASSERT_EQ(statewrapper.material_name(2), "aluminum2");

  // Create a multi-material state vector corresponding to volume fractions
  // of materials as shown in fig above. Because we will add another material
  // later in the test, we will make space for 4 materials

  double **vf_in = new double*[4];
  for (int i = 0; i < 4; i++)
    vf_in[i] = new double[9];

  double vfarr[4][9] = {{1.0, 0.5, 0.0, 1.0, 0.5,  0.0, 1.0, 0.5, 0.0},
                        {0.0, 0.5, 1.0, 0.0, 0.25, 0.5, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.0, 0.5, 1.0}};
  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      vf_in[m][c] = vfarr[m][c];
  for (int c = 0; c < 9; c++)
    vf_in[3][c] = 0.0;

  statewrapper.mat_add_celldata("volfrac",
                                (double const **) vf_in);
  

  // Create a multi-material state vector corresponding to density of the
  // materials (rho_0 = 10.0; rho_1 = 2.0; rho_2 = 2.0). Similar to the
  // volume fraction array, we will keep a spot for the fourth material
  double rho_in[4] = {10.0, 2.0, 2.0, 0.0};

  statewrapper.mat_add_celldata<double>("mat_density");

  // Then update the data one material at a time
  for (int m = 0; m < nmats; m++) {
    statewrapper.mat_add_celldata("mat_density", m, rho_in[m]);
  }

  
  // Expected cell density calculated as the volume-weighted average
  // of material densities in the cell

  std::vector<double> rhocell_exp(ncells, 0.0);
  for (int c = 0; c < ncells; c++)
    for (int m = 0; m < nmats; m++)
      rhocell_exp[c] += rho_in[m]*vf_in[m][c];

  // Verify against values from state vectors

  std::vector<double> rhocell(ncells, 0.0);
  for (int m = 0; m < nmats; m++) {
    double const *rhomatptr;
    statewrapper.mat_get_celldata("mat_density", m, &rhomatptr);
    double const *vfmatptr;
    statewrapper.mat_get_celldata("volfrac", m, &vfmatptr);
    std::vector<int> matcells;
    statewrapper.mat_get_cells(m, &matcells);

    int nmatcells = matcells.size();
    for (int i = 0; i < nmatcells; i++) {
      int c = matcells[i];
      rhocell[c] += rhomatptr[i]*vfmatptr[i];
    }
  }
      
  for (int c = 0; c < ncells; c++)
    ASSERT_NEAR(rhocell_exp[c], rhocell[c], 1.0e-12);

  // Now lets introduce a 4th material into the state and have it
  // occupy the top left corner of the mesh, pushing material 1 down
  // to the bottom left corner
  //
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    3   | 3  : 2 |    2   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*
  //     |    3   | 3  : 2 |    2   |
  //     |.............+............|
  //     |    0   |  0 : 1 |    1   |
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    0   |  0 : 1 |    1   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*

  std::vector<int> matcells3 = {3 , 4, 6, 7};
  statewrapper.add_material("steel2", matcells3);
  nmats = statewrapper.num_materials();

  // Material 3 is added to cells 3, 4, 6, 7 but material 0 is not removed
  // (only its volume fraction is zeroed out)
  cellmats[3].push_back(3);
  cellmats[4].push_back(3);
  cellmats[6].push_back(3);
  cellmats[7].push_back(3);

  for (int i = 0; i < 9; i++) {
    ASSERT_EQ(cellmats[i].size(), statewrapper.cell_get_num_mats(i));
    std::vector<int> cellmats2;
    statewrapper.cell_get_mats(i, &cellmats2);
    for (int j = 0; j < cellmats[i].size(); j++)
      ASSERT_EQ(cellmats[i][j], cellmats2[j]);
  }

  // Set the volume fractions material 3 in cells in the input vectors
  vf_in[3][3] = 0.5;
  vf_in[3][4] = 0.25;
  vf_in[3][6] = 1.0;
  vf_in[3][7] = 0.5;

  // Set the volume fractions material 3 in cells in the state vectors
  double *vfmatptr;
  statewrapper.mat_get_celldata("volfrac", 3, &vfmatptr);
  vfmatptr[statewrapper.cell_index_in_material(3, 3)] = 0.5;
  vfmatptr[statewrapper.cell_index_in_material(4, 3)] = 0.25;
  vfmatptr[statewrapper.cell_index_in_material(6, 3)] = 1.0;
  vfmatptr[statewrapper.cell_index_in_material(7, 3)] = 0.5;

  // Set the density of material 3 to be the same as material 0
  rho_in[3] = rho_in[0];

  statewrapper.mat_add_celldata("mat_density", 3, rho_in[3]);

  // Adjust the volume fractions of material 0 in the input vector
  vf_in[0][3] = 0.5;
  vf_in[0][4] = 0.25;
  vf_in[0][6] = vf_in[0][7] = 0.0;

  // Adjust the volume fractions of material 0 in the state vector (we
  // don't have rem_cells_from_material implemented yet)
  statewrapper.mat_get_celldata("volfrac", 0, &vfmatptr);
  vfmatptr[statewrapper.cell_index_in_material(3, 0)] = 0.5;
  vfmatptr[statewrapper.cell_index_in_material(4, 0)] = 0.25;
  vfmatptr[statewrapper.cell_index_in_material(6, 0)] = 0.0;
  vfmatptr[statewrapper.cell_index_in_material(7, 0)] = 0.0;

  // Since the new material that displaced the old material has the
  // same density, the expected cell densities should be the same as
  // before. Verify against values from state vectors

  for (int c = 0; c < ncells; c++) rhocell[c] = 0.0;
  for (int m = 0; m < nmats; m++) {
    double const *rhomatptr;
    statewrapper.mat_get_celldata("mat_density", m, &rhomatptr);
    double const *vfmatptr;
    statewrapper.mat_get_celldata("volfrac", m, &vfmatptr);
    std::vector<int> matcells;
    statewrapper.mat_get_cells(m, &matcells);


    int nmatcells = matcells.size();
    for (int i = 0; i < nmatcells; i++) {
      int c = matcells[i];
      rhocell[c] += rhomatptr[i]*vfmatptr[i];
    }
  }
      
  for (int c = 0; c < ncells; c++)
    ASSERT_NEAR(rhocell_exp[c], rhocell[c], 1.0e-12);

}  

  

