/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include "portage/driver/driver.h"
#include "portage/interpolate/interpolate_1st_order.h"

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"
#include "JaliStateVector.h"

#include "portage/support/portage.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/wonton/state/simple_state/simple_state_wrapper.h"
#include "portage/interpolate/test/simple_intersect_for_tests.h"

double TOL = 1e-12;


/// First order interpolation of constant cell-centered field in 2D

TEST(Interpolate_1st_Order, Cell_Ctr_Const_2D) {

  // Create simple meshes
  
  std::shared_ptr<Portage::Simple_Mesh> source_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Portage::Simple_Mesh> target_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 5, 5);

  // Create mesh wrappers
  
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells
  
  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Portage::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source, 1.25);
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));
  
  // Create state wrapper
  
  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes 
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<2>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<2>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points
  
  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);
      
  // Loop over target cells
      
  for (int c = 0; c < ncells_target; ++c) {
  
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects
		
    BOX_INTERSECT::intersection_moments<2>(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    // Pack the results into a vector of true Portage::Weights_t
    
    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    
    // Put the weights in final form
    
    sources_and_weights[c] = wtsvec;
  }
  
  // Now do it the Portage way
  
  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Portage::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

      
  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL), 
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c) ASSERT_DOUBLE_EQ(stdval, outvals[c]);
  
}


/// First order interpolation of linear cell-centered field in 2D

TEST(Interpolate_1st_Order, Cell_Ctr_Lin_2D) {

  // Create simple meshes
  
  std::shared_ptr<Portage::Simple_Mesh> source_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Portage::Simple_Mesh> target_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 2, 2);

  // Create mesh wrappers
  
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells
  
  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Portage::Simple_State source_state(source_mesh);

// Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    Portage::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));  

  // Create state wrapper
  
  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes 
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<2>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<2>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points
  
  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);
      
  // Loop over target cells
      
  for (int c = 0; c < ncells_target; ++c) {
  
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects
		
    BOX_INTERSECT::intersection_moments<2>(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    // Pack the results into a vector of true Portage::Weights_t
    
    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    
    // Put the weights in final form
    
    sources_and_weights[c] = wtsvec;
  }
  
  // Now do it the Portage way
  
  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Portage::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

      
  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL), 
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);
  
  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATED CORRECTLY
  // IN THIS TEST

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    Portage::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    stdvals[c] = cen[0]+cen[1];
  }

  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdvals[c], outvals[c]);
  
}

/*
/// First order interpolation of constant node-centered field in 2D

TEST(Interpolate_1st_Order, Node_Ctr_Const_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  const int nnodes_source =
      source_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int nnodes_target =
      target_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);


  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh);

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", source_mesh,
                                  Jali::Entity_kind::NODE,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  // Create Interpolation objects

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(source_state);

  // Gather the dual cell coordinates for source and target meshes for
  // intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_dualcell_coords(nnodes_target);

  // Because the meshes are rectangular we can get away with examining
  // the coordinates of the corners instead of the wedges

  // Also, because we will use only the bounding box of each cell to
  // do the search and intersection we can get away with adding all
  // the coordinates of the corners to list including duplicates

  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(nnodes_target);

  for (int n = 0; n < nnodes_source; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);


    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[n] = wtsvec;
  }

  // Interpolate from source to target mesh

  std::vector<double> outvals(nnodes_target);

  Portage::Interpolate_1stOrder<2, Portage::NODE,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("nodevars");

  Jali::Entity_ID_List const& targetnodes =
      target_mesh->nodes<Jali::Entity_type::ALL>();

  Portage::transform(targetnodes.begin(), targetnodes.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_NEAR(stdval, outvals[n], TOL);
}

*/

/// First order interpolation of constant cell-centered field in 3D

TEST(Interpolate_1st_Order, Cell_Ctr_Const_3D) {
 // Create simple meshes
  
  std::shared_ptr<Portage::Simple_Mesh> source_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  std::shared_ptr<Portage::Simple_Mesh> target_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  // Create mesh wrappers
  
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells
  
  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Portage::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source, 1.25);
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));
  
  // Create state wrapper
  
  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes 
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<3>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<3>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points
  
  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);
      
  // Loop over target cells
      
  for (int c = 0; c < ncells_target; ++c) {
  
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects
		
    BOX_INTERSECT::intersection_moments<3>(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    // Pack the results into a vector of true Portage::Weights_t
    
    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    
    // Put the weights in final form
    
    sources_and_weights[c] = wtsvec;
  }
  
  // Now do it the Portage way
  
  // Create Interpolation object

  Portage::Interpolate_1stOrder<3, Portage::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

      
  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL), 
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c) ASSERT_DOUBLE_EQ(stdval, outvals[c]);


}

/// First order interpolation of linear cell-centered field in 3D

TEST(Interpolate_1st_Order, Cell_Ctr_Lin_3D) {
  // Create simple meshes
  
  std::shared_ptr<Portage::Simple_Mesh> source_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  std::shared_ptr<Portage::Simple_Mesh> target_mesh =
    std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  // Create mesh wrappers
  
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells
  
  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Portage::Simple_State source_state(source_mesh);

// Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    Portage::Point<3> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));  

  // Create state wrapper
  
  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes 
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<3>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<3>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points
  
  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);
      
  // Loop over target cells
      
  for (int c = 0; c < ncells_target; ++c) {
  
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects
		
    BOX_INTERSECT::intersection_moments<3>(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    // Pack the results into a vector of true Portage::Weights_t
    
    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    
    // Put the weights in final form
    
    sources_and_weights[c] = wtsvec;
  }
  
  // Now do it the Portage way
  
  // Create Interpolation object

  Portage::Interpolate_1stOrder<3, Portage::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

      
  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL), 
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);
  
  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATED CORRECTLY
  // IN THIS TEST

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    Portage::Point<3> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    stdvals[c] = cen[0]+cen[1];
  }

  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdvals[c], outvals[c]);

}


/*

/// First order interpolation of constant node-centered field in 3D

TEST(Interpolate_1st_Order, Node_Ctr_Const_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);

  const int nnodes_source =
      source_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int nnodes_target =
      target_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh);

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", source_mesh,
                                  Jali::Entity_kind::NODE,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(source_state);


  // Interpolate from source to target mesh

  Jali::Entity_ID_List const& targetnodes =
      target_mesh->nodes<Jali::Entity_type::ALL>();

  std::vector<double> outvals(nnodes_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(nnodes_target);

  // Gather the dual cell coordinates for source and target meshes for
  // intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_dualcell_coords(nnodes_target);

  // Because the meshes are rectangular we can get away with examining
  // the coordinates of the corners instead of the wedges

  // Also, because we will use only the bounding box of each cell to
  // do the search and intersection we can get away with adding all
  // the coordinates of the corners to list including duplicates

  for (int n = 0; n < nnodes_source; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);


    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[n] = wtsvec;
  }

  // Create Interpolation object

  Portage::Interpolate_1stOrder<3, Portage::NODE,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("nodevars");

  Portage::transform(targetnodes.begin(), targetnodes.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_NEAR(stdval, outvals[n], TOL);
}
*/
