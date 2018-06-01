/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <sys/time.h>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <memory>

#include "portage/simple_mesh/simple_mesh.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/state/state_manager.h"
#include "portage/driver/mmdriver.h"

#include "read_material_data.h"
#include "user_field.h"

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/driver/write_to_gmv.h"
#endif

#ifdef ENABLE_MPI
#include <mpi.h>
#endif


/*! 
  @file simple_mesh_mm_app.cc 

  @brief A simple application that remaps multi-material fields
  between two meshes - the meshes can be internally generated
  rectangular meshes or externally read unstructured meshes

  This program is used to showcase our capabilities with various types
  of remap operations (e.g. interpolation order) on various types of
  meshes (2d or 3d; node-centered or zone-centered) with multiple
  materials. The material data (material volume fractions, material
  centroids and field data) are specified in a file with a .mat
  extension

*/

//////////////////////////////////////////////////////////////////////
using namespace Portage;

int print_usage() {
  std::cout << std::endl;
  std::cout << "Usage: portageapp \n\n";
  
  std::cout << "  --nsourcecells (NO DEFAULT): Internally generated rectangular "
            << "SOURCE mesh with num cells in each coord dir\n\n";

  std::cout << "  --ntargetcells (NO DEFAULT): Internally generated rectangular "
            << "TARGET mesh with num cells in each coord dir\n\n";

  std::cout << "  --mesh_min (default = 0.): " <<
      "coordinates (same in x, y, and z) of the lower corner of a mesh\n\n";

  std::cout << "  --mesh_max (default = 1.): " <<
      "coordinates (same in x, y, and z) of the upper corner of a mesh\n\n" ;

  std::cout << "  --conformal (default = y): 'y' means mesh boundaries match\n" <<
      "    If 'n', target mesh is shifted in all directions by half-domain width\n\n";

  std::cout << "  --material_file=filename: file containing volume fractions\n" <<
      "    and, optionally, material centroids for cells. The file lists the\n" <<     
			"    total number of materials and a numeric value indicating\n" <<
      "    if only volume fractions are given (1) or volume fractions and\n" <<
      "    centroids are given (2). Subsequent lines list the cell ID,\n" <<
      "    number of materials in the cell, the material IDs,\n" <<
      "    the volume fractions and the centroids for the materials.\n\n";

  std::cout << "  --material_fields: A comma separated list of quoted math expressions \n" <<
      "    expressed in terms of x, y and z following the syntax of \n" <<
      "    the expression parser package ExprTk \n" <<
      "    (http://www.partow.net/programming/exprtk/)\n"<<
  		"    The syntax is generally predictable, e.g. \"24*x + y*y\" or \"x + log(abs(y)+1)\"\n"<<
  		"    There must be as many expressions as there materials in the problem\n\n";

  std::cout << "--remap order (default = 1): " <<
      "order of accuracy of interpolation\n\n";

  std::cout << "--limiter (default = 0): " <<
      "slope limiter for a piecewise linear reconstrution\n\n";

  std::cout << "--output_meshes (default = y)\n";
  std::cout << "  If 'y', the source and target meshes are output with the " <<
      "remapped field attached as input.exo and output.exo. \n\n";

  std::cout << "--results_file=results_filename (default = output.txt)\n";
  std::cout << "  If a filename is specified, the target field values are " <<
      "output to the file given by 'results_filename' in ascii format\n\n";
  return 0;

}

// Run a remap between two meshes and return the L1 and L2 error norms
// with respect to the specified field. If a field was not specified and
// remap only volume fractions (and if specified, centroids)
void run(
	const Simple_Mesh& sourceMesh,
  const Simple_Mesh& targetMesh,
  const std::string material_filename,
  const std::vector<std::string> material_field_expressions,
  Portage::LimiterType limiter,
  int interp_order,
  std::string field_filename, 
  bool mesh_output,
  Portage::Entity_kind entityKind,
  double *L1_error,
  double *L2_error
	) {

  std::cout << "starting simple_mesh_mm_app...\n";

  // Wrappers for interfacing with the underlying mesh data structures.
  Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sourceMesh);
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(targetMesh);

  const int nsrccells = source_mesh_wrapper.num_owned_cells() +
      source_mesh_wrapper.num_ghost_cells(); // I don't believe there are ghosts
  const int ntarcells = target_mesh_wrapper.num_owned_cells();

  const int nsrcnodes = source_mesh_wrapper.num_owned_nodes() +
      source_mesh_wrapper.num_ghost_nodes(); // I don't believe there are ghosts
  const int ntarnodes = target_mesh_wrapper.num_owned_nodes();
  
	// declare local variables for the material data
  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;  // flattened 2D array
  std::vector<double> cell_mat_volfracs;  // flattened 2D array
  std::vector<Portage::Point<2>> cell_mat_centroids;  // flattened 2D array
 
  // Read volume fraction and centroid data from file 
  read_material_data<Wonton::Simple_Mesh_Wrapper, 2>(source_mesh_wrapper,
  	material_filename, cell_num_mats, cell_mat_ids, cell_mat_volfracs,
    cell_mat_centroids);
  
	// debug info, easy to remove
  if (0){
		std::cout<<"\n\nSource Mesh: \n"<<std::endl;
		std::cout<<"cell_num_mats: ";
		for (auto x: cell_num_mats) std::cout<<x<<" "; std::cout<<std::endl;
		                                                   
		std::cout<<"cell_mat_ids: ";
		for (auto x: cell_mat_ids) std::cout<<x<<" "; std::cout<<std::endl;
		                                                   
		std::cout<<"cell_mat_volfracs: ";
		for (auto x: cell_mat_volfracs) std::cout<<x<<" "; std::cout<<std::endl;
		                                                   
		std::cout<<"cell_mat_centroids:\n";
		for (auto x: cell_mat_centroids) std::cout<<x[0]<<" "<<x[1]<<"\n"; std::cout<<std::endl;
  }                                                  
  
  bool mat_centroids_given = (cell_mat_centroids.size() > 0);

  // Compute offsets into flattened arrays based on cell_num_mats 
  std::vector<int> offsets(nsrccells);
  offsets[0] = 0;
  for (int i = 1; i < nsrccells; i++)
    offsets[i] = offsets[i-1] + cell_num_mats[i-1];


  // Count the number of materials in the problem and gather their IDs  
	std::unordered_set<int> mat_ids{cell_mat_ids.begin(), cell_mat_ids.end()};
  int nmats{mat_ids.size()};


  // Spit out some information for the user
  
  std::cout << "Source mesh has " << nsrccells << " cells\n";
  std::cout << "Target mesh has " << ntarcells << " cells\n";
  if (material_field_expressions.size()) {
    std::cout << " Specified fields for materials are ";
    for (auto const& expr : material_field_expressions)
      std::cout << "    " << expr << ", ";
    std::cout << "\n";
    if (material_field_expressions.size() < nmats)
      std::cout << "Not all material fields are specified. Missing ones will be set to 0\n";
  }
  std::cout << "   Interpolation order is " << interp_order << "\n";
  if (interp_order == 2)
    std::cout << "   Limiter type is " << limiter << "\n";



  // Perform interface reconstruction for pretty pictures (optional)
  //auto interface_reconstructor =
  //    std::make_shared<Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
   //                                    Wonton::Simple_Mesh_Wrapper>>(source_mesh_wrapper);

	Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,Wonton::Simple_Mesh_Wrapper> interface_reconstructor{source_mesh_wrapper};
	
  // convert from Portage points to Tangram points
  std::vector<Tangram::Point<2>> 
  	Tcell_mat_centroids(cell_mat_centroids.begin(), cell_mat_centroids.end());
  	
  
  interface_reconstructor.set_volume_fractions(cell_num_mats,
                                                cell_mat_ids,
                                                cell_mat_volfracs,
                                                Tcell_mat_centroids);
                                                
  // do the interface reconstruction (go from volume fraction, centroid) to pure
  // material polygons
  interface_reconstructor.reconstruct();

  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> const&
      cellmatpoly_list = interface_reconstructor.cell_matpoly_ptrs();

  Tangram::write_to_gmv<Wonton::Simple_Mesh_Wrapper, 2>(source_mesh_wrapper,
                                                        nmats,
                                                        cell_num_mats,
                                                        cell_mat_ids,
                                                        cellmatpoly_list,
                                                        "source_ir.gmv");

  // Convert data from cell-centric to material-centric form as we
  // will need it for adding it to the state manager
  std::unordered_map<std::string,int> matnames;
  std::unordered_map<int, std::vector<int>> matcells;
  std::unordered_map<int, std::vector<double>> mat_volfracs;
  std::unordered_map<int, std::vector<Point<2>>>mat_centroids;
  for (int c = 0; c < nsrccells; c++) {
    int ibeg = offsets[c];
    for (int j = 0; j < cell_num_mats[c]; j++) {
      int m = cell_mat_ids[ibeg+j];
      matnames["mat"+std::to_string(m)]=m; //overwrites but no harm
      matcells[m].push_back(c);
      mat_volfracs[m].push_back(cell_mat_volfracs[ibeg+j]);
      if (mat_centroids_given)
        mat_centroids[m].push_back(cell_mat_centroids[ibeg+j]);
    }
  }
  
  // Native state managers for source and target
  // Note unlike portageapp_multimat_jali and others, I do not need to use a
  // wrapper since the statemanager is feature complete
	StateManager<Wonton::Simple_Mesh_Wrapper> source_state{
		source_mesh_wrapper, matnames, matcells};

	// add the volume fractions and centroids to the state manager
	source_state.add(std::make_shared<StateVectorMulti<>>("mat_volfracs",mat_volfracs));
	if (mat_centroids_given) 
		source_state.add(std::make_shared<StateVectorMulti<Point<2>>> ("mat_centroids", mat_centroids));
 
  // User specified fields on source
  std::vector<user_field_t> mat_fields(nmats);
  
  // only do this block if there are any user specified fields
  if (material_field_expressions.size()) {
  
   	// the mm state
  	std::unordered_map<int,std::vector<double>> user_field;
  
 		// loop over materials (there needs to be one field per material
    for (int m = 0; m < nmats; m++) {
    
    	// if we can't create the user field then die
      if (!mat_fields[m].initialize(2, material_field_expressions[m]))
      	throw std::runtime_error("Could not initialize user field: "
      		+material_field_expressions[m]);
      
      // number of cells for this material
      int nmatcells = matcells[m].size();
      
      // field values for this material
      std::vector<double> matData;
      
      // loop of cell indices for this material
      for (int c: matcells[m]) {
        Point<2> point;
        source_mesh_wrapper.cell_centroid(c, &point);
        matData.push_back(mat_fields[m](point));
      }
      
      // add the data to the user_field
      user_field[m]=matData;
    }
    
    // print block
    //std::cout<<"\nUser field values by material: "<<std::endl;
    //for (auto& kv: user_field){
    //	std::cout<<kv.first<<": ";
    //	for (auto x:kv.second) std::cout<<x<<" ";
    //	std::cout<<std::endl;
    //}
		
		// add the field to the state manager
		source_state.add(std::make_shared<StateVectorMulti<>>("cellmatdata",user_field));
		
		// from the state manager
    //std::cout<<"\nState Material User field values by material: "<<std::endl;
    //auto &cellmatdata= source_state.get<double>("cellmatdata")->get_data();
    //for (auto& kv: cellmatdata){
    //	std::cout<<kv.first<<": ";
    //	for (auto x:kv.second) std::cout<<x<<" ";
    //	std::cout<<std::endl;
    //}
		
		
  }

	// print what the state manager knows about
  std::cout<<"\nState manager fields: ";
  for (auto x: source_state.get_state_keys()) std::cout<<x.second<<" ";
  std::cout<<std::endl;

  // Add the materials into the target mesh without cell indices
  // The remap algorithm will figure out which cells contain which materials
	StateManager<Wonton::Simple_Mesh_Wrapper> target_state{target_mesh_wrapper, matnames};

	// Add the volume fractions, centroids and cellmatdata variables
  target_state.add(
  	std::make_shared<StateVectorMulti<>>(StateVectorMulti<>{"mat_volfracs"}));
  target_state.add(
  	std::make_shared<StateVectorMulti<Point<2>>>(StateVectorMulti<Point<2>>{"mat_centroids"}));

  // Register the variable name and interpolation order with the driver
  std::vector<std::string> remap_fields;

  // If the user specified some material fields, then add a placeholder for
  // them on the target side
  if (material_field_expressions.size()) {
  	target_state.add(
  		std::make_shared<StateVectorMulti<>>(StateVectorMulti<>{"cellmatdata"}));
    remap_fields.push_back("cellmatdata");
  } 

	// print what the state manager knows about
  //std::cout<<"\nTarget state manager fields: ";
  //for (auto x: target_state.get_state_keys()) std::cout<<x.second<<" ";
  //std::cout<<std::endl<<std::endl;


  if (interp_order == 1) {
    Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_1stOrder,
      2,
      Wonton::Simple_Mesh_Wrapper,
      StateManager<Wonton::Simple_Mesh_Wrapper>,
      Wonton::Simple_Mesh_Wrapper,
      StateManager<Wonton::Simple_Mesh_Wrapper>,
      Tangram::XMOF2D_Wrapper>
        d(source_mesh_wrapper, source_state,
          target_mesh_wrapper, target_state);
    d.set_remap_var_names(remap_fields);
    d.run(false);
  } else if (interp_order == 2) {
    Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_2ndOrder,
      2,
      Wonton::Simple_Mesh_Wrapper,
      StateManager<Wonton::Simple_Mesh_Wrapper>,
      Wonton::Simple_Mesh_Wrapper,
      StateManager<Wonton::Simple_Mesh_Wrapper>,
      Tangram::XMOF2D_Wrapper>
        d(source_mesh_wrapper, source_state,
          target_mesh_wrapper, target_state);
    d.set_remap_var_names(remap_fields, limiter);
    d.run(false);
  }
 
 	// make sure the volume fraction and centroid are available on the target mesh   

	// material dominant list of cells for each material
	auto target_material_cells=target_state.get_material_cells();
	
	// material dominant volume fractions and centroids
	std::unordered_map<int,vector<int>> map_target_cell_materials;
	std::unordered_map<int,vector<double>> map_target_cell_mat_volfracs;
	std::unordered_map<int,vector<Tangram::Point<2>>> map_target_cell_mat_centroids;
	
	// To calculate the cell dominant volume factions and cell centroids,
	// instead of using the reverse map of cell materials we go in the forward
	// direction again because we need to do the same thing for volume fractions and centroids.
	// If we used the cell materials, we would need a large number of finds. Going in the forward
	// direction, eliminates this at the cost of recomputing cell_mats, we get the
	// volume fraction and centroids for free.
	
	// loop over target material cells
	for (auto& kv:target_material_cells){
	
		// get the material
		const int m{kv.first};
		
		// get the vector of cell ids
		const vector<int> cells{kv.second};
		
		// get the volume fractions for this material
		//const std::shared_ptr<StateVectorMulti<double>> pvolfracs= target_state.get("mat_volfracs");
		//const auto& volfracs = pvolfracs->get_data().at(m);
		//const auto& volfracs = target_state.get<double>("mat_volfracs")->get_data()[m];
		
		// get the material centroids for this material
		const auto& volfracs = target_state.get<double>("mat_volfracs")->get_data()[m];		
		const auto& centroids = target_state.get<Point<2>>("mat_centroids")->get_data()[m];		
	
		// loop over the cell id's
		for (int i = 0; i<cells.size(); ++i){
		
			// get this cell id
			const int c=cells[i];
			
			// push the volume fraction and centroid data onto the maps
			map_target_cell_materials[c].push_back(m);
			map_target_cell_mat_volfracs[c].push_back(volfracs[i]);
			map_target_cell_mat_centroids[c].push_back(centroids[i]);
			
		}
	}
		
	// declare what Tangram needs (flat arrays
  std::vector<int> target_cell_num_mats;
  std::vector<int> target_cell_mat_ids;
  std::vector<double> target_cell_mat_volfracs;
  std::vector<Tangram::Point<2>> target_cell_mat_centroids;

	for (int c = 0; c<ntarcells; ++c){
		
		// get the vector of material ids,volume fraction and centroid
		// with thrust turned on, these need to be portage vectors, not std::vectors
		const Portage::vector<int> mats{map_target_cell_materials.at(c)};
		const Portage::vector<double> volfracs{map_target_cell_mat_volfracs.at(c)};
		const Portage::vector<Tangram::Point<2>> centroids{map_target_cell_mat_centroids.at(c)};
		
		// push the size onto the number of mats
		target_cell_num_mats.push_back(mats.size());
		
		// loop over the matid's
		for (int i=0; i<mats.size(); ++i){
		
			// get the material id
			int matid=mats[i];
			
			// push the material id, volume fraction and centroid
			target_cell_mat_ids.push_back(matid);
			target_cell_mat_volfracs.push_back(volfracs[i]);
			target_cell_mat_centroids.push_back(centroids[i]);
		}

	}

	// debug info, easy to remove
	if (0){
		std::cout<<"\n\nTarget Mesh: \n"<<std::endl;
		std::cout<<"target_cell_num_mats: ";
		for (auto x: target_cell_num_mats) std::cout<<x<<" "; std::cout<<std::endl;
		                                                   
		std::cout<<"target_cell_mat_ids: ";
		for (auto x: target_cell_mat_ids) std::cout<<x<<" "; std::cout<<std::endl;
		                                                   
		std::cout<<"target_cell_mat_volfracs: ";
		for (auto x: target_cell_mat_volfracs) std::cout<<x<<" "; std::cout<<std::endl;
		                                                   
		std::cout<<"target_cell_mat_centroids:\n";
		for (Tangram::Point<2>  x: target_cell_mat_centroids) std::cout<<x[0]<<" "<<x[1]<<"\n"; std::cout<<std::endl;
	}
	
  // Perform interface reconstruction on target mesh for pretty pictures
  // (optional)

  
  auto target_interface_reconstructor =
      std::make_shared<Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
                                       Wonton::Simple_Mesh_Wrapper>>(target_mesh_wrapper);

  target_interface_reconstructor->set_volume_fractions(target_cell_num_mats,
                                                target_cell_mat_ids,
                                                target_cell_mat_volfracs,
                                                target_cell_mat_centroids);
  target_interface_reconstructor->reconstruct();

  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> const&
      target_cellmatpoly_list = target_interface_reconstructor->cell_matpoly_ptrs();

  Tangram::write_to_gmv<Wonton::Simple_Mesh_Wrapper, 2>(target_mesh_wrapper,
                                                      nmats,
                                                      target_cell_num_mats,
                                                      target_cell_mat_ids,
                                                      target_cellmatpoly_list,
                                                      "target_ir.gmv");



  
  // Compute error

  if (material_field_expressions.size()) {
    double error, toterr = 0.0;
    double const * cellmatvals;
    double totvolume = 0.;
    for (int m = 0; m < nmats; m++) {
      target_state.mat_get_celldata<double>("cellmatdata", m, &cellmatvals);
      
      std::vector<int> matcells;
      target_state.mat_get_cells(m, &matcells);
      
      // Cell error computation
      Portage::Point<2> ccen;
      int nmatcells = matcells.size();
      for (int ic = 0; ic < nmatcells; ++ic) {
        int c = matcells[ic];
        target_mesh_wrapper.cell_centroid(c, &ccen);
        error = mat_fields[m](ccen) - cellmatvals[ic];
        
        if (!target_mesh_wrapper.on_exterior_boundary(Portage::Entity_kind::CELL, c)) {
          double cellvol = target_mesh_wrapper.cell_volume(c);
          totvolume += cellvol;
          *L1_error += fabs(error)*cellvol;
          *L2_error += error*error*cellvol;
        }
      }
    }
  }

  *L2_error = sqrt(*L2_error);
  std::printf("\n\nL1 NORM OF ERROR = %lf\n", *L1_error);
  std::printf("L2 NORM OF ERROR = %lf\n\n", *L2_error);
 
  // construct the field file name and open the file
  
  if (field_filename.length()) {
    
    std::ofstream fout(field_filename);
    fout << std::scientific;
    fout.precision(17);
    fout <<"User Field Remap: ";
    for (auto& exp : material_field_expressions) fout<<exp<< " ";
    fout << std::endl;
    
    // write out the values
    for (int m=0; m<nmats; ++m){
    	std::vector<int>cells = target_state.get_material_cells().at(m);
    	std::vector<double> data = target_state.get("cellmatdata")->get_data().at(m);
    	fout << "material: " << m << std::endl;
    	
    	for (int i=0; i<data.size(); ++i){
    		fout<<"  cell "<<cells[i]<<": "<<data[i]<<std::endl;
    	}
    }
  }
}

int main(int argc, char** argv) {

	#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);
	#endif


	// print usage if no args
  if (argc == 1) {
  	print_usage();
  	return 0;
  }
  
  struct timeval begin, end, diff;

	// declare simulation parameters
  int nsourcecells = 0, ntargetcells = 0;  // No default
  int dim = 2;
  double srclo = 0.0, srchi = 1.0;  // bounds of generated mesh in each dir
  bool conformal = true;
  std::string material_filename;
  std::vector<std::string> material_field_expressions;
  int interp_order = 1;
  bool mesh_output = true;
  int n_converge = 1;
  Portage::LimiterType limiter = Portage::LimiterType::NOLIMITER;
  Portage::Entity_kind entityKind = Portage::Entity_kind::CELL;
  double L1_error=0., L2_error=0.;

  std::string field_output_filename;  

  // Parse the input

  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    std::size_t len = arg.length();
    std::size_t keyword_beg = 2;
    std::size_t keyword_end = arg.find_first_of("=");
    std::string keyword = arg.substr(keyword_beg, keyword_end-keyword_beg);
    std::string valueword = arg.substr(keyword_end+1, len-(keyword_end+1));
    
    if (keyword == "nsourcecells")
      nsourcecells = stoi(valueword);
    else if (keyword == "ntargetcells")
      ntargetcells = stoi(valueword);
    else if (keyword == "mesh_min") {
      srclo = stof(valueword);
    } else if (keyword == "mesh_max") {
      srchi = stof(valueword);
    } else if (keyword == "conformal")
      conformal = (valueword == "y");
    else if (keyword == "material_file") {
      material_filename = valueword;
    }else if (keyword == "remap_order") {
      interp_order = stoi(valueword);
      assert(interp_order > 0 && interp_order < 3);
    }  else if (keyword == "material_fields") {
      // Expecting comma-separated quoted expressions
      std::string exprlist = valueword;
      std::size_t expr_beg = 0;
      std::size_t expr_end = 0;
      bool exprdone = false;
      std::string expr;
      while (!exprdone) {
        expr_end = exprlist.find_first_of(",", expr_beg);
        if (expr_end == std::string::npos) {
          exprdone = true;
          expr = exprlist.substr(expr_beg);
        } else {
          expr = exprlist.substr(expr_beg, expr_end-expr_beg);
          expr_beg = expr_end+1;
        }
        material_field_expressions.push_back(expr);
      }
    } else if (keyword == "limiter") {
      if (valueword == "barth_jespersen" || valueword == "BARTH_JESPERSEN")
        limiter = Portage::LimiterType::BARTH_JESPERSEN;
    } else if (keyword == "output_meshes") {
      mesh_output = (valueword == "y");
    } else if (keyword == "results_file") {
      field_output_filename = valueword;
    } else if (keyword == "convergence_study") {
      n_converge = stoi(valueword);
      if (n_converge <= 0) {
        std::cerr << "Number of meshes for convergence study should be greater than 0" << std::endl;
        throw std::exception();
      }
    } else if (keyword == "help") {
      print_usage();
    } else
      std::cerr << "Unrecognized option " << keyword << std::endl;
  }

  gettimeofday(&begin, 0);

	double mesh_offset = conformal ? 0. : (srchi-srclo)/ntargetcells/2.;
	std::cout<<"mesh_offset="<<mesh_offset<<std::endl;
	// define the source and target meshes 
  Simple_Mesh source_mesh{srclo, srclo, srchi, srchi, nsourcecells, nsourcecells};
  Simple_Mesh target_mesh{
  	srclo + mesh_offset, srclo+ mesh_offset, 
  	srchi + mesh_offset, srchi+ mesh_offset, 
  	ntargetcells, ntargetcells};
    
  // Now run the remap on the meshes and get back the L2 error
  run(
  	source_mesh, 
  	target_mesh, 
  	material_filename, 
  	material_field_expressions,         
  	limiter, 
  	interp_order, 
    field_output_filename,
    mesh_output, 
    entityKind,
    &L1_error,
    &L2_error
	);
	

  gettimeofday(&end, 0);
  timersub(&end, &begin, &diff);
  float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
  std::cout << "Remap Time: " << seconds << std::endl;
  
  #ifdef ENABLE_MPI
  MPI_Finalize();
	#endif

}


