/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_INTERPOLATE_GRADIENT_H_
#define SRC_INTERPOLATE_GRADIENT_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/Matrix.h"
#include "portage/support/lsfits.h"

#include "tangram/driver/driver.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"

#include "portage/interpolate/ir_dummy.h"

namespace Portage {

#define LIMITED_GRADIENT_BASE(ENTITY_KIND)				\
  Limited_Gradient_Base <D, ENTITY_KIND, MeshType, StateType, InterfaceReconstructorType>

#define LIMITED_GRADIENT_BASE_NIR(ENTITY_KIND)			\
  Limited_Gradient_Base <D, ENTITY_KIND, MeshType, StateType>

  template<int D, Entity_kind on_what, typename MeshType, typename StateType,
    template<class, int> class InterfaceReconstructorType =
    detail::DummyInterfaceReconstructor>
    class Limited_Gradient_Base {

  public:
 
#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
  Tangram::Driver<InterfaceReconstructorType, D, MeshType>;

  Limited_Gradient_Base(MeshType const & mesh, StateType const & state,
			std::shared_ptr<InterfaceReconstructor> ir) :
  Limited_Gradient_Base(mesh,state), interface_reconstructor_(ir) {}
#endif  

  Limited_Gradient_Base(MeshType const & mesh, StateType const & state)
  : mesh_(mesh),state_(state),vals_(nullptr) {
    assert(D == this->mesh_.space_dimension());
    assert(D == 2 || D == 3);
  }

  void set_material(int matid, std::string const var_name, LimiterType limiter_type) {
    matid_ = matid;
    var_name_ = var_name;
    limtype_ = limiter_type;
  }
 
  Limited_Gradient_Base & operator = (const Limited_Gradient_Base &) = delete;
  ~Limited_Gradient_Base() {}
 
  Vector<D> operator()(int entity_id) {
    std::cerr << "Limited gradient not implemented for this entity kind\n";
  }
 
  protected:
  LimiterType limtype_;
  std::string var_name_;
  int matid_;
  MeshType const & mesh_;
  StateType const & state_;
  double const *vals_;
  std::vector<int> cellids_;
  Field_type field_type_;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
  };

  ///////////////////////////////////////////////////////////////////////////////

  /*! @class Limited_Gradient gradient.h
    @brief Compute limited gradient of a field or components of a field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
    @tparam on_what An enum type which indicates different entity types
  */

  template<int D, Entity_kind on_what, typename MeshType,typename StateType,
    template<class, int> class InterfaceReconstructorType =
    detail::DummyInterfaceReconstructor>

    class Limited_Gradient : public LIMITED_GRADIENT_BASE(on_what) {
  public:

  /*! @brief Constructor
    @param[in] mesh  Mesh class than one can query for mesh info
    @param[in] state A state manager class that one can query for field info
    @param[in] on_what An enum that indicates what type of entity the field is on
    @param[in] var_name Name of field for which the gradient is to be computed
    @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)
   
    @todo must remove assumption that field is scalar
  */

#ifdef HAVE_TANGRAM
  Limited_Gradient(MeshType const & mesh, 
		   StateType const & state,
		   std::string const var_name,  // TODO: remove
		   LimiterType limiter_type,    // TODO: remove
		   std::shared_ptr<typename LIMITED_GRADIENT_BASE(on_what)::InterfaceReconstructor> ir)
  : LIMITED_GRADIENT_BASE(on_what) (mesh, state, ir) {}
#endif

  Limited_Gradient(MeshType const & mesh, 
		   StateType const & state,   
		   std::string const var_name, // TODO: remove
		   LimiterType limiter_type)   // TODO: remove
  : LIMITED_GRADIENT_BASE(on_what) (mesh, state) {}
  };

  ///////////////////////////////////////////////////////////////////////////////

  /*! @class Limited_Gradient<MeshType,StateType,NODE> gradient.h
    @brief Specialization of limited gradient class for @c node-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
  */

#define LIMITED_GRADIENT_NIR(ENTITY_KIND)			\
  Limited_Gradient <D, ENTITY_KIND, MeshType, StateType>
 
  template<int D, typename MeshType, typename StateType>
    class Limited_Gradient<D, NODE, MeshType, StateType>
    : public LIMITED_GRADIENT_BASE_NIR(NODE) {

 public:

  /*! @brief Constructor
    @param[in] mesh  Mesh class than one can query for mesh info
    @param[in] state A state manager class that one can query for field info
    @param[in] var_name Name of field for which the gradient is to be computed
    @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)

    @todo must remove assumption that field is scalar
  */
   
  Limited_Gradient(MeshType const & mesh, StateType const & state,
		   std::string const var_name,
		   LimiterType limiter_type)
    : LIMITED_GRADIENT_BASE_NIR(NODE)(mesh, state) {
      int nnodes = this->mesh_.num_entities(NODE);
      node_neighbors_.resize(nnodes);
      Portage::for_each(this->mesh_.begin(NODE), this->mesh_.end(NODE),
			[this](int n) { this->mesh_.dual_cell_get_node_adj_cells(
					n, ALL, &(node_neighbors_[n])); } );
    }

    /* Collect and keep the list of neighbors for each CELL as it may
       be expensive to go to the mesh layer and collect this data for
       each cell during the actual gradient calculation */
    void set_material(int matid,
	      std::string const var_name, 
	      LimiterType limiter_type) {
      this->set_material(matid,var_name,limiter_type);
      this->state_.mesh_get_data(NODE, this->var_name_, &this->vals_);
    }

    Vector<D> operator() (int nodeid);
   
  private:
    std::vector<std::vector<int>> node_neighbors_;
  };

  // @brief Limited gradient functor implementation for NODE

  template<int D, typename MeshType, typename StateType>
    Vector<D> LIMITED_GRADIENT_NIR(NODE)::operator() (int nodeid) {

    assert(this->vals_);

    double phi = 1.0;
    Vector<D> grad;

    std::vector<int> const & nbrids = node_neighbors_[nodeid];  
    std::vector<Point<D>> nodecoords(nbrids.size()+1);
    std::vector<double> nodevalues(nbrids.size()+1);
    this->mesh_.node_get_coordinates(nodeid, &(nodecoords[0]));
    nodevalues[0] = this->vals_[nodeid];
  
    int i = 1;
    for (auto const & nbrnode : nbrids) {
      this->mesh_.node_get_coordinates(nbrnode, &nodecoords[i]);
      nodevalues[i] = this->vals_[nbrnode];
      i++;
    }
  
    grad = ls_gradient(nodecoords, nodevalues);

    if (this->limtype_ == BARTH_JESPERSEN &&
	!this->mesh_.on_exterior_boundary(NODE, nodeid)) {  // No limiting on boundary
    
      // Min and max vals of function (cell centered vals) among neighbors
      double minval = this->vals_[nodeid];
      double maxval = this->vals_[nodeid];
    
      for (auto const & val : nodevalues) {
	minval = std::min(val, minval);
	maxval = std::max(val, maxval);
      }
    
      // Find the min and max of the reconstructed function in the cell
      // Since the reconstruction is linear, this will occur at one of
      // the nodes of the cell. So find the values of the reconstructed
      // function at the nodes of the cell
    
      double nodeval = this->vals_[nodeid];
    
      std::vector<Point<D>> dualcellcoords;
      this->mesh_.dual_cell_get_coordinates(nodeid, &dualcellcoords);
    
      for (auto const & coord : dualcellcoords) {
	Vector<D> vec = coord-nodecoords[0];
	double diff = dot(grad, vec);
	double extremeval = (diff > 0.0) ? maxval : minval;
	double phi_new = (diff == 0.0) ? 1 : (extremeval-nodeval)/diff;
	phi = std::min(phi_new, phi);
      }
    }

    return phi*grad;
  }

  ///////////////////////////////////////////////////////////////////////////////

  /*! @class Limited_Gradient<MeshType,StateType,CELL> gradient.h
    @brief Specialization of limited gradient class for @c cell-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
  */

#define LIMITED_GRADIENT(ENTITY_KIND)					\
  Limited_Gradient <D, ENTITY_KIND, MeshType, StateType, InterfaceReconstructorType>
 
  template<int D, typename MeshType, typename StateType,
    template<class, int> class InterfaceReconstructorType>
    class Limited_Gradient<D, CELL, MeshType, StateType, InterfaceReconstructorType>
    : public LIMITED_GRADIENT_BASE(CELL) {
 
 public:

#ifdef HAVE_TANGRAM 
  Limited_Gradient(MeshType const & mesh, 
		   StateType const & state, 
		   std::string const var_name, // TODO: remove
		   LimiterType limiter_type,   // TODO: remove
		   std::shared_ptr<typename LIMITED_GRADIENT_BASE(CELL)::InterfaceReconstructor> ir) 
    : LIMITED_GRADIENT_BASE(CELL) (mesh, state, ir) { 
      this->var_name_=var_name;      // TODO: remove
      this->limtype_=limiter_type;   // TODO: remove

      // Collect and keep the list of neighbors for each CELL as it may
      // be expensive to go to the mesh layer and collect this data for
      // each cell during the actual gradient calculation

      cell_neighbors_.resize(this->mesh_.num_entities(CELL));
      Portage::for_each(this->mesh_.begin(CELL), this->mesh_.end(CELL), 
			[this](int c) { this->mesh_.cell_get_node_adj_cells(
					 c, ALL, &(cell_neighbors_[c])); } );
    }
#endif
  
  Limited_Gradient(MeshType const & mesh, 
		   StateType const & state, 
		   std::string const var_name, // TODO: remove
		   LimiterType limiter_type)   // TODO: remove
    : LIMITED_GRADIENT_BASE_NIR(CELL) (mesh, state) { 
      this->var_name_=var_name;      // TODO: remove
      this->limtype_=limiter_type;   // TODO: remove
    }
  
    void set_material(int matid,
	      std::string const var_name, 
	      LimiterType limiter_type) {

      // Extract the field data from the statemanager
      this->field_type_ = this->state_.field_type(CELL, this->var_name_);
      if (this->field_type_ == Field_type::MESH_FIELD) {
	this->state_.mesh_get_data(CELL, this->var_name_, &this->vals_);
      } else {
	this->state_.mat_get_celldata(this->var_name_, this->matid_, &this->vals_);
      }
      assert(this->vals_);
    }

    Vector<D> operator() (int cellid);
  
  private:
    std::vector<std::vector<int>> cell_neighbors_;
  };

  // @brief Implementation of Limited_Gradient functor for CELLs

  template<int D, typename MeshType, typename StateType, 
    template<class, int> class InterfaceReconstructorType>
    Vector<D> LIMITED_GRADIENT(CELL)::operator()(int cellid) {

    double phi = 1.0;
    Vector<D> grad;

    std::vector<int> nbrids{cellid}; // Include cell where grad is needed as first element
    nbrids.insert(std::end(nbrids), 
		  std::begin(cell_neighbors_[cellid]), 
		  std::end(cell_neighbors_[cellid]));
    std::vector<Point<D>> ls_coords;
    std::vector<double> ls_vals;
    
    // Loop over cell where grad is needed and its neighboring cells
    for (auto nbrid_g : nbrids) {

      // Regardless of Tangram availablity, could have material or mesh data
#ifdef HAVE_TANGRAM
    
      // In order to access its field values for material data (non-mesh data), 
      // need local compressed neighbor cellid (nbrid_l); for mesh data,
      // local cellid (nbrid_l) and global cellid (nbrid_g) are equal
      int nbrid_l = (this->field_type_ == Field_type::MESH_FIELD) ?
	nbrid_g : this->state_.cell_index_in_material(nbrid_g, this->matid_);
      
      // In the case of material-data, need to check that neighbors contain material of interest
      // (i.e. have valid local id); in the case of mesh data, this is always true,
      // since local cellid (nbrid_l) is equal to global cellid (nbrid_g)
      if (nbrid_l >= 0) { // cell_index_in_material can return -1

	// Get cell's cellmatpoly
	Tangram::CellMatPoly<D> const& cellmatpoly =
	  this->interface_reconstructor_->cell_matpoly_data(nbrid_g);
	
	std::vector<int> cellmats;
	  this->state_.cell_get_mats(nbrid_g, &cellmats);
	if (cellmats.size() != 1){ // Multi-material cell
	  // Collect all the matpolys in this cell for the material of interest
	  std::vector<Tangram::MatPoly<D>> matpolys =
	    cellmatpoly.get_matpolys(this->matid_);
	    
	  // If there are multiple matpolys in this cell for the material of interest,
	  // aggregate moments to compute new centroid
	  for (int ipoly=0; ipoly<matpolys.size(); ipoly++) {
	    Portage::Point2 cell_centroid;
	    this->mesh_.cell_centroid(nbrid_g, &cell_centroid);
	    ls_coords.push_back(cellmatpoly.matpoly_centroid(ipoly));
	    break; // TODO: Instead of cutting out after the first matpoly,
	    // Get matpoly moments directly using new interface in Tangram,
	    // aggregate, and use to compute overall material centroid
	  }
	    
	  // Populate least squares vectors with centroid for material
	  // of interest and field value in the current cell for that material
	  ls_vals.push_back(this->vals_[nbrid_l]);
	    
	} else { // Single material cell

	  // Ensure that the single material is the material of interest
	  if (cellmats[0] == this->matid_) {
	    // Get the cell-centered value for this material
	    Point<D> centroid;
	    this->mesh_.cell_centroid(nbrid_g, &centroid);
	    ls_coords.push_back(centroid);
	    ls_vals.push_back(this->vals_[nbrid_l]);
	  }
	}
#else
	// If we get here, we must have mesh data which is cell-centered 
	// and not dependent on material, so just get the centroid and value
	Point<D> centroid;
	this->mesh_.cell_centroid(nbrid_g, &centroid);
	ls_coords.push_back(centroid);
	ls_vals.push_back(this->vals_[nbrid_g]);
#endif
      }
    }
    grad = ls_gradient(ls_coords, ls_vals);

    // Limit the gradient to enforce monotonicity preservation
    if (this->limtype_ == BARTH_JESPERSEN &&
	!this->mesh_.on_exterior_boundary(CELL, cellid)) {  // No limiting on boundary
      
      phi = 1.0;
      
      // Min and max vals of function (cell centered vals) among neighbors
      // and the cell itself
      /// @todo: must remove assumption the field is scalar
    
      double minval = ls_vals[0];
      double maxval = ls_vals[0]; 
      double cellcenval = ls_vals[0];

      // Find min and max values among all neighbors (exlude the first element 
      // in nbrids because it corresponds to the cell itself, not a neighbor)
      for (int i = 1; i < nbrids.size(); ++i) {
	minval = std::min(ls_vals[i], minval);
	maxval = std::max(ls_vals[i], maxval);
      }
      
      // Find the min and max of the reconstructed function in the cell
      // Since the reconstruction is linear, this will occur at one of
      // the nodes of the cell. So find the values of the reconstructed
      // function at the nodes of the cell

      int dim = this->mesh_.space_dimension();
      std::vector<Point<D>> cellcoords;
      this->mesh_.cell_get_coordinates(cellid, &cellcoords);
      
      for (auto coord : cellcoords) {
	Vector<D> vec = coord-ls_coords[0];
	double diff = dot(grad, vec);
	double extremeval = (diff > 0.0) ? maxval : minval;
	double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff;
	phi = std::min(phi_new, phi);
      }
    }
    
    // Limited gradient is phi*grad
    return phi*grad;
  }
}  // namespace Portage

#endif  // SRC_INTERPOLATE_GRADIENT_H_
