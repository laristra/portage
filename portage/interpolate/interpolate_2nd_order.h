/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_
#define PORTAGE_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_

#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <iostream>
#include <utility>
#include <vector>

#include "portage/interpolate/gradient.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/support/portage.h"


namespace Portage {

/*!
  @class Interpolate_2ndOrder interpolate_2nd_order.h
  @brief Interpolate_2ndOrder does a 2nd order interpolation of scalars
  @tparam MeshType The type of the mesh wrapper used to access mesh info
  @tparam StateType The type of the state manager used to access data.
  @tparam OnWhatType The type of entity-based data we wish to interpolate;
  e.g. does it live on nodes, cells, edges, etc.

  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

  [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
  (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
  Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
  pp. 305-321, 1987.

  @todo Template on variable type (YES)
*/


template<int D, 
         Entity_kind on_what, 
         typename SourceMeshType, 
         typename TargetMeshType, 
         typename StateType,
         template<class, int> class InterfaceReconstructorType =
          DummyInterfaceReconstructor>
class Interpolate_2ndOrder {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType>;
#endif


 public:

/*!
  @brief Constructor with interface reconstructor
  @param[in] source_mesh The mesh wrapper used to query source mesh info
  @param[in] target_mesh The mesh wrapper used to query target mesh info
  @param[in] source_state The state-manager wrapper used to query field info
  @param[in] ir The interface_reconstructor used to query matpoly's on source mesh
*/
#ifdef HAVE_TANGRAM
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state, 
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir),
      interp_var_name_("VariableNameNotSet"),
      limiter_type_(NOLIMITER),
      source_vals_(nullptr) {}
#endif

  /*!
    @brief Constructor without interface reconstructor 
    @param[in] source_mesh The mesh wrapper used to query source mesh info
    @param[in] target_mesh The mesh wrapper used to query target mesh info
    @param[in] source_state The state-manager wrapper used to query field info
  */

  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      limiter_type_(NOLIMITER),
      source_vals_(nullptr) {}

  /// Copy constructor (disabled)
  //  Interpolate_2ndOrder(const Interpolate_2ndOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_2ndOrder & operator = (const Interpolate_2ndOrder &) = delete;

  /// Destructor
  ~Interpolate_2ndOrder() {}

  /// Set the material we are operating on. 

  int set_material(int m) {
    matid_ = m;
  } //set_material

  /// Set the name of the interpolation variable and the limiter type

  void set_interpolation_variable(std::string const & interp_var_name,
                                  LimiterType limiter_type = NOLIMITER) {
    interp_var_name_ = interp_var_name;
    limiter_type_ = limiter_type;

    // Extract the field data from the statemanager
    field_type_ = source_state_.field_type(on_what, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(on_what, interp_var_name,
                                  &source_vals_);
    else
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
    
    
    // Compute the limited gradients for the field 
  #ifdef HAVE_TANGRAM  
    Limited_Gradient<D, on_what, SourceMeshType, StateType, InterfaceReconstructorType>
        limgrad(source_mesh_, source_state_, interp_var_name_, limiter_type_, 
                interface_reconstructor_);
   limgrad.set_material(matid_, interp_var_name_, limiter_type_);      
  #else      
    Limited_Gradient<D, on_what, SourceMeshType, StateType>
        limgrad(source_mesh_, source_state_);
  #endif 

   // Get the correct number of source cells for which the gradient has to be computed
   #ifdef HAVE_TANGRAM
    std::vector<int> cellids;
    source_state_.mat_get_cells(matid_, &cellids);
    int nentities =  cellids.size();
   #else
    int nentities = source_mesh_.end(on_what)-source_mesh_.begin(on_what);
   #endif   

    gradients_.resize(nentities);
    
    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform
  #ifdef HAVE_TANGRAM
    Portage::transform(cellids.begin(), cellids.end(),
                       gradients_.begin(), limgrad);
  #else 
    Portage::transform(source_mesh_.begin(on_what), source_mesh_.end(on_what),
                       gradients_.begin(), limgrad);
  #endif  
  } //set_interpolation_variable 


  /*!
    @brief Functor to do the actual interpolate calculation
    @param[in] cells_and_weights A pair of two vectors
    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order interpolate, only the first element (or zero'th moment)
    of the weights vector (i.e. the volume of intersection) is used. Source
    entities may be repeated in the list if the intersection of a target entity
    and a source entity consists of two or more disjoint pieces
    @param[in] targetCellID The index of the target cell.

    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.
    @todo must remove assumption that field is scalar
  */

  double operator() (int const targetCellID,
                     std::vector<Weights_t> const & sources_and_weights) const {
    // not implemented for all types - see specialization for cells and nodes

    std::cerr << "Interpolation operator not implemented for this entity type"
              << std::endl;
  }

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string interp_var_name_;
  LimiterType limiter_type_;
  double const * source_vals_;

  // Portage::vector is generalization of std::vector and
  // Portage::Vector<D> is a geometric vector
  Portage::vector<Vector<D>> gradients_;

  int matid_;
  Field_type field_type_;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
}; // class Interpolate_2ndOrder




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief 2nd order interpolate class specialization for cells
  @param[in] cells_and_weights Pair containing vector of contributing source
  cells and vector of contribution weights
*/

template<int D,
         typename SourceMeshType, 
         typename TargetMeshType, 
         typename StateType,
         template<class, int> class InterfaceReconstructorType>
class Interpolate_2ndOrder<D, 
                           CELL, 
                           SourceMeshType, 
                           TargetMeshType, 
                           StateType, 
                           InterfaceReconstructorType> {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType>;
#endif

 public:

// Constructor with interface reconstructor
#ifdef HAVE_TANGRAM
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state, 
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir),
      interp_var_name_("VariableNameNotSet"),
      limiter_type_(NOLIMITER),
      source_vals_(nullptr) {}
#endif

//Constructor without interface reconstructor
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      limiter_type_(NOLIMITER),
      source_vals_(nullptr) {}

  
  /// Set the name of the interpolation variable and the limiter type

  void set_interpolation_variable(std::string const & interp_var_name,
                                  LimiterType limiter_type = NOLIMITER) {

    interp_var_name_ = interp_var_name;
    limiter_type_ = limiter_type;

    // Extract the field data from the statemanager
    field_type_ = source_state_.field_type(CELL, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(CELL, interp_var_name, &source_vals_);
    else
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
    
    // Compute the limited gradients for the field 
  #ifdef HAVE_TANGRAM  
    Limited_Gradient<D, CELL, SourceMeshType, StateType, InterfaceReconstructorType>
        limgrad(source_mesh_, source_state_, interp_var_name_, limiter_type_, 
                interface_reconstructor_);
     limgrad.set_material(matid_, interp_var_name_, limiter_type_);      
  #else      
    Limited_Gradient<D, CELL, SourceMeshType, StateType>
        limgrad(source_mesh_, source_state_, interp_var_name_, limiter_type_);
  #endif 

   // Get the correct number of source cells for which the gradient has to be computed
   #ifdef HAVE_TANGRAM
    std::vector<int> cellids;
    source_state_.mat_get_cells(matid_, &cellids);
    int nentities =  source_state_.mat_get_num_cells(matid_);
   #else
    int nentities = source_mesh_.num_entities(CELL);
   #endif  

    gradients_.resize(nentities);
    
    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform
  #ifdef HAVE_TANGRAM
    Portage::transform(cellids.begin(), cellids.end(),
                       gradients_.begin(), limgrad);
  #else 
    Portage::transform(source_mesh_.begin(CELL), source_mesh_.end(CELL),
                       gradients_.begin(), limgrad);
  #endif  
  } //set_interpolation_variable


  /// Copy constructor (disabled)
  //  Interpolate_2ndOrder(const Interpolate_2ndOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_2ndOrder & operator = (const Interpolate_2ndOrder &) = delete;

  /// Destructor
  ~Interpolate_2ndOrder() {}


  /// Set the material we are operating on (MM interpolate not implemented yet)

  int set_material(int m) {
    matid_ = m ;
  } //set_material

  /*!
    @brief   Functor to do the 2nd order interpolation of cell values

    @param[in] targetCellID The index of the target cell.

    @param[in] sources_and_weights Vector of source mesh entities and
    corresponding weight vectors.  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order interpolation, only the first element (or zero'th moment) of
    the weights vector (i.e. the volume of intersection) is
    used. Source entities may be repeated in the list if the
    intersection of a target entity and a source entity consists of
    two or more disjoint pieces

    @todo must remove assumption that field is scalar
  */

  double operator() (int const targetCellID,
                     std::vector<Weights_t> const & sources_and_weights) const
  {
     int nsrccells = sources_and_weights.size();
     if (!nsrccells) {
#ifdef DEBUG
     std::cerr << "WARNING: No source cells contribute to target cell." <<
        std::endl;
#endif
      return 0.0;
    }

    double totalval = 0.0;
    double wtsum0 = 0.0;

    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)

    /// @todo Should use zip_iterator here but I am not sure I know how to

    double vol = target_mesh_.cell_volume(targetCellID);
    int nsummed = 0;

    //Loop over source cells
    for (int j = 0; j < nsrccells; ++j) {

      //Get source cell and the intersection weights 
      int srccell = sources_and_weights[j].entityID;
      std::vector<double> xsect_weights = sources_and_weights[j].weights;
      double xsect_volume = xsect_weights[0];

      double eps = 1e-12;
      if (xsect_volume/vol <= eps) continue;  // no intersection

      // Obtain source cell centroid
      Point<D> srccell_centroid;

#ifdef HAVE_TANGRAM
      int nmats = source_state_.cell_get_num_mats(srccell);
      std::vector<int> cellmats;
      source_state_.cell_get_mats(srccell, &cellmats);
      
      if (!nmats || (nmats == 1 && cellmats[0] == matid_)) 
      { //pure cell
        source_mesh_.cell_centroid(srccell, &srccell_centroid);
      } 
     else 
      { //multi-material cell
         assert(interface_reconstructor_ != nullptr);  // cannot be nullptr
        
         if (std::find(cellmats.begin(), cellmats.end(), matid_) !=
              cellmats.end()) 
         { // mixed cell containing this material 

          //Obtain matpoly's for this material
          Tangram::CellMatPoly<D> const& cellmatpoly =
              interface_reconstructor_->cell_matpoly_data(srccell);
          std::vector<Tangram::MatPoly<D>> matpolys =
              cellmatpoly.get_matpolys(matid_);

          int cnt = 0;
          for (int k = 0; k < D; k++) srccell_centroid[k]=0;

          //Compute centroid of all matpoly's 
          for (int j = 0; j < matpolys.size(); j++)
          {
            std::vector<Tangram::Point<D>> tpnts = matpolys[j].points();
            cnt += tpnts.size();
            for (int i = 0; i < tpnts.size(); i++)
               for (int k = 0; k < D; k++)
                 srccell_centroid[k] += tpnts[i][k];
          }
          srccell_centroid = srccell_centroid/cnt; 
       }
     }
#else
     source_mesh_.cell_centroid(srccell, &srccell_centroid);
#endif
    
   //Computer intersection centroid
    Point<D> xsect_centroid;
    for (int i = 0; i < D; ++i)
      xsect_centroid[i] = xsect_weights[1+i]/xsect_volume;  // (1st moment)/vol

   //Get correct source cell index
   int srcindex;
#ifdef HAVE_TANGRAM
   srcindex = source_state_.cell_index_in_material(srccell, matid_);
#else
   srcindex = srccell;
#endif  

    Vector<D> gradient = gradients_[srcindex];
    Vector<D> vec = xsect_centroid - srccell_centroid;
    double val = source_vals_[srcindex] + dot(gradient, vec);
    val *= xsect_volume;
    totalval += val;
    wtsum0 += xsect_volume;
    nsummed++;
  }

  // Normalize the value by the volume of the intersection of the target cells
  // with the source mesh. This will do the right thing for single-material
  // and multi-material remap (conservative and constant preserving) if there
  // is NO mismatch between source and target mesh boundaries. IF THERE IS A
  // MISMATCH, THIS WILL PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE. THEN
  // WE HAVE TO DO A SEMI-LOCAL OR GLOBAL REPAIR.

  if (nsummed)
    totalval /= wtsum0;
  else
    totalval = 0.0;

#ifdef DEBUG
  static bool first = true;
  if (first && fabs((vol-wtsum0)/vol) > 1.0e-10) {
    std::cerr <<
        "WARNING: Meshes may be mismatched in the neighborhood of cell " <<
        targetCellID << " in the target mesh (and maybe other places too)\n";
    first = false;
  }
#endif

    return totalval;

} //operator()

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string interp_var_name_;
  LimiterType limiter_type_;
  double const * source_vals_;

  // Portage::vector is generalization of std::vector and
  // Portage::Vector<D> is a geometric vector
  Portage::vector<Vector<D>> gradients_;

  int matid_;
  Field_type field_type_;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
}; //2nd order interpolate class specialization for cells


//////////////////////////////////////////////////////////////////////////////
/*!
  @brief 2nd order interpolate class specialization for nodes
*/
template<int D,
         typename SourceMeshType, 
         typename TargetMeshType, 
         typename StateType,
         template<class, int> class InterfaceReconstructorType>
class Interpolate_2ndOrder<D, 
                           NODE, 
                           SourceMeshType, 
                           TargetMeshType, 
                           StateType,
                           InterfaceReconstructorType> {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType>;
#endif

 public:

// Constructor with interface reconstructor
#ifdef HAVE_TANGRAM
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state, 
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir),
      interp_var_name_("VariableNameNotSet"),
      limiter_type_(NOLIMITER),
      source_vals_(nullptr) {}
#endif

//Constructor without interface reconstructor
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      limiter_type_(NOLIMITER),
      source_vals_(nullptr) {}

  /// Copy constructor (disabled)
  //  Interpolate_2ndOrder(const Interpolate_2ndOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_2ndOrder & operator = (const Interpolate_2ndOrder &) = delete;

  /// Destructor
  ~Interpolate_2ndOrder() {}


  /// Set the name of the interpolation variable and the limiter type

  void set_interpolation_variable(std::string const & interp_var_name,
                                  LimiterType limiter_type = NOLIMITER) {

    interp_var_name_ = interp_var_name;
    limiter_type_ = limiter_type;

    // Extract the field data from the statemanager
    field_type_ = source_state_.field_type(NODE, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(NODE, interp_var_name,
                                  &source_vals_);
    else {
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
      std::cerr << "Cannot remap NODE-centered multi-material data" << "\n";
    }
    

    // Compute the limited gradients for the field
    
// #ifdef HAVE_TANGRAM  
//    Limited_Gradient<D, NODE, SourceMeshType, StateType, InterfaceReconstructorType>
//        limgrad(source_mesh_, source_state_, interp_var_name_, limiter_type_, 
//          interface_reconstructor_);
//  #else 
   Limited_Gradient<D, NODE, SourceMeshType, StateType>
        limgrad(source_mesh_, source_state_, interp_var_name_, limiter_type_);
//  #endif
  
    int nentities = source_mesh_.end(NODE)-source_mesh_.begin(NODE);
    gradients_.resize(nentities);
    
    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)
    
    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform

    Portage::transform(source_mesh_.begin(NODE), source_mesh_.end(NODE),
                       gradients_.begin(), limgrad);
  }


  /// Set the material we are operating on. 

  int set_material(int m) {
    matid_ = m;
  } //set_material

  /*!
    @brief Functor to do the 2nd order interpolation of node values
    @param[in] sources_and_weights      A pair of two vectors

    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.

    @c sources_and_weights.second() is the
    @param[in] targetCellID The index of the target cell.

    @param[in] sources_and_weights Vector of source mesh entities and
    corresponding weight vectors.  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order interpolation, only the first element (or zero'th moment) of
    the weights vector (i.e. the volume of intersection) is
    used. Source entities may be repeated in the list if the
    intersection of a target entity and a source entity consists of
    two or more disjoint pieces

    @todo must remove assumption that field is scalar
  */

  double operator() (const int targetNodeID,
                     std::vector<Weights_t> const & sources_and_weights) const
  {
     int nsrcnodes = sources_and_weights.size();
     if (!nsrcnodes) {
#ifdef DEBUG
      std::cerr << "WARNING: No source nodes contribute to target node." <<
          std::endl;
#endif
      return 0.0;
    }

    double totalval = 0.0;
    double wtsum0 = 0.0;

    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)

    /// @todo Should use zip_iterator here but I am not sure I know how to

    double vol = target_mesh_.dual_cell_volume(targetNodeID);
    int nsummed = 0;
    for (int j = 0; j < nsrcnodes; ++j) {
      int srcnode = sources_and_weights[j].entityID;
      std::vector<double> xsect_weights = sources_and_weights[j].weights;
      double xsect_volume = xsect_weights[0];

      double eps = 1e-12;
      if (xsect_volume/vol <= eps) continue;  // no intersection

      // note: here we are getting the node coord, not the centroid of
      // the dual cell
      Point<D> srcnode_coord;
      source_mesh_.node_get_coordinates(srcnode, &srcnode_coord);

      Point<D> xsect_centroid;
      for (int i = 0; i < D; ++i)
        // (1st moment)/(vol)
        xsect_centroid[i] = xsect_weights[1+i]/xsect_volume;

      Vector<D> gradient = gradients_[srcnode];
      Vector<D> vec = xsect_centroid - srcnode_coord;
      double val = source_vals_[srcnode] + dot(gradient, vec);
      val *= xsect_volume;
      totalval += val;
      wtsum0 += xsect_volume;
      nsummed++;
    }

    // Normalize the value by volume of the target dual cell

    // Normalize the value by the volume of the intersection of the target cells
    // with the source mesh. This will do the right thing for single-material
    // and multi-material remap (conservative and constant preserving) if there
    // is NO mismatch between source and target mesh boundaries. IF THERE IS A
    // MISMATCH, THIS WILL PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE. THEN
    // WE HAVE TO DO A SEMI-LOCAL OR GLOBAL REPAIR.
  
    if (nsummed)
      totalval /= wtsum0;
    else
      totalval = 0.0;

#ifdef DEBUG
    static bool first = true;
    if (first && fabs((vol-wtsum0)/vol) > 1.0e-10) {
      std::cerr <<
          "WARNING: Meshes may be mismatched in the neighborhood of node " <<
          targetNodeID << " in the target mesh (and maybe other places too) \n";
      first = false;
    }
#endif


    return totalval;
} //operator()

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string interp_var_name_;
  LimiterType limiter_type_;
  double const * source_vals_;

  // Portage::vector is generalization of std::vector and
  // Portage::Vector<D> is a geometric vector
  Portage::vector<Vector<D>> gradients_;

  int matid_;
  Field_type field_type_;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
};


}  // namespace Portage

#endif  // PORTAGE_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_
