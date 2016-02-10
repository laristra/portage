/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

#include<algorithm>
#include<vector>
#include<iterator>
#include<sys/time.h>

#include "portage/support/portage.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/search/search_kdtree2.h"
#include "portage/search/search_kdtree3.h"
#include "portage/intersect/intersectClipper.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/remap/remap_1st_order.h"
#include "portage/remap/remap_2nd_order.h"

/*!
  @file driver.h
  @brief Example driver for mapping between two Jali meshes.

  This should serve as a good example for how to write your own driver routine
  and datastructures.
 */

namespace Portage {

/*! 
    @class MeshWrapperDual "driver.h"
    @brief Wrapper for dual mesh.

    Utilizes a Jali_Mesh_Wrapper to the original mesh, but treats
    the nodes of the original mesh as the centroids of the dual mesh.
 */
class MeshWrapperDual { // cellid is the dual cell (i.e. node) id
 public:

  /*!
    @brief Constructor of a wrapper to a 2d mesh.
    @param[in] w Jali_Mesh_Wrapper to original mesh.
  */
  MeshWrapperDual(const Jali_Mesh_Wrapper &w) : w_(w) {}

  /*!
    @brief Get the spatial dimensions of the mesh.
    @return The spatial dimension of the mesh.
  */
  int space_dimension() const { return w_.space_dimension(); }


  /*!
    @brief Get the number of cells on this processor.
    @return The number of cells owned on this processor.
  */
  int num_owned_cells() const { return w_.num_owned_nodes(); }

  /*!
    @brief Get the number of ghost (at domain boundaries @e and processor boundaries)
    cells for this processor.
    @return The number of ghost cells on this processor.
  */
  int num_ghost_cells() const { return w_.num_ghost_nodes(); }

  /*!
    @brief Gets the coordinates of the cell centroid in the dual mesh.
    @param[in] dualcellid The dual cell id (i.e. the node id in the original mesh).
    @param[in,out] xylist The list of (x,y) coordinate pairs for each node of the
    dual cell given by @c dualcellid
   */
  void cell_get_coordinates(int const dualcellid,
                            std::vector<std::pair<double,double> > *xylist) const {
    w_.dual_cell_get_coordinates(dualcellid, xylist);
  }

  /*!
    @brief Gets the coordinates of the cell centroid in the dual mesh.
    @param[in] dualcellid The dual cell id (i.e. the node id in the original
    mesh).
    @param[in,out] xyzlist The list of (x,y,z) coordinate tuples for each node
    of the dual cell given by @c dualcellid.
   */
  void cell_get_coordinates(int const dualcellid,
			    std::vector<std::tuple<double,double,double> > *xyzlist) const {
    w_.dual_cell_get_coordinates(dualcellid, xyzlist);
  }

  /*!
    @brief Get an iterator to the start of the vector of @e entity-type objects
    in the dual mesh.
    @param[in] entity Which type of data do you want to iterate over (e.g.
    @c CELL, @c NODE, etc.)
    @return An iterator pointing to the beginning of the list of @c entity data
    for this mesh.
   */
  counting_iterator begin(Entity_kind const entity) const {
    if (entity == NODE) return w_.begin(CELL);
    return w_.begin(NODE);
  }
  
  /*!
    @brief Get an iterator to the end of the vector of @e entity-type objects in
    the dual mesh.
    @param[in] entity Which type of data do you want to iterate over (e.g.
    @c CELL, @c NODE, etc.)
    @return An iterator pointing to the end of the list of @c entity data for
    this mesh.
   */
  counting_iterator end(Entity_kind const entity) const {
    if (entity == NODE) return w_.end(CELL);
    return w_.end(NODE);
  }

  /*!
    @brief Gets the coordinates of the cell centroid in the dual mesh.
    @param[in] dualcellid The dual cell id (i.e. the node id in the original
    mesh).
    @return The list of (x,y) coordinate pairs for each node of the
    dual cell given by @c dualcellid
    @todo Remove this in favor of @c cell_get_coordinates() ?
   */
  std::vector<std::pair<double, double> > cellToXY(int dualcellID) const{
    std::vector<std::pair<double, double> > cellPoints;
    cell_get_coordinates(dualcellID, &cellPoints);
    return cellPoints;
  }

  /*!
    @brief Get the IDs of all cells that share a node with the specified cell
    <em> of the original mesh </em>.

    Sharing of cells is determined from the Parallel_type (e.g. @c OWNED,
    @c GHOST, @c ALL ).
    
    @param[in] dualcellID The cell ID for which you would like to find the
    neighbors.
    @param[in] ptype The type of data you want (@c OWNED, @c GHOST, @c ALL)
    @param[in,out] adjcells List of IDs of adjacent cells.
   */
  void cell_get_node_adj_cells(int const dualcellID,
                               Parallel_type const ptype,
                               std::vector<int> *adjcells) const {
    w_.node_get_cell_adj_nodes(dualcellID,ptype,adjcells);
  }
  
  /*!
    @brief Get the IDs of all cells that share a node with the specified cell of
    the <em> dual mesh </em>.

    Sharing of cells is determined from the Parallel_type (e.g. @c OWNED,
    @c GHOST, @c ALL).
    
    @param[in] dualnodeID The cell ID for which you would like to find the
    neighbors.
    @param[in] ptype The type of data you want (@c OWNED, @c GHOST, @c ALL)
    @param[in,out] adjnodes List of IDs of adjacent cells.

    @todo Clarify this wrt to @c MeshWrapperDual::cell_get_node_adj_cells()
   */
  void dual_cell_get_node_adj_cells(int const dualnodeID,
                                    Parallel_type const ptype,
                                    std::vector<int> *adjnodes) const {
    w_.cell_get_node_adj_cells(dualnodeID,ptype,adjnodes);
  }

  /*!
    @brief Get the coordinates of the centroid of a given dual mesh cell ID.
    @param[in] dualcellID ID of the cell in the dual mesh.
    @param[in,out] centroid (x,y,z) coordinates of the cell center (for 3d).
   */
  void cell_centroid(int const dualcellID, std::vector<double> *centroid) const {
    w_.dual_cell_centroid(dualcellID, centroid); 
  }

  /*!
    @brief Get the coordinates of the centroid of a given dual mesh cell ID.
    @param[in] dualcellID ID of the cell in the dual mesh.
    @param[in,out] centroid (x,y,z) coordinates of the cell center (for 3d).
    @todo Clarify this wrt to @c MeshWrapperDual::cell_centroid().
   */
  void dual_cell_centroid(int const dualnodeID, std::vector<double> *centroid) const {
    w_.cell_centroid(dualnodeID, centroid);
  }
  /*!
    @brief Get the coordinates of the points that make up the wedge.

    A wedge corresponds to a cell center, a face center, and two points that
    share an edge.

    @param[in] dualcellid ID of the cell in the dual mesh.
    @param[in,out] wcoords (x,y,z) coordinates of each of the four points that
    comprise the tetrahedron that is the wedge.
   */
  void wedges_get_coordinates(int const dualcellid,
			      std::vector<std::array<std::array<double, 3>, 4> > *wcoords) const {
    w_.dual_wedges_get_coordinates(dualcellid, wcoords);
  }

private:
    const Jali_Mesh_Wrapper &w_;
};


// Forward definition
template <typename SearchType, typename IsectType, typename RemapType>
struct composerFunctor;


/*!
    @class Driver "driver.h"
    @brief Driver provides the API to mapping from one mesh to another.
    @tparam Mesh_Wrapper A lightweight wrapper to a specific mesh implementation
    that provides certain functionality.  See Jali_Mesh_Wrapper for an example.
 */
template <class Mesh_Wrapper>
class Driver
{
 public:
  
  /*!
    @brief Constructor for running the remap driver.
    @param[in] sourceMesh A @c Mesh_Wrapper to the source mesh.
    @param[in] sourceState A state manager for the data that lives on the source mesh.
    @param[in] targetMesh A @c Mesh_Wrapper to the target mesh.
    @param[in,out] targetState A state manager for the data that will be mapped to
    the target mesh.
   */
  Driver(Entity_kind remapEntity, 
         Mesh_Wrapper const & sourceMesh, 
         Jali_State_Wrapper const & sourceState,           
         Mesh_Wrapper const & targetMesh,
         Jali_State_Wrapper & targetState) 
      : remap_entity_(remapEntity), 
        source_mesh_(sourceMesh), source_state_(sourceState),
        target_mesh_(targetMesh), target_state_(targetState),
        remap_order_(1), dim_(sourceMesh.space_dimension())
  {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
  }
  
  /// Copy constructor (disabled)
  Driver(const Driver &) = delete;
  
  /// Assignment operator (disabled)
  Driver & operator = (const Driver &) = delete;
  
  /// Destructor
  ~Driver() {}
  
  /*!
    @brief Specify the names of the variables to be remapped
    @param[in] remap_var_names_in A list of variable names of the variables to remap
    from the source mesh to the target mesh.
   */
  void set_remap_var_names(std::vector<std::string> &remap_var_names_in) {
    remap_var_names_ = remap_var_names_in;
  }
  

  /*!
    @brief Get the names of the variables to be remapped.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> remap_var_names() {
    return remap_var_names_;
  }
  
  /// Set the order of accuracy of remap
  void set_remap_order(unsigned int const order) {
    remap_order_ = order;
  }
  
  /*!
    @brief Get the order of accuracy of remap
    @return The order of accuracy for the remap.
  */
  unsigned int remap_order() {
    return remap_order_;
  }

  /*!
    @brief Get the dimensionality of the meshes.
    @return The dimensionality of the meshes.
   */
  unsigned int dim() {
    return dim_;
  }
  
  /*!
    @brief This method calls the search, intersect, and remap routines
    needed to map one mesh to another.

    Most of the heavy lifting in this routine is via a Portage::transform() over
    the cells in the target mesh, applying a custom composerFunctor() that specifies
    how the search, intersect, and remap calculations should be performed.
  */
  void run()
  {
    std::printf("in Driver::run()...\n");

    // FIXME: most of this is duplicated for 3d...
    if (dim() == 2) {
      // 2d
    
      // Get an instance of the desired search algorithm type
      const SearchKDTree2<Mesh_Wrapper,Mesh_Wrapper>
        search(source_mesh_, target_mesh_);
    
      // Get an instance of the desired intersect algorithm type
      const IntersectClipper<Mesh_Wrapper, Mesh_Wrapper>
        intersect{source_mesh_, target_mesh_};
    
      int numTargetCells = target_mesh_.num_owned_cells();
      std::cout << "Number of target cells in target mesh "
		<< numTargetCells << std::endl;
    
      // Ask for a StateVector with the name remap_var_names_[0] to be
      // added to the targetState_. If it is already present, the
      // existing StateVector reference is returned. If its not present,
      // it is added. This logic needs to be reversed. The find function
      // should add it if it is not found (if so requested).
    
      double *target_field_raw = NULL;
      target_state_.get_data(remap_entity_,remap_var_names_[0],&target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      // Check order of accuracy of remap and do the appropriate thing

      if (remap_order() == 1) {

	std::cout << "Remapping variable " << remap_var_names_[0]
		  << " using a 1st order accurate algorithm" << std::endl;

	// Eventually put this in a loop over remap variable names as well
	// Assume for now that we are only doing cell-based remap
	const Remap_1stOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> 
          remap(source_mesh_, source_state_, remap_entity_, remap_var_names_[0]);
      
	composerFunctor<SearchKDTree2<Mesh_Wrapper,Mesh_Wrapper>,
	                IntersectClipper<Mesh_Wrapper, Mesh_Wrapper>,
	                Remap_1stOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> >
          composer(&search, &intersect, &remap, remap_var_names_[0]);

	// This populates targetField with the doubles returned from
	// the final remap For some reason (perhaps some code outside
	// portage is pulling in a "using thrust" command), the
	// compiler is not able to disambiguate Portage::transform and
	// thrust::transform here. So, be explicit that we want
	// Portage::transform

#ifdef ENABLE_PROFILE
        __itt_resume();
#endif

	struct timeval begin, end, diff;
	gettimeofday(&begin, 0);

	Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
			   (counting_iterator)(target_mesh_.end(CELL)),
			   target_field,composer);
      
#ifdef ENABLE_PROFILE
        __itt_pause();
#endif

	gettimeofday(&end, 0);
	timersub(&end, &begin, &diff);
	float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
	std::cout << "Transform Time: " << seconds << std::endl;
	
      } // done first order remap
      else {
      
	if (remap_order() != 2) {
	  std::cerr << "Remap order can be 1 or 2 only. "
		    << "Doing 2nd order remap" << std::endl;
	}
	
	std::cout << "Remapping variable " << remap_var_names_[0]
		  << " using a 2nd order accurate algorithm" << std::endl;
      
	/// @todo Eventually put this in a loop over remap variable names as well
	// Assume for now that we are only doing cell-based remap
	const Remap_2ndOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> 
	  remap(source_mesh_, source_state_, remap_entity_, remap_var_names_[0],
                NOLIMITER);
          
	composerFunctor<SearchKDTree2<Mesh_Wrapper,Mesh_Wrapper>,
	  IntersectClipper<Mesh_Wrapper, Mesh_Wrapper>,
	  Remap_2ndOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> >
          composer(&search, &intersect, &remap, remap_var_names_[0]);
      
	// This populates targetField with the doubles returned from
	// the final remap For some reason (perhaps some code outside
	// portage is pulling in a "using thrust" command), the
	// compiler is not able to disambiguate Portage::transform and
	// thrust::transform here. So, be explicit that we want
	// Portage::transform
      
#ifdef ENABLE_PROFILE
        __itt_resume();
#endif

	struct timeval begin, end, diff;
	gettimeofday(&begin, 0);

	Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
			   (counting_iterator)(target_mesh_.end(CELL)),
			   target_field,composer);

#ifdef ENABLE_PROFILE
        __itt_pause();
#endif

	gettimeofday(&end, 0);
	timersub(&end, &begin, &diff);
	float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
	std::cout << "Transform Time: " << seconds << std::endl;
      } // done remap order test 

    }   // done 2d test
    else {
      // 3d
      
      // Get an instance of the desired search algorithm type
      const SearchKDTree3<Mesh_Wrapper,Mesh_Wrapper>
        search(source_mesh_, target_mesh_);
    
      // Get an instance of the desired intersect algorithm type
      const IntersectR3D<Mesh_Wrapper, Mesh_Wrapper>
        intersect{source_mesh_, target_mesh_};
    
      int numTargetCells = target_mesh_.num_owned_cells();
      std::cout << "Number of target cells in target mesh "
		<< numTargetCells << std::endl;
    
      // Ask for a StateVector with the name remap_var_names_[0] to be
      // added to the targetState_. If it is already present, the
      // existing StateVector reference is returned. If its not present,
      // it is added. This logic needs to be reversed. The find function
      // should add it if it is not found (if so requested).
    
      double *target_field_raw = NULL;
      target_state_.get_data(remap_entity_,remap_var_names_[0],&target_field_raw);
      Portage::pointer<double> target_field(target_field_raw);

      // Check order of accuracy of remap and do the appropriate thing

      if (remap_order() == 1) {

	std::cout << "Remapping variable " << remap_var_names_[0]
		  << " using a 1st order accurate algorithm" << std::endl;

	// Eventually put this in a loop over remap variable names as well
	// Assume for now that we are only doing cell-based remap
	const Remap_1stOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> 
          remap(source_mesh_, source_state_, remap_entity_, remap_var_names_[0]);
      
          
	composerFunctor<SearchKDTree3<Mesh_Wrapper,Mesh_Wrapper>,
	                IntersectR3D<Mesh_Wrapper, Mesh_Wrapper>,
	                Remap_1stOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> >
          composer(&search, &intersect, &remap, remap_var_names_[0]);

	// This populates targetField with the doubles returned from
	// the final remap For some reason (perhaps some code outside
	// portage is pulling in a "using thrust" command), the
	// compiler is not able to disambiguate Portage::transform and
	// thrust::transform here. So, be explicit that we want
	// Portage::transform

#ifdef ENABLE_PROFILE
	__itt_resume();
#endif

	struct timeval begin, end, diff;
	gettimeofday(&begin, 0);

	Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
			   (counting_iterator)(target_mesh_.end(CELL)),
			   target_field,composer);

#ifdef ENABLE_PROFILE
	__itt_pause();
#endif

	gettimeofday(&end, 0);
	timersub(&end, &begin, &diff);
	float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
	std::cout << "Transform Time: " << seconds << std::endl;
      } // done first order remap
      else {
      
      	if (remap_order() != 2)
      	  std::cerr << "Remap order can be 1 or 2 only. "
		    << "Doing 2nd order remap" <<  std::endl;

      	std::cout << "Remapping variable " << remap_var_names_[0]
		  << " using a 2nd order accurate algorithm" << std::endl;
      
      	// Eventually put this in a loop over remap variable names as well
      	// Assume for now that we are only doing cell-based remap
      	const Remap_2ndOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind>
          remap(source_mesh_, source_state_, remap_entity_, remap_var_names_[0],
                NOLIMITER);
      
          
      	composerFunctor<SearchKDTree3<Mesh_Wrapper,Mesh_Wrapper>,
      	                IntersectR3D<Mesh_Wrapper, Mesh_Wrapper>,
      	                Remap_2ndOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> >
          composer(&search, &intersect, &remap, remap_var_names_[0]);
      
      	// This populates targetField with the doubles returned from
      	// the final remap For some reason (perhaps some code outside
      	// portage is pulling in a "using thrust" command), the
      	// compiler is not able to disambiguate Portage::transform and
      	// thrust::transform here. So, be explicit that we want
      	// Portage::transform

#ifdef ENABLE_PROFILE
	__itt_resume();
#endif

	struct timeval begin, end, diff;
	gettimeofday(&begin, 0);

	Portage::transform((counting_iterator)(target_mesh_.begin(CELL)),
			   (counting_iterator)(target_mesh_.end(CELL)),
			   target_field,composer);

#ifdef ENABLE_PROFILE
	__itt_pause();
#endif

        gettimeofday(&end, 0);
	timersub(&end, &begin, &diff);
	float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
	std::cout << "Transform Time: " << seconds << std::endl;
      } // done remap test
    } // done 3d test
  } // run()

private:
  
    Mesh_Wrapper  const & source_mesh_;
    Mesh_Wrapper const & target_mesh_;
    Jali_State_Wrapper const & source_state_;
    Jali_State_Wrapper & target_state_;
    std::vector<std::string> remap_var_names_;
    Entity_kind remap_entity_;
    unsigned int remap_order_;
    unsigned int dim_;

}; // class Driver


/*!
  @struct composerFunctor "driver.h"
  @brief This functor is used inside a Portage::transform() inside Driver::run() to
  actually do the search, intersect, and remap calculations.
  @tparam SearchType The type of search method (e.g. SearchSimple or SearchKDTree3).
  @tparam IsectType The type of intersect method (e.g. IntersectClipper).
  @tparam RemapType The type of remap method (e.g. Remap_1stOrder or Remap_2ndOrder).
 */
template <typename SearchType, typename IsectType, typename RemapType>
struct composerFunctor
{
    const SearchType* search_;         ///< search method (e.g. SearchSimple)
    const IsectType* intersect_;       ///< intersect method (e.g. IntersectClipper)
    const RemapType* remap_;           ///< remap method (e.g. Remap_2ndOrder)
    const std::string remap_var_name_; ///< variable name to remap
    //----------------------------------------

    /*!
      @brief Constructor.
      @param[in] searcher The search method to use (e.g. SearchSimple)
      @param[in] intersecter The intersect method to use (e.g. IntersectClipper)
      @param[in] remapper The remap method to use (e.g. Remap_2ndOrder)
      @param[in] remap_var_name The name of the variable to remap
     */
    composerFunctor(const SearchType* searcher, const IsectType* intersecter, 
                    const RemapType* remapper, const std::string remap_var_name)
			: search_(searcher), intersect_(intersecter), remap_(remapper), 
              remap_var_name_(remap_var_name) { }

     /*!
       @brief Operator for making this struct a functor

       This is called from within a Portage::transform() operation that iterates
       over the cells in a target mesh.  

       @param[in] targetCellindex The cell ID in the target mesh that this functor
       is currently operating on.

       @return Value of the field @c remap_var_name in the target mesh cell with ID
       @c targetCellIndex.
      */
    double operator()(int const targetCellIndex)
    {

        // Search for candidates and return their cells indices
        std::vector<int> candidates;
        search_->search(targetCellIndex, &candidates);

        // Intersect target cell with cells of source mesh and return the
        // moments of intersection
        std::vector<std::vector<std::vector<double> > > 
                moments(candidates.size());
        for (int i=0;i<candidates.size();i++)
            moments[i] = (*intersect_)(candidates[i], targetCellIndex);

        // Compute new value on target cell based on source mesh
        // values and intersection moments
		
        // Each cell-cell intersection can result in multiple
        // disjointed pieces if one of the cells in non-convex.
        // therefore, there may be more than one set of moments per
        // cell pair. Transform the 3 nested std::vector form to 2
        // nested std::vector form with duplicate candidate entries if
        // need be 
        int nalloc = 0;
        for (int i = 0; i < candidates.size(); ++i) {
            nalloc += moments[i].size(); // number of moment sets generated by
            //                           // intersection of target cell with
            //                           // candidate cell i
        }
        std::vector<int> candidates_dup(nalloc);
        std::vector< std::vector<double> > remap_moments(nalloc);
  
        int ninserted = 0;
        for (int i = 0; i < candidates.size(); ++i) {
            std::vector< std::vector<double> > & candidate_moments = moments[i];
            int num_moment_sets = candidate_moments.size();
            for (int j = 0; j < num_moment_sets; j++) {
                candidates_dup[ninserted] = candidates[i]; // repeated as needed 
                remap_moments[ninserted] = candidate_moments[j];
                ++ninserted;
            }
        }
        
        std::pair< std::vector<int> const &, 
                   std::vector< std::vector<double> > const & >
                source_cells_and_weights(candidates_dup,remap_moments);
        
        double remappedValue = (*remap_)(source_cells_and_weights);
        
        return remappedValue;

    }
};  // struct composerFunctor

} // namespace Portage

#endif // DRIVER_H

