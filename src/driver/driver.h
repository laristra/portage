/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

#include<algorithm>
#include<vector>
#include<iterator>

#include "portage/support/portage.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersectClipper.h"
#include "portage/remap/remap_1st_order.h"

namespace Portage {

/*! 
    \class MeshWrapperDual driver.h
    \brief Wrapper for dual mesh
 */

class MeshWrapperDual {
public:
    MeshWrapperDual(const Jali_Mesh_Wrapper &w) : w_(w) {}
    int num_owned_cells() const { return w_.num_owned_nodes(); }
    int num_ghost_cells() const { return w_.num_ghost_nodes(); }
    void cell_get_coordinates(int const cellid,
            std::vector<std::pair<double,double> > *xylist) const {
        w_.dual_cell_get_coordinates(cellid, xylist);
    }

  counting_iterator begin(Entity_kind const entity) const {
        if (entity == NODE) return w_.begin(CELL);
        return w_.begin(NODE);
    }

  counting_iterator end(Entity_kind const entity) const {
        if (entity == NODE) return w_.end(CELL);
        return w_.end(NODE);
    }
    std::vector<std::pair<double, double> > cellToXY(int cellID) const{
        std::vector<std::pair<double, double> > cellPoints;
        cell_get_coordinates(cellID, &cellPoints);
        return cellPoints;
    }

private:
    const Jali_Mesh_Wrapper &w_;
};


// Forward definition
template <typename SearchType, typename IsectType, typename RemapType>
struct composerFunctor;


/*!
    \class Driver driver.h
    \brief Driver provides the API to mapping from one mesh to another.
 */

template <class Mesh_Wrapper>
class Driver
{
  public:
  
    //! Constructor - takes in wrapper classes for source/target mesh and state

    Driver(Entity_kind remapEntity, 
           Mesh_Wrapper const & sourceMesh, 
           Jali_State_Wrapper const & sourceState,           
           Mesh_Wrapper const & targetMesh,
           Jali_State_Wrapper & targetState) 
            : remap_entity_(remapEntity), 
              source_mesh_(sourceMesh), source_state_(sourceState),
              target_mesh_(targetMesh), target_state_(targetState) 
    {}
  
    //! Copy constructor (disabled)
    Driver(const Driver &) = delete;
  
    //! Assignment operator (disabled)
    Driver & operator = (const Driver &) = delete;
  
    //! Destructor
    ~Driver() {}
  
  
    // Specify the names of the variables to be remapped

    void set_remap_var_names(std::vector<std::string> &remap_var_names_in) {
        remap_var_names_ = remap_var_names_in;
    }

    // Get the names of the variables to be remapped

    std::vector<std::string> remap_var_names() {
        return remap_var_names_;
    }

    /*!
      \brief This method calls the search, intersect, and remap routines
      needed to map one mesh to another.
    */
    void run()
    {
        std::printf("in Driver::run()...\n");
       
        // Get an instance of the desired search algorithm type
        const SearchKDTree<Mesh_Wrapper,Mesh_Wrapper>
                search(source_mesh_, target_mesh_);
        
        // Get an instance of the desired intersect algorithm type
        const IntersectClipper<Mesh_Wrapper, Mesh_Wrapper>
                intersect{source_mesh_, target_mesh_};
        
        // Eventually put this in a loop over remap variable names as well
        // Assume for now that we are only doing cell-based remap
        const Remap_1stOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind>
                remap(source_mesh_, source_state_, remap_entity_, remap_var_names_[0]);
        
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

        // Create a cellIndices vector and populates with a sequence of
        // ints starting at 0.  
        composerFunctor<SearchKDTree<Mesh_Wrapper,Mesh_Wrapper>,
            IntersectClipper<Mesh_Wrapper, Mesh_Wrapper>,
            Remap_1stOrder<Mesh_Wrapper,Jali_State_Wrapper,Entity_kind> >
                composer(&search, &intersect, &remap, remap_var_names_[0]);

        // This populates targetField with the doubles returned from
        // the final remap For some reason (perhaps some code outside
        // portage is pulling in a "using thrust" command), the
        // compiler is not able to disambiguate Portage::transform and
        // thrust::transform here. So, be explicit that we want
        // Portage::transform

        Portage::transform((counting_iterator)target_mesh_.begin(CELL),
                  (counting_iterator)target_mesh_.end(CELL),
                  target_field,composer);

    }


private:
  
    Mesh_Wrapper  const & source_mesh_;
    Mesh_Wrapper const & target_mesh_;
    Jali_State_Wrapper const & source_state_;
    Jali_State_Wrapper & target_state_;
    std::vector<std::string> remap_var_names_;
    Entity_kind remap_entity_;

}; // class Driver

// This functor is used inside a std::transform inside Driver::run

// Should we move composerFunctor into the .cc file? - RVG

template <typename SearchType, typename IsectType, typename RemapType>
struct composerFunctor
{
    const SearchType* search_;
    const IsectType* intersect_;
    const RemapType* remap_;
    const std::string remap_var_name_;
    //----------------------------------------
   
    composerFunctor(const SearchType* searcher, const IsectType* intersecter, 
                    const RemapType* remapper, const std::string remap_var_name)
			: search_(searcher), intersect_(intersecter), remap_(remapper), 
              remap_var_name_(remap_var_name) { }

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

