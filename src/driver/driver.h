/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

#include<algorithm>
#include<vector>
#include<iterator>

#include "Mesh.hh"   // Jali mesh header
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

/*!
    \class Driver driver.h
    \brief Driver provides the API to mapping from one mesh to another.
 */

// When we get Portage::Entity_kind, Portage::Entity_ID definitions
// into a header file in the appropriate place, we can make everything
// in the method free of Jali references except the Jali_Mesh_Wrapper
// and Jali_State_Wrapper references

namespace Portage {

class Driver
{
  public:
  
    //! Constructor - takes in wrapper classes for source/target mesh and state

    Driver(Jali_Mesh_Wrapper const & sourceMesh, 
           Jali_State_Wrapper const & sourceState,           
           Jali_Mesh_Wrapper const & targetMesh,
           Jali_State_Wrapper & targetState) 
            : source_mesh_(sourceMesh), source_state_(sourceState),
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
    void run();


private:
  
    Jali_Mesh_Wrapper  const & source_mesh_;
    Jali_Mesh_Wrapper const & target_mesh_;
    Jali_State_Wrapper const & source_state_;
    Jali_State_Wrapper & target_state_;
    std::vector<std::string> remap_var_names_;

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
		
        // RVG - Need to reconcile how moments are returned and
        // how remap expects them - for now create a dummy vector
        // that conforms to what remap functor expects assuming
        // that the intersection of each pair of cells results in
        // only domain of intersection and that we only care about
        // the 0th order moments (area)
        
        std::vector<double> remap_moments(candidates.size(),0.0);
        for (int i = 0; i < candidates.size(); ++i) {
            std::vector< std::vector<double> > & candidate_moments = moments[i];
            // BUG: candidate_moments might be empty; sum over it? amh
            std::vector<double> & piece_moments = candidate_moments[0];
            remap_moments[i] = piece_moments[0];
        }
        
        std::pair< std::vector<int> const &, std::vector<double> const & >
                source_cells_and_weights(candidates,remap_moments);
        
        double remappedValue = (*remap_)(source_cells_and_weights);
        
        return remappedValue;
    }
};  // struct composerFunctor

} // namespace Portage

#endif // DRIVER_H

/*--------------------------------------------------------------------------~-*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:nil
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
