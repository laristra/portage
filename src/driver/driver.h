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
#include "portage/state/state.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

/*!
    \class Driver driver.h
    \brief Driver provides the API to mapping from one mesh to another.
 */

namespace Portage {

class Driver
{
  public:
  
    //! Constructor uses established meshes for remap
    Driver(Jali::Mesh const & sourceMesh, State const & sourceState,
           Jali::Mesh const & targetMesh, State & targetState) 
            : source_mesh_wrapper_(Jali_Mesh_Wrapper(sourceMesh)), 
              sourceMesh_(sourceMesh),
              sourceState_(sourceState),
              target_mesh_wrapper_(Jali_Mesh_Wrapper(targetMesh)), 
              targetMesh_(targetMesh),
              targetState_(targetState) 
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
  
    Jali_Mesh_Wrapper  const & source_mesh_wrapper_;
    Jali_Mesh_Wrapper const & target_mesh_wrapper_;
    Jali::Mesh const & sourceMesh_;
    Jali::Mesh const & targetMesh_;
    State const & sourceState_;
    State & targetState_;
    std::vector<std::string> remap_var_names_;

}; // class Driver

// This functor is used inside a std::transform inside Driver::run
template <typename SearchType, typename IsectType, typename RemapType,
          typename SourceMeshWrapper, typename TargetMeshWrapper>
struct composerFunctor
{
        const SearchType* search_;
		const IsectType* intersect_;
		const RemapType* remap_;
		const SourceMeshWrapper & sourceMesh_;
		const TargetMeshWrapper & targetMesh_;
		const std::string remap_var_name_;
		//----------------------------------------
		composerFunctor(const SearchType* s, const IsectType* i, 
						const RemapType* r,
						const SourceMeshWrapper & sourceMesh, 
						const TargetMeshWrapper & targetMesh,
						const std::string remap_var_name)
			: search_(s), intersect_(i), remap_(r), sourceMesh_(sourceMesh), 
			  targetMesh_(targetMesh), remap_var_name_(remap_var_name) { }


		double operator()(int const targetCellIndex)
		{
			// Search for candidates and return their cells indices
            std::vector<int> candidates;
			search_->search(targetCellIndex, &candidates);

			std::vector<std::vector<std::vector<double> > > moments(candidates.size());
            for (int i=0;i<candidates.size();i++){
                moments[i] = (*intersect_)(candidates[i], targetCellIndex);             
            }

			// Remap
			
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

			// Compute the remap value from all the candidate cells and weights

			double remappedValue = (*remap_)(source_cells_and_weights);

			return remappedValue;
		}
	};

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
