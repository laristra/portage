/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

#include<algorithm>
#include<vector>

#include "Mesh.hh"
#include "portage/state/state.h"

/*!
    \class Driver driver.h
    \brief Driver provides the API to mapping from one mesh to another.
 */

namespace Portage {

class Driver
{
public:
  
	//! Constructor uses established meshes for remap
	Driver(Jali::Mesh const & sourceMesh,State const & sourceState,
		   Jali::Mesh const & targetMesh, State & targetState) 
		: sourceMesh_(sourceMesh), sourceState_(sourceState),
		  targetMesh_(targetMesh), targetState_(targetState) {}
  
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

	template <typename SearchType, typename IsectType, typename RemapType>
	struct composerFunctor
	{
		const SearchType* s_;
		const IsectType* i_;
		const RemapType* r_;
		composerFunctor(const SearchType* s, const IsectType* i, const RemapType* r)
			: s_(s), i_(i), r_(r) { }

		
		{
			
		}

		double operator()(Jali::Entity_ID const targetCellIndex)
		{
			// Search for candidates and return their cell indices
			Jali::Entity_ID_List* candidates;
			s_->search(targetCellIndex, candidates);

			// Intersect wants a vector of pairs of (x,y) coordinates for
			// each candidate.
			std::vector<std::vector<std::pair<double,double> > > 
				candidatesNodesCoords(candidates->size());
		    std::transform(candidates->begin(), candidates->end(),
			    		   candidatesNodesCoords.begin(),
						   [&](const Jali::Entity_ID cellID) -> std::vector<std::pair<double,double> >
						   {
							   std::vector<JaliGeometry::Point> cellNodes;
							   targetMesh_->cell_get_coordinates(cellID, &cellNodes);
							   std::vector<std::pair<double,double> > cellNodeCoords(cellNodes.size());
							   std::transform(cellNodes.begin(), cellNodes.end(),
											  cellNodeCoords.begin(),
											  [](const JaliGeometry::Point node)
											  {
												  return std::make_pair(node.x(), node.y());
											  });
							   return cellNodeCoords;});
	}



private:
  
	Jali::Mesh const & sourceMesh_;
	Jali::Mesh const & targetMesh_;
	State const & sourceState_;
	State & targetState_;
	std::vector<std::string> remap_var_names_;

}; // class Driver


} // namespace Portage

#endif // DRIVER_H

/*--------------------------------------------------------------------------~-*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
