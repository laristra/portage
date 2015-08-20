/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

#include<algorithm>
#include<vector>
#include<iterator>

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


private:
  
	Jali::Mesh const & sourceMesh_;
	Jali::Mesh const & targetMesh_;
	State const & sourceState_;
	State & targetState_;
	std::vector<std::string> remap_var_names_;

}; // class Driver

	template <typename SearchType, typename IsectType, typename RemapType>
	struct composerFunctor
	{
		const SearchType* s_;
		const IsectType* i_;
		const RemapType* r_;
		// CMM: This seems redundant...
		const Jali::Mesh* sourceMesh_;
		const Jali::Mesh* targetMesh_;
		const std::string remap_var_name_;
		//----------------------------------------
		composerFunctor(const SearchType* s, const IsectType* i, const RemapType* r,
						const Jali::Mesh* sourceMesh, const Jali::Mesh* targetMesh,
						const std::string remap_var_name)
			: s_(s), i_(i), r_(r), sourceMesh_(sourceMesh), targetMesh_(targetMesh),
		      remap_var_name_(remap_var_name) { }

		// CMM: this should probably be somewhere else - perhaps part of Jali?
		struct pointToXY
		{
			pointToXY() { }
			std::pair<double,double> operator()(const JaliGeometry::Point point)
			{
				return std::make_pair(point.x(), point.y());
			}
		};

		double operator()(Jali::Entity_ID const targetCellIndex)
		{
			std::cout << "target cell index " << targetCellIndex << std::endl;
			// Search for candidates and return their cell indices
			Jali::Entity_ID_List candidates;
			s_->search(targetCellIndex, &candidates);

			// Get the target cell's (x,y) coordinates from the Jali Point datastructure
			// CMM: do I really need this? without explicit size of targetCellPoints vector
			//      the compiler either complains or the executable err's
			Jali::Entity_ID_List nodes;
			int numnodes;
			targetMesh_->cell_get_nodes(targetCellIndex, &nodes);
			numnodes = nodes.size();
			std::vector<JaliGeometry::Point> targetCellPoints(numnodes);
			targetMesh_->cell_get_coordinates(targetCellIndex, &targetCellPoints);
			std::vector<std::pair<double, double> > targetCellCoords(numnodes);
			std::transform(targetCellPoints.begin(), targetCellPoints.end(),
						   targetCellCoords.begin(),pointToXY());

			// Intersect routine wants candidates' node coordinates
			// First, get the Jali Points for each candidate cells
			std::vector<std::vector<JaliGeometry::Point> > candidateCellsPoints(candidates.size());
			std::transform(candidates.begin(), candidates.end(),
						   candidateCellsPoints.begin(),
						   // given a candidate cell in the sourceMesh, get its Points
						   [&](Jali::Entity_ID candidateCellIndex) -> std::vector<JaliGeometry::Point>
						   {
							   Jali::Entity_ID_List nodes;
							   sourceMesh_->cell_get_nodes(candidateCellIndex, &nodes);
							   std::vector<JaliGeometry::Point> ret(nodes.size());
							   sourceMesh_->cell_get_coordinates(candidateCellIndex, &ret);
							   return ret;
						   }
						   );
			// Now get the (x,y) coordinates from each Point for each candidate cell
			std::vector<std::vector<std::pair<double, double> > > 
				candidateCellsCoords(candidates.size());
			std::transform(candidateCellsPoints.begin(), candidateCellsPoints.end(),
						   candidateCellsCoords.begin(),
						   [&](std::vector<JaliGeometry::Point> points) ->
						   std::vector<std::pair<double, double> >
						   {
							   std::vector<std::pair<double, double> > ret(points.size());
							   std::transform(points.begin(), points.end(),
											  ret.begin(),
											  pointToXY());
							   return ret;
						   });

			// Calculate the intersection of each candidate with the target Cell
			// For each polygon/polygon intersection, IntersectClipper returns a 
			// std::vector<std::vector<double>>
			std::vector<std::vector<std::vector<double> > > moments(candidates.size());
			// CMM: To do a std::transform here instead, we need to use std::tuple's (I think)
			//      both here and in IntersectClipper::operator()
			for (int i = 0; i < candidateCellsCoords.size(); i++)
				{
					// awkward syntax...
					moments[i] = (*i_)(candidateCellsCoords[i], targetCellCoords);
				}

			// Remap
			double remappedValue = r_->remap(remap_var_name_, targetCellIndex, candidates, moments);
						   
			
			return 0.0;
																				   
		}
	};




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
