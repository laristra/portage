/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SEARCH_H
#define SEARCH_H

#include "Mesh.hh"

/*!
    \class Search search.h
    \brief Search provides...
 */

namespace Portage {

	class Search
	{
	public:

		//! Default constructor
		Search() { };

		//! Constructor with Meshes
		Search(const Jali::Mesh* sourceMesh, const Jali::Mesh* targetMesh)
			: sourceMesh_(sourceMesh), targetMesh_(targetMesh) {}

		//! Copy constructor (disabled)
		Search(const Search &) = delete;

		//! Assignment operator (disabled)
		Search & operator = (const Search &) = delete;

		//! Destructor
		~Search() {}

		/*!
		  \brief This method does...

		  \param inputMesh a pointer to a Jali Mesh that I pass in...
		  \param targetMesh a pointer to a Jali Mesh that I pass in...

		  This method does something useful...
		*/
		void search(const Jali::Entity_ID cellId, Jali::Entity_ID_List* candidates) const;

	private:

		// Aggregate data members
		const Jali::Mesh* sourceMesh_;
		const Jali::Mesh* targetMesh_;
		//    double val_;

	}; // class Search

} // namespace Portage

#endif // SEARCH_H

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
