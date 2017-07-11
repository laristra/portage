// $Header: /cvsroot/loreli1/loreli/src/estimator/pairs.hh,v 1.5 2008/02/06 19:44:01 gad Exp $

#ifndef pairs_INCLUDED
#define pairs_INCLUDED

#include <vector>
#include <list>
#include <memory>

#include "pile.hh"
#include "lretypes.hh"

namespace lre {

  //\///////////////////////////////////////////////////////////////////////////
  // pair finding functions
  //\///////////////////////////////////////////////////////////////////////////

  /// types of pair finder
  enum contain_type{CELLS, SORT, HASHX, HASHY};

  /// Neighbor finding based on containment: driver function
  std::shared_ptr<std::vector<std::list<ulong>>> PairsFind(
      const vpile &x, const vpile &y, const vpile &h,
		  const contain_type type=CELLS, const bool half_pairs=false);

}

#endif
