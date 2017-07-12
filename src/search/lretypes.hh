// $Header: /cvsroot/loreli1/loreli/src/estimator/types.hh,v 1.1 2006/06/21 00:02:41 gad Exp $

#ifndef types_INCLUDED
#define types_INCLUDED

#include <valarray>

namespace lre {
  using std::valarray;

  typedef unsigned int uint;                                  ///< convenience definition
  typedef unsigned long ulong;                                ///< convenience definition
  typedef valarray< ulong > vulong;                           ///< convenience definition
  typedef valarray< valarray< ulong > > vvulong;              ///< convenience definition
  typedef valarray< valarray< valarray< ulong > > > vvvulong; ///< convenience definition
}

#endif
