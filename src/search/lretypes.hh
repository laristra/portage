/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.
*/

#ifndef types_INCLUDED
#define types_INCLUDED

#include <valarray>

namespace Portage {
namespace Meshfree {
namespace Pairs {
  using std::valarray;

  typedef unsigned int uint;                                  ///< convenience definition
  typedef unsigned long ulong;                                ///< convenience definition
  typedef valarray< ulong > vulong;                           ///< convenience definition
  typedef valarray< valarray< ulong > > vvulong;              ///< convenience definition
  typedef valarray< valarray< valarray< ulong > > > vvvulong; ///< convenience definition
}
}
}

#endif
