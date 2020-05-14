/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef PORTAGE_INTERPOLATE_INTERPOLATE_RND_H_
#define PORTAGE_INTERPOLATE_INTERPOLATE_RND_H_

#include "wonton/support/wonton.h"

#include "portage/support/portage.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"

namespace Portage {

// A struct templated on interpolation order so that interpolate
// classes can be referred to in a uniform manner
//
// Use as
// Portage::Interpolate_NthOrder<O>::Interpolate<blah,blah...>
// where O is a compile time parameter
//
// If Interpolate becomes a dependent name, one may have to prefix it
// with the keyword 'template', like so
// Portage::Interpolate_NthOrder<D>::template Interpolate

template <int order>
struct Interpolate_NthOrder;

template <> struct Interpolate_NthOrder<1> {
  template<int D,
           Wonton::Entity_kind ONWHAT,
           class SourceMeshType,
           class TargetMeshType,
           class SourceStateType,
           class TargetStateType = SourceStateType,
           typename T = double,
           template<class, int, class, class>
             class InterfaceReconstructorType = DummyInterfaceReconstructor,
           class MatPoly_Splitter = void,
           class MatPoly_Clipper = void,
           class CoordSys = Wonton::DefaultCoordSys>
  using Interpolate = Portage::Interpolate_1stOrder<D, ONWHAT,
                                                    SourceMeshType,
                                                    TargetMeshType,
                                                    SourceStateType,
                                                    TargetStateType,
                                                    T,
                                                    InterfaceReconstructorType,
                                                    MatPoly_Splitter,
                                                    MatPoly_Clipper,
                                                    CoordSys>;
};

template <> struct Interpolate_NthOrder<2> {
  template<int D,
           Wonton::Entity_kind ONWHAT,
           class SourceMeshType,
           class TargetMeshType,
           class SourceStateType,
           class TargetStateType = SourceStateType,
           typename T = double,
           template<class, int, class, class>
             class InterfaceReconstructorType = DummyInterfaceReconstructor,
           class MatPoly_Splitter = void,
           class MatPoly_Clipper = void,
           class CoordSys = Wonton::DefaultCoordSys>
  using Interpolate = Portage::Interpolate_2ndOrder<D, ONWHAT,
                                                    SourceMeshType,
                                                    TargetMeshType,
                                                    SourceStateType,
                                                    TargetStateType,
                                                    T,
                                                    InterfaceReconstructorType,
                                                    MatPoly_Splitter,
                                                    MatPoly_Clipper,
                                                    CoordSys>;
};

}

#endif
