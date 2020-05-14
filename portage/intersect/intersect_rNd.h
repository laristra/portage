/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef PORTAGE_INTERSECT_INTERSECT_RND_H_
#define PORTAGE_INTERSECT_INTERSECT_RND_H_

#include "wonton/support/wonton.h"

#include "portage/support/portage.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"

namespace Portage {

// A struct templated on spatial dimension so that R3D based
// intersection classes can be referred to in a uniform manner
//
// Use as
// Portage::Intersect_RND<D>::Intersect<blah,blah,blah...>
//
// If Intersect becomes a dependent name, one may have to prefix it
// with the keyword 'template', like so
// Portage::Intersect_RND<D>::template Intersect

template <int dim> struct IntersectRND;

template <> struct IntersectRND<2> {
  template<Wonton::Entity_kind ONWHAT,
           class SourceMeshType,
           class SourceStateType,
           class TargetMeshType,
           template<class, int, class, class> class InterfaceReconstructorType =
           DummyInterfaceReconstructor,
           class MatPoly_Splitter = void,
           class MatPoly_Clipper = void>
  using Intersect = Portage::IntersectR2D<ONWHAT,
                                          SourceMeshType,
                                          SourceStateType,
                                          TargetMeshType,
                                          InterfaceReconstructorType,
                                          MatPoly_Splitter,
                                          MatPoly_Clipper>;
};

template <> struct IntersectRND<3> {
  template<Wonton::Entity_kind ONWHAT,
           class SourceMeshType,
           class SourceStateType,
           class TargetMeshType,
           template<class, int, class, class> class InterfaceReconstructorType =
           DummyInterfaceReconstructor,
           class MatPoly_Splitter = void,
           class MatPoly_Clipper = void>
  using Intersect = Portage::IntersectR3D<ONWHAT,
                                          SourceMeshType,
                                          SourceStateType,
                                          TargetMeshType,
                                          InterfaceReconstructorType,
                                          MatPoly_Splitter,
                                          MatPoly_Clipper>;
};

}

#endif
