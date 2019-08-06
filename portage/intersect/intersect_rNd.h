/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef PORTAGE_INTERSECT_INTERSECT_RND_H_
#define PORTAGE_INTERSECT_INTERSECT_RND_H_

#include "portage/support/portage.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"

namespace Portage {

template <int dim> struct IntersectRND;

template <> struct IntersectRND<2> {
  template<Wonton::Entity_kind ONWHAT,
           class SourceMeshType, class SourceStateType, class TargetMeshType,
           template<class, int, class, class> class InterfaceReconstructorType,
           class MatPoly_Splitter, class MatPoly_Clipper>
  using Intersect = Portage::IntersectR2D<ONWHAT,
                                          SourceMeshType, SourceStateType,
                                          TargetMeshType,
                                          InterfaceReconstructorType,
                                          MatPoly_Splitter, MatPoly_Clipper>;
};

template <> struct IntersectRND<3> {
  template<Wonton::Entity_kind ONWHAT,
           class SourceMeshType, class SourceStateType, class TargetMeshType,
           template<class, int, class, class> class InterfaceReconstructorType,
           class MatPoly_Splitter, class MatPoly_Clipper>
  using Intersect = Portage::IntersectR3D<ONWHAT,
                                          SourceMeshType, SourceStateType,
                                          TargetMeshType,
                                          InterfaceReconstructorType,
                                          MatPoly_Splitter, MatPoly_Clipper>;
};

}

#endif
