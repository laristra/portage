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

template <int order> struct Interpolate_NthOrder;

template <> struct Interpolate_NthOrder<1> {
  template<int D,
           Wonton::Entity_kind ONWHAT,
           class SourceMeshType,
           class TargetMeshType,
           class SourceStateType,
           template<class, int, class, class> class InterfaceReconstructorType,
           class MatPoly_Splitter,
           class MatPoly_Clipper,
           class CoordSys>
  using Interpolate = Portage::Interpolate_1stOrder<D, ONWHAT,
                                                    SourceMeshType,
                                                    SourceStateType,
                                                    TargetMeshType,
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
           template<class, int, class, class> class InterfaceReconstructorType,
           class MatPoly_Splitter,
           class MatPoly_Clipper,
           class CoordSys>
  using Interpolate = Portage::Interpolate_2ndOrder<ONWHAT,
                                                    SourceMeshType,
                                                    TargetMeshType,
                                                    SourceStateType,
                                                    InterfaceReconstructorType,
                                                    MatPoly_Splitter,
                                                    MatPoly_Clipper,
                                                    CoordSys>;
};

}

#endif
