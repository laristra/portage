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

template <int dim,
          Wonton::Entity_kind ONWHAT,
          class SourceMeshType,
          class SourceStateType,
          class TargetMeshType,
          template<class, int, class, class> class InterfaceReconstructorType =
          DummyInterfaceReconstructor,
          class MatPoly_Splitter = void,
          class MatPoly_Clipper = void>
class IntersectRnD {
 public:
  using Intersector =
      typename std::conditional<dim == 2,
                                IntersectR2D<ONWHAT,
                                             SourceMeshType,
                                             SourceStateType,
                                             TargetMeshType,
                                             InterfaceReconstructorType,
                                             MatPoly_Splitter,
                                             MatPoly_Clipper>,
                                IntersectR3D<ONWHAT,
                                             SourceMeshType,
                                             SourceStateType,
                                             TargetMeshType,
                                             InterfaceReconstructorType,
                                             MatPoly_Splitter,
                                             MatPoly_Clipper>
                                >::type;

#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, dim,
                      SourceMeshType,
                      MatPoly_Splitter, MatPoly_Clipper>;

  IntersectRnD(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor> ir)
      : intersector_(source_mesh, source_state, target_mesh, num_tols, ir) {}

#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectRnD(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols)
      : intersector_(source_mesh, source_state, target_mesh, num_tols) {}

  /// \brief Set the source mesh material that we have to intersect against
  inline
  void set_material(int m) {
    intersector_.set_material(m);
  }

  /// \brief Intersect control volume of a target entity with control volumes of a set of source entities
  /// \param[in] tgt_entity  Entity of target mesh to intersect
  /// \param[in] src_entities Entities of source mesh to intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  inline
  std::vector<Weights_t>
  operator() (const int tgt_entity, std::vector<int> const& src_entities) const {
    return intersector_(tgt_entity, src_entities);
  }

 private:
  Intersector intersector_;
};

}

#endif
