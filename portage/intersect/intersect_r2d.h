/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef PORTAGE_INTERSECT_INTERSECT_R2D_H_
#define PORTAGE_INTERSECT_INTERSECT_R2D_H_

#include <array>
#include <stdexcept>
#include <vector>
#include <algorithm>

// portage includes
extern "C" {
#include "wonton/intersect/r3d/r2d.h"
}
#include "portage/support/portage.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/intersect/intersect_polys_r2d.h"

#ifdef HAVE_TANGRAM
#include "tangram/driver/CellMatPoly.h"
#include "tangram/driver/driver.h"
#include "tangram/support/MatPoly.h"
#endif 

#include "wonton/support/Point.h"

namespace Portage {

///
/// \class IntersectR2D  2-D intersection algorithm


template <Entity_kind on_what, class SourceMeshType,
          class SourceStateType, class TargetMeshType,
          template<class, int, class, class> class InterfaceReconstructorType =
          DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void>
class IntersectR2D {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor2D =
      Tangram::Driver<InterfaceReconstructorType, 2, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef HAVE_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               std::shared_ptr<InterfaceReconstructor2D> ir)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh) {}

  /// \brief Set the source mesh material that we have to intersect against

  int set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect control volume of a target entity with control volumes of a set of source entities
  /// \param[in] tgt_entity  Entity of target mesh to intersect
  /// \param[in] src_entities Entities of source mesh to intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  std::vector<Weights_t>
  operator() (const int tgt_entity, const std::vector<int> src_entities) const {
    std::cerr << "IntersectR3D not implemented for this entity type" <<
        std::endl;
  }

  IntersectR2D() = delete;

  /// Assignment operator (disabled)
  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
  int matid_ = -1;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
};  // class IntersectR2D



//////////////////////////////////////////////////////////////////////////////
// Specialization of Intersect2D class for cells

template <class SourceMeshType, class SourceStateType,
          class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter, class Matpoly_Clipper>
class IntersectR2D<Entity_kind::CELL, SourceMeshType, SourceStateType, TargetMeshType,
                   InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper> {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor2D =
      Tangram::Driver<InterfaceReconstructorType, 2, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef HAVE_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               std::shared_ptr<InterfaceReconstructor2D> ir)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh) {}

  /// \brief Set the source mesh material that we have to intersect against

  int set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect target cell with a set of source cell
  /// \param[in] tgt_entity  Cell of target mesh to intersect
  /// \param[in] src_entities List of source cells to intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  std::vector<Weights_t>
  operator() (const int tgt_cell, const std::vector<int> src_cells) const {
    std::vector<Wonton::Point<2>> target_poly;
    targetMeshWrapper.cell_get_coordinates(tgt_cell, &target_poly);

    int nsrc = src_cells.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_cells[i];

      Weights_t & this_wt = sources_and_weights[ninserted];
      this_wt.entityID = s;

#ifdef HAVE_TANGRAM
      int nmats = sourceStateWrapper.cell_get_num_mats(s);
      std::vector<int> cellmats;
      sourceStateWrapper.cell_get_mats(s, &cellmats);

      if (!nmats || (matid_ == -1) || (nmats == 1 && cellmats[0] == matid_)) {
        // ---------- Intersection with pure cell ---------------
        // nmats == 0 -- no materials ==> single material
        // matid_ == -1 -- intersect with mesh not a particular material
        // nmats == 1 && cellmats[0] == matid -- intersection with pure cell
        //                                       containing matid

        std::vector<Wonton::Point<2>> source_poly;
        sourceMeshWrapper.cell_get_coordinates(s, &source_poly);

        this_wt.weights = intersect_polys_r2d(source_poly, target_poly);

      } else {  // multi-material case
        // How can I check that I didn't get DummyInterfaceReconstructor
        //
        // static_assert(InterfaceReconstructorType<SourceMeshType, 2> !=
        //               DummyInterfaceReconstructor<SourceMeshType, 2>);

        assert(interface_reconstructor != nullptr);  // cannot be nullptr

        if (std::find(cellmats.begin(), cellmats.end(), matid_) !=
            cellmats.end()) {
          // mixed cell containing this material - intersect with
          // polygon approximation of this material in the cell
          // (obtained from interface reconstruction)

          Tangram::CellMatPoly<2> const& cellmatpoly =
              interface_reconstructor->cell_matpoly_data(s);
          std::vector<Tangram::MatPoly<2>> matpolys =
              cellmatpoly.get_matpolys(matid_);

          this_wt.weights.resize(3, 0.0);
          for (int j = 0; j < matpolys.size(); j++) {
            std::vector<Wonton::Point<2>> tpnts = matpolys[j].points();
            std::vector<Wonton::Point<2>> source_poly;
            source_poly.reserve(tpnts.size());
            for (auto const & p : tpnts) source_poly.push_back(p);

            std::vector<double> momvec = intersect_polys_r2d(source_poly,
                                                             target_poly);
            for (int k = 0; k < 3; k++)
              this_wt.weights[k] += momvec[k];
          }
        }
      }
#else
      std::vector<Wonton::Point<2>> source_poly;
      sourceMeshWrapper.cell_get_coordinates(s, &source_poly);

      this_wt.weights = intersect_polys_r2d(source_poly, target_poly);
#endif

      // Increment if vol of intersection > 0; otherwise, allow overwrite
      if (this_wt.weights.size() && this_wt.weights[0] > 0.0)
        ninserted++;
    }

    sources_and_weights.resize(ninserted);
    return sources_and_weights;
  }

  IntersectR2D() = delete;

  /// Assignment operator (disabled)
  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
  int matid_ = -1;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
};  // class IntersectR2D




//////////////////////////////////////////////////////////////////////////////
// Specialization of Intersect2D class for nodes

template <class SourceMeshType, class SourceStateType,
          class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter, class Matpoly_Clipper>
class IntersectR2D<Entity_kind::NODE, SourceMeshType, SourceStateType, TargetMeshType,
                   InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper> {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor2D =
      Tangram::Driver<InterfaceReconstructorType, 2, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef HAVE_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               std::shared_ptr<InterfaceReconstructor2D> ir)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir) {}
#endif


  /// Constructor WITHOUT interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh) {}

  /// \brief Set the source mesh material that we have to intersect against

  int set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect control volume of a target node with control volumes of
  /// a set of source nodes
  /// \param[in] tgt_node  Target mesh node whose control volume we consider
  /// \param[in] src_nodes List of source nodes whose control volumes we will intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  std::vector<Weights_t>
  operator() (const int tgt_node, const std::vector<int> src_nodes) const {
    std::vector<Wonton::Point<2>> target_poly;
    targetMeshWrapper.dual_cell_get_coordinates(tgt_node, &target_poly);

    int nsrc = src_nodes.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_nodes[i];
      std::vector<Wonton::Point<2>> source_poly;
      sourceMeshWrapper.dual_cell_get_coordinates(s, &source_poly);

      Weights_t & this_wt = sources_and_weights[ninserted];
      this_wt.entityID = s;
      this_wt.weights = intersect_polys_r2d(source_poly, target_poly);

      // Increment if vol of intersection > 0; otherwise, allow overwrite
      if (this_wt.weights.size() && this_wt.weights[0] > 0.0)
        ninserted++;
    }

    sources_and_weights.resize(ninserted);
    return sources_and_weights;
  }

  IntersectR2D() = delete;

  /// Assignment operator (disabled)

  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
  int matid_ = -1;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
};  // class IntersectR2D

} // namespace Portage

#endif // PORTAGE_INTERSECT_INTERSECT_R2D_H_
