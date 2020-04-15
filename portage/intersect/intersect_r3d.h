/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERSECT_INTERSECT_R3D_H_
#define PORTAGE_INTERSECT_INTERSECT_R3D_H_

#include <array>
#include <stdexcept>
#include <vector>
#include <algorithm>

// portage includes
extern "C" {
#include "wonton/intersect/r3d/r3d.h"
}
#include "portage/support/portage.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/intersect/intersect_polys_r3d.h"

#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/driver/CellMatPoly.h"
#include "tangram/driver/driver.h"
#include "tangram/support/MatPoly.h"
#endif

namespace Portage {

///
/// \class IntersectR3D  3-D intersection algorithm
///
/// In this routine, we will utilize the fact that R3d can intersect a
/// non-convex polyhedron with planar faces with a convex polyhedron
/// (more precisely, R3D can clip a non-convex polyhedron with a set of
/// planes). We are given a target polyhedron with possibly non-planar
/// faces to be intersected with a list of source polyhedra again with
/// possibly non-planar faces. We will convert the target polyhedron
/// into a set of convex polyhedra using a symmetric tetrahedral
/// decomposition (24 tets for a hex). We will convert each source
/// polyhedron into a faceted non-convex polyhedron where each facet is
/// a triangle and therefore planar.
///
/// If this class is being adapted for use with a different intersector
/// and it can only intersect convex polyhedra, both target and source
/// polyhedra may have to be decomposed into tets
///

template <Entity_kind on_what, class SourceMeshType,
          class SourceStateType, class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType =
          DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void>
class IntersectR3D {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor3D =
      Tangram::Driver<InterfaceReconstructorType, 3, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:
#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR3D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor3D> ir,
               bool rectangular_mesh = false)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir),
        rectangular_mesh_(rectangular_mesh), num_tols_(num_tols) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR3D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               bool rectangular_mesh = false)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), rectangular_mesh_(rectangular_mesh),
        num_tols_(num_tols) {}

  /// \brief Set the source mesh material that we have to intersect against

  void set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect a control volume of a target_entity with control volumes of a set of source entities
  /// \param[in] tgt_entity Entity of target mesh to intersect
  /// \param[in] src_entities Entity of source cells to intersect against
  /// \return vector of Weights_t structure containing moments of intersection
  ///

  std::vector<Weights_t> operator() (int tgt_entity, std::vector<int> const& src_entities) const {
    throw std::runtime_error("not implemented for this entity type");
  }



  IntersectR3D() = delete;

  /// Assignment operator (disabled)
  IntersectR3D & operator = (const IntersectR3D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor3D> interface_reconstructor;
#endif
  bool rectangular_mesh_ = false;
  int matid_ = -1;
  NumericTolerances_t num_tols_ {};
};


// Specialization of Intersect3D class for cells

template <class SourceMeshType, class SourceStateType,
          class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter, class Matpoly_Clipper>
class IntersectR3D<Entity_kind::CELL, SourceMeshType, SourceStateType, TargetMeshType,
                   InterfaceReconstructorType,
                   Matpoly_Splitter, Matpoly_Clipper> {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor3D =
      Tangram::Driver<InterfaceReconstructorType, 3, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR3D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor3D> ir,
               bool rectangular_mesh = false)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir),
        rectangular_mesh_(rectangular_mesh), num_tols_(num_tols) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR3D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               bool rectangular_mesh = false)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), rectangular_mesh_(rectangular_mesh),
        num_tols_(num_tols) {}


  /// \brief Set the source mesh material that we have to intersect against

  void set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect a cell with a set of candidate cells
  /// \param[in] tgt_cell cell of target mesh to intersect
  /// \param[in] src_cells list of source cells to intersect against
  /// \return vector of Weights_t structures containing moments of intersection
  ///

  std::vector<Weights_t> operator() (const int tgt_cell,
                                     const std::vector<int>& src_cells) const {

    std::vector<std::array<Point<3>, 4>> target_tet_coords;

    // We should avoid any decomposition for cells of a rectangular
    // mesh but for now we will decompose the target all the time

    targetMeshWrapper.decompose_cell_into_tets(tgt_cell, &target_tet_coords,
                                               rectangular_mesh_);

    // CAN MAKE THIS INTO A THRUST::TRANSFORM CALL
    int nsrc = src_cells.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_cells[i];

      Weights_t & this_wt = sources_and_weights[ninserted];
      this_wt.entityID = s;

#ifdef PORTAGE_HAS_TANGRAM
      int nmats = sourceStateWrapper.cell_get_num_mats(s);
      std::vector<int> cellmats;
      sourceStateWrapper.cell_get_mats(s, &cellmats);

      if (!nmats || (matid_ == -1) || (nmats == 1 && cellmats[0] == matid_)) {
        // ---------- Intersection with pure cell ---------------
        // nmats == 0 -- no materials ==> single material
        // matid_ == -1 -- intersect with mesh not a particular material
        // nmats == 1 && cellmats[0] == matid -- intersection with pure cell
        //                                       containing matid

        facetedpoly_t srcpoly;
        sourceMeshWrapper.cell_get_facetization(s, &srcpoly.facetpoints,
                                                &srcpoly.points);

        this_wt.weights = intersect_polys_r3d(srcpoly, target_tet_coords,
                                              num_tols_);

      } else if (std::find(cellmats.begin(), cellmats.end(), matid_) !=
                 cellmats.end()) {
        // mixed cell containing this material - intersect with
        // polygon approximation of this material in the cell
        // (obtained from interface reconstruction)

        Tangram::CellMatPoly<3> const& cellmatpoly =
            interface_reconstructor->cell_matpoly_data(s);
        std::vector<Tangram::MatPoly<3>> matpolys =
            cellmatpoly.get_matpolys(matid_);

        this_wt.weights.resize(4,0.0);
        for (const auto& matpoly : matpolys) {
          facetedpoly_t srcpoly = get_faceted_matpoly(matpoly);

          std::vector<double> momvec = intersect_polys_r3d(srcpoly,
                                          target_tet_coords, num_tols_);
          for (int k = 0; k < 4; k++)
            this_wt.weights[k] += momvec[k];
        }

      }
#else
      facetedpoly_t srcpoly;
      sourceMeshWrapper.cell_get_facetization(s, &srcpoly.facetpoints,
                                              &srcpoly.points);

      this_wt.weights = intersect_polys_r3d(srcpoly, target_tet_coords,
                                            num_tols_);
#endif
      // Increment if vol of intersection > 0; otherwise, allow overwrite
      if (!this_wt.weights.empty() && this_wt.weights[0] > 0.0)
        ninserted++;
    }

    sources_and_weights.resize(ninserted);
    return sources_and_weights;
  }


  IntersectR3D() = delete;

  /// Assignment operator (disabled)
  IntersectR3D & operator = (const IntersectR3D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor3D> interface_reconstructor;
#endif
  bool rectangular_mesh_ = false;
  int matid_ = -1;
  NumericTolerances_t num_tols_ {};
};


// Specialization of IntersectR3D class for nodes

template <class SourceMeshType, class SourceStateType,
          class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter, class Matpoly_Clipper>
class IntersectR3D<Entity_kind::NODE, SourceMeshType, SourceStateType, TargetMeshType,
                   InterfaceReconstructorType,
                   Matpoly_Splitter, Matpoly_Clipper> {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor3D =
      Tangram::Driver<InterfaceReconstructorType, 3, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR3D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor3D> ir,
               bool rectangular_mesh = false)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir),
        rectangular_mesh_(rectangular_mesh), num_tols_(num_tols) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR3D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               bool rectangular_mesh = false)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), rectangular_mesh_(rectangular_mesh),
        num_tols_(num_tols) {}

  /// \brief Set the source mesh material that we have to intersect against

  void set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect a control volume corresponding to a target node
  /// with a set of control volumes corresponding to candidate source
  /// nodes
  /// \param[in] tgt_node   Target mesh node whose control volume we consider
  /// \param[in] src_nodes  List of source nodes whose control volumes we will intersect against
  /// \return vector of Weights_t structures containing moments of intersection
  ///

  std::vector<Weights_t> operator() (const int tgt_node,
                                     const std::vector<int>& src_nodes) const {

    // for debug
    Point<3> tgtxyz;
    targetMeshWrapper.node_get_coordinates(tgt_node, &tgtxyz);

    std::vector<std::array<Point<3>, 4>> target_tet_coords;

    // We should avoid any decomposition for duall cells of a
    // rectangular mesh but for now we will decompose the target all
    // the time

    targetMeshWrapper.dual_wedges_get_coordinates(tgt_node, &target_tet_coords);


    // CAN MAKE THIS INTO A THRUST TRANSFORM CALL
    int nsrc = src_nodes.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_nodes[i];

      facetedpoly_t srcpoly;
      sourceMeshWrapper.dual_cell_get_facetization(s, &srcpoly.facetpoints,
                                                   &srcpoly.points);

      Weights_t & this_wt = sources_and_weights[ninserted];
      this_wt.entityID = s;
      this_wt.weights = intersect_polys_r3d(srcpoly, target_tet_coords,
                                            num_tols_);

      // Increment if vol of intersection > 0; otherwise, allow overwrite
      if (!this_wt.weights.empty() && this_wt.weights[0] > 0.0)
        ninserted++;
    }

    sources_and_weights.resize(ninserted);
    return sources_and_weights;
  }


  IntersectR3D() = delete;

  /// Assignment operator (disabled)
  IntersectR3D & operator = (const IntersectR3D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor3D> interface_reconstructor;
#endif
  bool rectangular_mesh_ = false;
  int matid_ = -1;
  NumericTolerances_t num_tols_ {};
};  // class IntersectR3D


}  // namespace Portage

#endif  // PORTAGE_INTERSECT_INTERSECT_R3D_H_
