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
#ifdef DEBUG
#include <sstream>
#endif

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

using Wonton::Point;
using Wonton::Vector;
using Wonton::SIDE;
using Wonton::WEDGE;
using Wonton::ALL;

///
/// \class IntersectR3D  3-D intersection algorithm
///
/// In this routine, we will utilize the fact that R3d can intersect a
/// non-convex polyhedron with planar faces with a convex polyhedron
/// (more precisely, R3D can clip a non-convex polyhedron with a set
/// of planes). We are given a target polyhedron with possibly
/// non-planar faces to be intersected with a list of source polyhedra
/// again with possibly non-planar faces. We will convert the target
/// polyhedron into a set of convex polyhedra using a symmetric
/// tetrahedral decomposition (24 tets for a hex) or a custom 5/6 tet
/// decomposition for rectangular meshes. We will convert each source
/// polyhedron into a faceted non-convex polyhedron where each facet
/// is a triangle and therefore planar.
///
/// If this class is being adapted for use with a different intersector
/// and it can only intersect convex polyhedra, both target and source
/// polyhedra may have to be decomposed into tets
///

// There are few criteria that could be used to evaluate whether the
// polygons/polyhedra are suitable for use in the intersection using
// R2D/R3D's clipping routines:
//
// 1. "VALIDITY" OF SOURCE CELL: The source polyhedron (from a cell, a
// dual cell or a material polyhedron) can be non-convex since R2D/R3D
// can clip non-convex cells. We may, however, want to verify that the
// source polyhedron is "valid" before computing with it (although the
// way R2D/R3D are written and the interpolate routines are written,
// there may be nothing to preclude us from computing with negative
// intersection volumes). We can test the validity of a source
// polyhedron by (a) checking that its volume computed using
// divergence theorem (by integrating over its boundary) is positive
// (b) verifying that it is star-convex i.e. allows a decomposition
// into positive volume simplices.
//
// 2. STAR-CONVEXITY OF TARGET CELL: The target cell can be non-convex
// but it MUST BE STAR-CONVEX, i.e. allow a decomposition into
// positive volume simplices. Since each simplex is convex, one can
// clip the source cell with the planes of each simplex to get the
// intersection volume between the simplex and the source cell. The
// intersection volumes from all the simplices of the decomposition
// can be summed up to get the total intersection volume.
//
//
// Summary: We will use the star-convexity check for both the source
// and the target polyhedra as it's the more stringent check. For
// cells or dual cells, this decomposition is already available in the
// form of wedges. Material polygons/polyhedra from interface
// reconstruction, on the other hand, are tested using the divergence
// theorem since they already undergo strict checks in the interface
// reconstruction procedure. Also, in 3D they are guaranteed to be
// convex since a non-convex cell is decomposed into simplices and
// then sliced by the interface plane.



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

#ifdef DEBUG
    if (rectangular_mesh_) {  // may have a special decomposition into 5/6 tets
      int ntets = target_tet_coords.size();
      for (int i = 0; i < ntets; i++) {
        // Let's check the volume of the target
        const Point<3> p0 = target_tet_coords[i][0];
        const Point<3> p1 = target_tet_coords[i][1];
        const Point<3> p2 = target_tet_coords[i][2];
        const Point<3> p3 = target_tet_coords[i][3];
        
        Vector<3> v0 = p1-p0;
        Vector<3> v1 = p2-p0;
        Vector<3> cpvec = cross(v0,v1);
        Vector<3> v2 = p3-p0;
        double wvol = dot(v2, cpvec)/6.0;
        if (wvol < 0.0) {
          // if a tet in target cell decomposition is negative, throw an error
          std::stringstream sstr;
          sstr << "In intersect_r3d.h: " <<
              " Tet in decomposition of target cell " << tgt_cell <<
              " has negative volume " << wvol << "\n";
          double cvol = targetMeshWrapper.cell_volume(tgt_cell);
          if (cvol <= 0.0)
            sstr << "The cell (vol = " << cvol << ") is inside out or degenerate";
          else
            sstr << "The cell (vol = " << cvol << ") has positive volume";
          throw std::runtime_error(sstr.str());
        }
      }
    } else {  // just check the sides of the cell (AuxMeshTopology always enables sides)
      if (targetMeshWrapper.num_entities(SIDE, ALL) == 0) {
        std::stringstream sstr;
        sstr << "In intersect_r3d:\n" <<
            " Decomposition of cells into sides not available." <<
            " Cannot check validity of input.\n" <<
            " Request wedges in mesh wrapper creation or make sure your " <<
            " mesh framework supports cell_get_sides and side_volume methods";
        throw std::runtime_error(sstr.str());
      }
    
      std::vector<int> sides;
      targetMeshWrapper.cell_get_sides(tgt_cell, &sides);

      for (auto const& sd : sides) {
        double svol = targetMeshWrapper.side_volume(sd);
        if (svol < 0.0) {
          std::stringstream sstr;
          sstr << "In intersect_r3d.h: " <<
              " Tet in decomposition of target cell " << tgt_cell <<
              " has negative volume " << svol << "\n";
          double cvol = targetMeshWrapper.cell_volume(tgt_cell);
          if (cvol <= 0.0)
            sstr << "The cell (vol = " << cvol << ") is inside out or degenerate";
          else
            sstr << "The cell (vol = " << cvol << ") may be highly non-convex";
          throw std::runtime_error(sstr.str());
        }
      }
    }
#endif

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

#ifdef DEBUG
        // Check validity of source cell (unfortunately may be
        // repeated if it is an intersection candidate for multiple
        // target cells)
        
        if (sourceMeshWrapper.num_entities(SIDE, ALL) == 0) {
          std::stringstream sstr;
          sstr << "In intersect_r3d:" <<
              " Decomposition of cells into sides not available." <<
              " Cannot check validity of input.\n" <<
              " Request sides in mesh wrapper creation or make sure your " <<
              " mesh framework supports cell_get_sides and side_volume methods";
          throw std::runtime_error(sstr.str());
        }
        
        std::vector<int> sides;
        sourceMeshWrapper.cell_get_sides(s, &sides);

        for (auto const& sd: sides) {
          double svol = sourceMeshWrapper.side_volume(sd);
          if (svol < 0.0) {
            std::stringstream sstr;
            sstr << "In intersect_polys_r3d.h: " <<
                "Tet in decomposition of source cell " << s <<
                " has negative volume " << svol << "\n ";
            double cvol = sourceMeshWrapper.cell_volume(s);
            if (cvol <= 0.0)
              sstr << "The cell (vol = " << cvol << ") is inside out or degenerate";
            else
              sstr << "The cell (vol = " << cvol << ") may be highly non-convex";
            throw std::runtime_error(sstr.str());            
          }
        }        
#endif

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

#ifdef DEBUG
          // Lets check the volume of the source material polyhedron

          std::vector<double> smom = matpoly.moments();
          if (smom[0] < 0.0) {
            std::stringstream sstr;
            sstr << "In intersect_polys_r3d.h: " <<
                "Material polygon for material " << matid_ << " in cell " <<
                s << " has negative volume " << smom[0] << "\n";
            throw std::runtime_error(sstr.str());
          }
#endif
          
          facetedpoly_t srcpoly = get_faceted_matpoly(matpoly);
          std::vector<double> momvec = intersect_polys_r3d(srcpoly,
                                          target_tet_coords, num_tols_);
          for (int k = 0; k < 4; k++)
            this_wt.weights[k] += momvec[k];
        }

      }
#else  // No Tangram
      facetedpoly_t srcpoly;
      sourceMeshWrapper.cell_get_facetization(s, &srcpoly.facetpoints,
                                              &srcpoly.points);
      
#ifdef DEBUG
      // Check validity of source cell (unfortunately may be
      // repeated if it is an intersection candidate for multiple
      // target cells)
      
      if (sourceMeshWrapper.num_entities(SIDE, ALL) == 0) {
        std::stringstream sstr;
        sstr << "In intersect_r3d:" <<
            " Decomposition of cells into sides not available." <<
            " Cannot check validity of input.\n" <<
            " Request sides in mesh wrapper creation or make sure your " <<
            " mesh framework supports cell_get_sides and side_volume methods";
        throw std::runtime_error(sstr.str());
      }
    
      std::vector<int> sides;
      sourceMeshWrapper.cell_get_sides(s, &sides);
      
      for (auto const& sd: sides) {
        double svol = sourceMeshWrapper.side_volume(sd);
        if (svol < 0.0) {
          std::stringstream sstr;
          sstr << "In intersect_polys_r3d.h: " <<
              "Tet in decomposition of SOURCE cell " << s <<
              " has negative volume " << svol << "\n ";
          double cvol = sourceMeshWrapper.cell_volume(s);
          if (cvol <= 0.0)
            sstr << "The cell (vol = " << cvol << ") is inside out or degenerate";
          else
            sstr << "The cell (vol = " << cvol << ") may be highly non-convex";
          throw std::runtime_error(sstr.str());            
        }
      }        
#endif

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

    // We should avoid any decomposition for dual cells of a
    // rectangular mesh but for now we will decompose the target all
    // the time

    targetMeshWrapper.dual_wedges_get_coordinates(tgt_node, &target_tet_coords);

#ifdef DEBUG
    if (targetMeshWrapper.num_entities(WEDGE, ALL) == 0) {
      std::stringstream sstr;
      sstr << "In intersect_r3d:" <<
          " Wedge decomposition of cells not requested." <<
          " Cannot check validity of input.\n" <<
          " Request wedges in mesh wrapper creation or make sure your " <<
          " mesh framework supports node_get_wedges and wedge_volume methods";
      throw std::runtime_error(sstr.str());
    }
    
    // Lets check if all the wedges in the dual cell are valid
    std::vector<int> wedges;
    targetMeshWrapper.node_get_wedges(tgt_node, Entity_type::ALL, &wedges);

    for (auto const& w : wedges) {
      double wvol = targetMeshWrapper.wedge_volume(w);
      if (wvol < 0.0) {
        std::stringstream sstr;
        sstr << "In intersect_r3d.h: " <<
            " Tet in decomposition of dual cell of target node " << tgt_node <<
            " has negative volume " << wvol << "\n";
        double dvol = targetMeshWrapper.dual_cell_volume(tgt_node);
        if (dvol <= 0.0)
          sstr << "The cell (vol = " << dvol << ") is inside out or degenerate";
        else
          sstr << "The cell (vol = " << dvol << ") may be highly non-convex";
        throw std::runtime_error(sstr.str());
      }
    }
#endif
    

    // CAN MAKE THIS INTO A THRUST TRANSFORM CALL
    int nsrc = src_nodes.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_nodes[i];

      facetedpoly_t srcpoly;
      sourceMeshWrapper.dual_cell_get_facetization(s, &srcpoly.facetpoints,
                                                   &srcpoly.points);

      
#ifdef DEBUG
      // Lets check that the wedges in the source dual cell have
      // positive volume

      if (sourceMeshWrapper.num_entities(WEDGE, ALL) == 0) {
        std::stringstream sstr;
        sstr << "In intersect_r3d:" <<
            " Wedge decomposition of cells not requested." <<
            " Cannot check validity of input.\n" <<
            " Request wedges in mesh wrapper creation or make sure your " <<
            " mesh framework supports node_get_wedges and wedge_volume methods";
        throw std::runtime_error(sstr.str());
      }
      
      std::vector<int> wedges;
      sourceMeshWrapper.node_get_wedges(s, Entity_type::ALL, &wedges);

      for (auto const& w: wedges) {
        double wvol = sourceMeshWrapper.wedge_volume(w);
        if (wvol < 0.0) {
          std::stringstream sstr;
          sstr << "In intersect_polys_r3d.h: " <<
              "Tet in decomposition of dual cell of source node " << s <<
              " has negative volume " << wvol << "\n";
          double dvol = targetMeshWrapper.dual_cell_volume(s);
          if (dvol <= 0.0)
            sstr << "The cell (vol = " << dvol << ") is inside out or degenerate";
          else
            sstr << "The cell (vol = " << dvol << ") may be highly non-convex";
          throw std::runtime_error(sstr.str());
        }
      }
#endif

      
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
