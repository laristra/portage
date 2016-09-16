/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_MESH_WRAPPER_H_
#define JALI_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <utility>

#include "Mesh.hh"                      // Jali mesh header

#include "portage/wrappers/mesh/AuxMeshTopology.h"
#include "portage/support/portage.h"
#include "portage/support/Point.h"


namespace Portage {

/*!
  \class Jali_Mesh_Wrapper jali_mesh_wrapper.h
  \brief Jali_Mesh_Wrapper implements mesh methods for Jali

  Jali_Mesh_Wrapper implements methods required for Portage mesh
  queries for the Jali mesh infrastructure. It is derived from a
  helper class, called AuxMeshTopology that provides the
  side/wedge/corner functionality. The helper class should be
  templated on this mesh wrapper itself because it relies on the mesh
  wrapper to answer some queries about the basic topology through a
  specific interface.  This uses a template pattern called Curiously
  Recurring Template Pattern (CRTP) which allows this kind of mutual
  invocation. If the mesh wrapper has the auxiliary topology already,
  it can choose not to be derived from the helper class. See
  https://en.m.wikipedia.org/wiki/Curiously_recurring_template_pattern
*/

//! helper function:  convert Jali point to Portage point
template <long D>
Point<D> toPortagePoint(const JaliGeometry::Point& jp) {
  Point<D> pp;
  assert(jp.dim() == D);
  for (int d = 0; d < D; ++d)
    pp[d] = jp[d];
  return pp;
}


class Jali_Mesh_Wrapper : public AuxMeshTopology<Jali_Mesh_Wrapper> {
 public:

  //! Constructor
  //!
  //! It is possible to construct this class with all, some or none of
  //! the auxiliary entities requested (usually to save memory)

  explicit Jali_Mesh_Wrapper(Jali::Mesh const & mesh,
                             bool request_sides = true,
                             bool request_wedges = true,
                             bool request_corners = true) :
      jali_mesh_(mesh),
      AuxMeshTopology<Jali_Mesh_Wrapper>(request_sides, request_wedges,
                                         request_corners)
  {
    build_aux_entities();  // base class (AuxMeshTopology) method that has
    //                     // to be called here and not in the constructor
    //                     // of the base class because it needs access to
    //                     // methods in this class which in turn need access
    //                     // to its member variables. But these member vars
    //                     // don't get initialized until the base class is
    //                     // constructed
  }

  //! Copy constructor
  Jali_Mesh_Wrapper(Jali_Mesh_Wrapper const & inmesh) :
      jali_mesh_(inmesh.jali_mesh_)
  {}

  //! Assignment operator (Deleted)
  Jali_Mesh_Wrapper & operator=(Jali_Mesh_Wrapper const & inmesh) = delete;

  //! Empty destructor
  ~Jali_Mesh_Wrapper() {}

  //! Dimension of space or mesh points
  int space_dimension() const {
    return jali_mesh_.space_dimension();
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED);
  }

  //! Number of owned faces in the mesh
  int num_owned_faces() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::FACE,
                                   Jali::Entity_type::PARALLEL_OWNED);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED);
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_GHOST);
  }

  //! Number of ghost faces in the mesh
  int num_ghost_faces() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::FACE,
                                   Jali::Entity_type::PARALLEL_GHOST);
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_GHOST);
  }

  //! Get the type of the cell - PARALLEL_OWNED or PARALLEL_GHOST
  //! Assumes a 1-1 correspondence between integer values of the
  //! enum types to avoid switch statements

  Portage::Entity_type cell_get_type(int const cellid) const {
    static Portage::Entity_type jali2portage_type[5] =
        {DELETED, PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL};
    Jali::Entity_type etype =
        jali_mesh_.entity_get_type(Jali::Entity_kind::CELL, cellid);
    return jali2portage_type[static_cast<int>(etype)];
  }

  //! Get the element type of a cell - TRI, QUAD, POLYGON, TET, HEX,
  //! PRISM OR POLYHEDRON

  Portage::Element_type cell_get_element_type(int const cellid) const {
    static Portage::Element_type jali2portage_elemtype[9] =
        {UNKNOWN_TOPOLOGY, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX,
         POLYHEDRON};
    Jali::Cell_type ctype = jali_mesh_.cell_get_type(cellid);
    return jali2portage_elemtype[static_cast<int>(ctype)];
  }

  //! Get cell faces and the directions in which they are used
  void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const {
    std::vector<Jali::dir_t> fdirs;

    jali_mesh_.cell_get_faces_and_dirs(cellid, cfaces, &fdirs);

    cfdirs->resize(fdirs.size());
    for (int i = 0; i < fdirs.size(); ++i)
      (*cfdirs)[i] = static_cast<int>(fdirs[i]);
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int cellid, std::vector<int> *nodes) const {
    jali_mesh_.cell_get_nodes(cellid, nodes);
  }


  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {
    jali_mesh_.cell_get_node_adj_cells(cellid, (Jali::Entity_type) ptype,
                                       adjcells);
  }

  //! \brief Get "adjacent" nodes of given node
  //!
  //! Get "adjacent" nodes of given node - nodes that share a common
  //! cell with given node
  void node_get_cell_adj_nodes(int const nodeid,
                               Entity_type const ptype,
                               std::vector<int> *adjnodes) const {
    adjnodes->clear();

    Jali::Entity_ID_List nodecells;
    jali_mesh_.node_get_cells(nodeid, (Jali::Entity_type) ptype, &nodecells);

    for (auto const& c : nodecells) {
      Jali::Entity_ID_List cellnodes;
      jali_mesh_.cell_get_nodes(c, &cellnodes);

      for (auto const& n : cellnodes) {
        if (n == nodeid) continue;
        if (std::find(adjnodes->begin(), adjnodes->end(), n) == adjnodes->end())
          adjnodes->emplace_back(n);
      }
    }
  }

  //! Get nodes of a face
  void face_get_nodes(int const faceid, std::vector<int> *fnodes) const {
    jali_mesh_.face_get_nodes(faceid, fnodes);
  }

  //! Get global id
  int get_global_id(int const id, Entity_kind const kind) const {
    return jali_mesh_.GID(id, (Jali::Entity_kind)kind);
  }

  //! Virtual and local addresses are equivalent in non-distributed case
  int virtual_to_local(int virtualId) const { return virtualId; }

  //! coords of a node
  template <long D>
  void node_get_coordinates(int const nodeid, Point<D>* pp) const {
    JaliGeometry::Point jp;
    jali_mesh_.node_get_coordinates(nodeid, &jp);
    assert(jp.dim() == D);
    *pp = toPortagePoint<D>(jp);
  }

  //! coords of nodes of a cell
  template <long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Point<D>> *pplist) const {
    assert(jali_mesh_.space_dimension() == D);

    std::vector<JaliGeometry::Point> jplist;
    jali_mesh_.cell_get_coordinates(cellid, &jplist);

    pplist->resize(jplist.size());
    // This cast appears necessary for proper template deduction
    std::transform(jplist.begin(), jplist.end(), pplist->begin(),
                   (Point<D>(*)(const JaliGeometry::Point&))toPortagePoint<D>);
  }

  /// \brief Centroid of a cell
  //
  // Return the centroid of a cell

  template <long D>
  void cell_centroid(Jali::Entity_ID cellid, Point<D> *ccentroid) const {
    *ccentroid = toPortagePoint<D>(jali_mesh_.cell_centroid(cellid));
  }

  /// \brief Centroid of a face
  //
  // Return the centroid of a cell face

  template <long D>
  void face_centroid(Jali::Entity_ID faceid, Point<D> *fcentroid) const {
    *fcentroid = toPortagePoint<D>(jali_mesh_.face_centroid(faceid));
  }

  //! Cell area/volume
  double cell_volume(int cellID) const {
    return jali_mesh_.cell_volume(cellID);
  }

 private:
  Jali::Mesh const & jali_mesh_;

};  // class Jali_Mesh_Wrapper


}  // end namespace Portage

#endif // JALI_MESH_WRAPPER_H_
