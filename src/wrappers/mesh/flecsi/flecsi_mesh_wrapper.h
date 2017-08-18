/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef FLECSI_MESH_WRAPPER_H_
#define FLECSI_MESH_WRAPPER_H_

namespace std
{
    typedef decltype(nullptr) nullptr_t;
}


#include <cassert>
#include <algorithm>
#include <array>

#include "flecsi/specializations/burton/burton.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"

using mesh_t = flecsi::burton_mesh_t;

using real_t = flecsi::mesh_t::real_t;
using vertex_t = flecsi::mesh_t::vertex_t;
using point_t = flecsi::mesh_t::point_t;

namespace Portage {

//! helper function to convert to Portage::Point
template <long D>
Portage::Point<D> toPortagePoint(const point_t &fp) {
  assert(fp.size() == D);
  Portage::Point<D> pp;
  for (auto i = 0; i < D; ++i)
    pp[i] = fp[i];
  return pp;
}

/*!
  @brief Helper routine to make a cartesian 2d grid in FleCSI
  @param[in] xmin,xmax,ymin,ymax The extents of the mesh
  @param[in] ncellsx,ncellsy The number of cells in each direction
  @param[in,out] mesh The flecsi::burton_mesh_t mesh object
*/
void make_mesh_cart2d(const real_t xmin, const real_t xmax,
                      const real_t ymin, const real_t ymax,
                      const int ncellsx, const int ncellsy,
                      mesh_t & mesh) {
  // grid spacing
  auto dx = (xmax - xmin) / real_t(ncellsx);
  auto dy = (ymax - ymin) / real_t(ncellsy);
  
  //  mesh_t mesh;
  auto num_verts = (ncellsx+1)*(ncellsy+1);
  // this initializes storage for the number of vertices
  mesh.init_parameters(num_verts);
  
  // create the verts
  std::vector<vertex_t*> verts;
  for (auto j = 0; j < ncellsy+1; ++j) {
    for (auto i = 0; i < ncellsx+1; ++i) {
      auto vert = mesh.create_vertex(
          {xmin + dx*real_t(i), ymin + dy*real_t(j)});
      verts.push_back(vert);
    }
  }
  
  // create the cells
  auto ncellsx1 = ncellsx+1;
  for (auto j = 0; j < ncellsy; ++j) {
    for (auto i = 0; i < ncellsx; ++i) {
      // go over verts in counter-clockwise fashion
      auto c = mesh.create_cell({verts[i + j*ncellsx1],
              verts[i + 1 + j*ncellsx1],
              verts[i + 1 + (j + 1)*ncellsx1],
              verts[i + (j + 1)*ncellsx1]});
    }
  }
  
  // this setups up the faces, edges, wedges, etc and connectivity info
  mesh.init();
  
}

/*!
  @class Flecsi_Mesh_Wrapper flecsi_mesh_wrapper.h
  @brief Flecsi_Mesh_Wrapper implements mesh methods for Flecsi

  Flecsi_Mesh_Wrapper implements methods required for Portage mesh
  queries for the Flecsi mesh infrastructure
*/

class Flecsi_Mesh_Wrapper {
 public:

  //! Constructor
  explicit Flecsi_Mesh_Wrapper(mesh_t & mesh) :
      flecsi_mesh_(mesh)
  {}

  //! Copy constructor
  Flecsi_Mesh_Wrapper(Flecsi_Mesh_Wrapper const & inmesh) :
      flecsi_mesh_(inmesh.flecsi_mesh_)
  {}

  //! Assignment operator (disabled)
  Flecsi_Mesh_Wrapper & operator=(Flecsi_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flecsi_Mesh_Wrapper() {}


  //! Dimension of space or mesh points
  int space_dimension() const {
    return flecsi_mesh_.dimension();
  }

  //! Cell area/volume
  double cell_volume(int const cellID) const {
    if (space_dimension() > 2)
      assert(false && "FleCSI 3D not implemented");
    return flecsi_mesh_.cells()[cellID]->area();
  }

  //! Dual cell area/volume
  double dual_cell_volume(int const nodeid) const {
    if (space_dimension() > 2)
      assert(false && "FleCSI 3D not implemented");
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    double vol = 0.0;
    for (auto corner : flecsi_mesh_.corners(thisNode))
      vol += corner->area();
    return vol;
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return flecsi_mesh_.num_cells();
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return flecsi_mesh_.num_vertices();
  }

  //! Number of owned faces in the mesh
  int num_owned_faces() const {
    return 0;  // Only 2d is implemented and no faces available in 2D FleCSI
  }

  // NOTE: I don't know where to get the ghosts yet, so setting this to 0.
  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return 0;
  }

  // NOTE: I don't know where to get the ghosts yet, so setting this to 0.
  //! Number of ghost faces in the mesh
  int num_ghost_faces() const {
    return 0;
  }

  // NOTE: I don't know where to get the ghosts yet, so setting this to 0.
  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return 0;
  }

  //! Number of items of given entity
  int num_entities(Entity_kind const entity, Entity_type const etype=Entity_type::ALL) const {
    switch(entity) {
      case NODE :
        return flecsi_mesh_.num_vertices();
        break;
      case EDGE :
        return flecsi_mesh_.num_edges();
        break;
        // This case does not work in FleCSI because only 2d is implemented and
        // a 2d face is an edge
      /* case FACE : */
      /*   return flecsi_mesh_.num_faces(); */
      /*   break; */
      case CELL :
        return flecsi_mesh_.num_cells();
        break;
      case WEDGE :
        return flecsi_mesh_.num_wedges();
        break;
      case CORNER :
        return flecsi_mesh_.num_corners();
        break;
      default :
        assert(false && "Error: invalid entity request in num_entities");
    }
  }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity,
                          Entity_type const etype = Entity_type::ALL) const {
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity,
                        Entity_type const etype=Entity_type::ALL) const {
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_entities(entity, etype));
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int cellid, std::vector<int> *nodes) const {
    auto thisCell = flecsi_mesh_.cells()[cellid];
    nodes->clear();
    for (auto v : flecsi_mesh_.vertices(thisCell))
      nodes->emplace_back(v.global_id());
  }

  //! Get cell faces and the directions in which they are used
  void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const {

    // Do nothing - faces not represented in 2D flecsi
    
    cfaces->clear();
    cfdirs->clear();

  }

  //! Get nodes of a face
  void face_get_nodes(int const faceid, std::vector<int> *fnodes) const {
    // Do nothing - faces not represented in 2D flecsi

    fnodes->clear();
  }

  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {
    // I think this is supposed to work, but doesn't:
    // flecsi_mesh_.cells(mycell);
    auto thisCell = flecsi_mesh_.cells()[cellid];
    adjcells->clear();
    // Loop over all nodes of this cell
    for (auto node : flecsi_mesh_.vertices(thisCell)) {
      // Loop over all cells associated with this node
      for (auto cell : flecsi_mesh_.cells(node)) {
        if (cell != thisCell)
          adjcells->emplace_back(cell.global_id());
      }
    }
  }

  //! @brief Get "adjacent" nodes of given node
  //!
  //! Get "adjacent" nodes of given node - nodes that share a common
  //! cell with given node
  void node_get_cell_adj_nodes(int const nodeid,
                               Entity_type const ptype,
                               std::vector<int> *adjnodes) const {
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    adjnodes->clear();
    // Loop over cells associated with this node
    for (auto cell : flecsi_mesh_.cells(thisNode)) {
      // Loop over nodes of this cell
      for (auto node : flecsi_mesh_.vertices(cell)) {
        if (thisNode != node)
          adjnodes->emplace_back(node.global_id());
      }
    }
  }

  //! @brief Get adjacent "dual cells" of a given "dual cell"
  void dual_cell_get_node_adj_cells(int const nodeid,
                                    Entity_type const ptype,
                                    std::vector<int> *adjnodes) const {
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    adjnodes->clear();
    // Loop over cells associated with this node
    for (auto cell : flecsi_mesh_.cells(thisNode)) {
      // Loop over the nodes associated with this cell
      for (auto node : flecsi_mesh_.vertices(cell)) {
        if (thisNode != node)
          adjnodes->emplace_back(node.global_id());
      }
    }
  }

  /*!
    @brief Templated version of coords of a node
    @tparam D The dimension of the mesh.
    @param[in] nodeid The ID of the node.
    @param[in,out] pp The Portage::Point object containing the coordinate
    information.

    Note: FleCSI Burton specialization doesn't currently fully support 1 or 3D.
  */
  template <long D>
  void node_get_coordinates(int const nodeid, Portage::Point<D>* pp) const {
    auto thisVert = flecsi_mesh_.vertices()[nodeid];
    *pp = toPortagePoint<D>(thisVert->coordinates());
  }


  /*!
    @brief Templated version of coodinates of the nodes of a cell.
    @tparam D The dimension of the mesh.
    @param[in] cellid The ID of the cell.
    @param[in,out] pplist The vector of Portage::Point objects containing the
    coordinates of a node.  The length of the vector is equal to the number of
    nodes in the cell with ID @c cellid.

    Note: FleCSI Burton specialization doesn't currently fully support 1 or 3D.
  */
  template <long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Portage::Point<D>> *pplist)  const {
    // Get this cell object
    auto thisCell = flecsi_mesh_.cells()[cellid];

    // Loop over vertices of this cell to get their coordinates
    pplist->clear();
    auto theseVerts = flecsi_mesh_.vertices(thisCell);
    for (auto v : theseVerts)
      pplist->emplace_back(toPortagePoint<D>(v->coordinates()));
  }


  /*!
    @brief 2D version of coords of nodes of a dual cell
    @param[in] nodeid The ID of the node or dual cell in the dual mesh.
    @param[in,out] pplist The vector of Portage::Point objects containing the
    coordinates of a node in the dual mesh / cell in the regular mesh.  The
    length of the vector is equal to the number of nodes in the dual mesh cell
    with ID @c nodeid.

    The vertices are ordered CCW. For node @c nodeid not on a
    boundary, the vector @c pplist starts with a random vertex, but it is still
    ordered CCW. Use the dual_cell_coordinates_canonical_rotation() function to
    rotate the @c pplist into a canonical (unique) form.

    @TODO worry about boundary cases
  */
  void dual_cell_get_coordinates(int const nodeid,
                                 std::vector<Portage::Point2 > *pplist) const {
    assert(space_dimension() == 2);
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    pplist->clear();
    // Loop over the corners associated with this node
    for (auto corner : flecsi_mesh_.corners(thisNode)) {
      std::vector<int> wedgeIds;
      // The two wedges in this corner
      for (auto wedge : flecsi_mesh_.wedges(corner))
        wedgeIds.emplace_back(wedge.global_id());
      order_wedges_ccw(&wedgeIds);
      // Get the coordinates of each wedgeId, given the new ordering
      for (auto wId : wedgeIds) {
        auto w = flecsi_mesh_.wedges()[wId];
        // push back centroid and then edge midpoint
        auto cc = w->cell()->centroid();
        pplist->emplace_back(toPortagePoint<2>(cc));//cc[0], cc[1]);
        auto nc = w->edge()->midpoint();
        pplist->emplace_back(toPortagePoint<2>(nc));//nc[0], nc[1]);
      }
    }
  }

  void order_wedges_ccw(std::vector<int> *wedgeids) const {
    assert(wedgeids->size() == 2);
    auto firstWedge = flecsi_mesh_.wedges()[(*wedgeids)[0]];
    // NOTE: This mimics the order of Jali in 2d: node, face/edge, cell
    std::vector<Portage::Point2> wcoords;
    auto nc = firstWedge->vertex()->coordinates();
    wcoords.emplace_back(nc[0], nc[1]);
    auto ec = firstWedge->edge()->midpoint();
    wcoords.emplace_back(ec[0], ec[1]);
    auto cc = firstWedge->cell()->centroid();
    wcoords.emplace_back(cc[0], cc[1]);

    // Ensure (*wedgeids)[0] is the first wedge in a CCW direction
    if (not ccw(wcoords[0], wcoords[1], wcoords[2]))
      std::swap((*wedgeids)[0], (*wedgeids)[1]);
  }

  // Returns true if the three 2D points (p1, p2, p3) are a counter-clockwise
  // turn, otherwise returns false (corresponding to clockwise or collinear)
  bool ccw(const Portage::Point2 p1,
           const Portage::Point2 p2,
           const Portage::Point2 p3) const {
    return (p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]) > 0;
  }

  std::vector<Portage::Point2> cellToXY(int const cellID) const {
    std::vector<Portage::Point2> cellPoints;
    cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }

  // NOTE: This ASSUMES 3D - the "3" is for the 3 coordinates of a spatial point
  // and the "4" is for the four points that make up a wedge (tet) in 3d
  // NOTE: FleCSI doesn't have 3D yet!!!
  void wedges_get_coordinates(int const cellID,
                              std::vector<std::array<Portage::Point3, 4>>
                              *wcoords) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented");
  }

  // Get the simplest possible decomposition of a 3D cell into tets.
  // NOTE: FleCSI doesn't have 3D yet!!!
  void decompose_cell_into_tets(int const cellID,
                                std::vector<std::array<Portage::Point3, 4>>
                                *tcoords, const bool planar_hex) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented");
  }

  // NOTE: FleCSI doesn't have 3D yet!!!
  //! 3D version of coords of nodes of a dual cell
  // Input is the node ID 'nodeid', and it returns the vertex coordinates of
  // the dual cell around this node in `pplist`.  The vertices are NOT ordered
  // in any particular way
  void dual_cell_get_coordinates(int const nodeid,
        			 std::vector<Portage::Point3> *pplist) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented");
  }

  // NOTE: This ASSUMES 3D - the "3" is for the 3 coordinates of a spatial point
  // and the "4" is for the four points that make up a wedge (tet) in 3d
  // NOTE: FleCSI doesn't have 3D yet!!!
  // Get the coordinates of the wedges of the dual mesh
  void dual_wedges_get_coordinates(int const nodeID,
                                   std::vector<std::array<Portage::Point3, 4>>
                                   *wcoords) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented yet");
  }


  /*!
    @brief Centroid of a cell.
    @param[in] cellid The ID of the cell.
    @param[in,out] centroid The vector of coordinates of the cell @c cellid's
    centroid.  The length of the vector is equal to the dimension of the mesh.
  */

  template<long D>
  void cell_centroid(int const cellid,
                     Point<D> *centroid) const {
    auto thisCell = flecsi_mesh_.cells()[cellid];
    *centroid = toPortagePoint<D>(thisCell->centroid());
  }

  //! Get global id
  int get_global_id(int const id, Entity_kind const kind) const {
    return id;
  }

  /*!
    @brief Centroid of a dual cell.
    @param[in] nodeid The ID of the node in the normal mesh / cell in the dual
    mesh.
    @param[in,out] centroid The vector of coordinates of the node in the normal
    mesh / the cell in the dual mesh with ID @c nodeid.  The length of the
    vector is equal to the dimension of the mesh.

    @TODO NOTE: THIS IS ASSUMED TO BE THE NODE COORDINATE BECAUSE
    THE NODAL VARIABLES LIVE THERE, BUT FOR DISTORTED GRIDS, THE
    NODE COORDINATED MAY NOT BE THE CENTROID OF THE DUAL CELL
   */
  template<long D>
  void dual_cell_centroid(int const nodeid,
                          Point<D> *centroid) const {
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    *centroid = toPortagePoint<D>(thisNode->coordinates());
  }

 private:
  mesh_t & flecsi_mesh_;

};  // class Flecsi_Mesh_Wrapper

}  // end namespace Portage

#endif  // FLECSI_MESH_WRAPPER_H_
