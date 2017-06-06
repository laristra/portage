/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsale/geom/shapes/geometric_shapes.h"
#include "flecsale/utils/errors.h"

// library includes
#include <portage/support/portage.h>
#include <portage/support/Point.h>
#include <portage/wonton/mesh/AuxMeshTopology.h>

// system includes
#include <map>
#include <memory>
#include <vector>

namespace flecsale {
namespace mesh {

namespace portage {

////////////////////////////////////////////////////////////////////////////////
/// \brief A utility function to convert a type to a portage point
////////////////////////////////////////////////////////////////////////////////
// @{
template<
  typename T, 
  template<typename,std::size_t> class A
>
Portage::Point<3> make_point( const A<T,3> & a )
{
  return { a[0], a[1], a[2] };
}

template<
  typename T, 
  template<typename,std::size_t> class A
>
Portage::Point<2> make_point( const A<T,2> & a )
{
  return { a[0], a[1] };
}

template<typename T>
Portage::Point<1> make_1d_point( T && a )
{
  return { std::forward<T>(a)[0] };
}

template<typename T>
Portage::Point<2> make_2d_point( T && a )
{
  return { std::forward<T>(a)[0], std::forward<T>(a)[1] };
}

template<typename T>
Portage::Point<3> make_3d_point( T && a )
{
  return { 
    std::forward<T>(a)[0], std::forward<T>(a)[1], std::forward<T>(a)[2] 
  };
}
// @}


} // namespace 

////////////////////////////////////////////////////////////////////////////////
/// \brief The base class for the portage flecsi mesh wrapper.
///
/// This is used to hide some static initialization for this header-only
/// implementation.
////////////////////////////////////////////////////////////////////////////////

template< typename M >
struct portage_mesh_base_t {

  //! \brief The mesh type
  using mesh_t = M;

  //! the gometric shape 
  using shape_t = typename mesh_t::shape_t;
  //! \brief the portage element type
  using element_type_t = Portage::Element_type;

  //! \brief a map type for mapping mesh shapes to portage shapes
  using shape_map_t = std::map< shape_t, element_type_t >;

  //! \brief a mapping between mesh shapes and portage
  static const shape_map_t shapes_to_portage;
};

// Static initialization of portage_mesh_base_t::shapes_to_portage
template< typename M >
const typename portage_mesh_base_t<M>::shape_map_t 
  portage_mesh_base_t<M>::shapes_to_portage = 
{ 
  { 
    shape_t::triangle, 
    element_type_t::TRI 
  },
  { 
    shape_t::quadrilateral, 
    element_type_t::QUAD 
  },
  { 
    shape_t::polygon, 
    element_type_t::POLYGON 
  },
  { 
    shape_t::tetrahedron, 
    element_type_t::TET
  },
  { 
    shape_t::hexahedron, 
    element_type_t::HEX 
  },
  { 
    shape_t::polyhedron, 
    element_type_t::POLYHEDRON 
  }
};


////////////////////////////////////////////////////////////////////////////////
///  \brief Implements a mesh wrapper for Portage mesh queries.
////////////////////////////////////////////////////////////////////////////////
template< typename M >
class portage_mesh_t : public Portage::AuxMeshTopology<portage_mesh_t<M>> 
{

public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the auxiliary class
  using portage_mesh_aux_t = Portage::AuxMeshTopology<portage_mesh_t<M>>;

  //! \brief The mesh type
  using mesh_t = M;
  //! \brief the size type
  using size_t = typename mesh_t::size_t;
  //! \brief the real type
  using real_t = typename mesh_t::real_t;
  //! the gometric shape 
  using shape_t = typename mesh_t::shape_t;

  //! \brief The entity kind type
  using entity_kind_t = Portage::Entity_kind;
  //! \brief The entity type 
  using entity_type_t = Portage::Entity_type;
  //! \brief the portage point type
  using point_t = Portage::Point< mesh_t::num_dimensions >;
  //! \brief the portage element type
  using element_type_t = Portage::Element_type;

  //! \brief a map type for mapping mesh shapes to portage shapes
  using shape_map_t = std::map< shape_t, element_type_t >;

  //! \brief the 2d portage point type
  using point_1d_t = Portage::Point< 1 >;
  //! \brief the 2d portage point type
  using point_2d_t = Portage::Point< 2 >;
  //! \brief the 3d portage point type
  using point_3d_t = Portage::Point< 3 >;

  //============================================================================
  // Public Static Data
  //============================================================================

  //! \brief the map between 
  static const shape_map_t shapes_to_portage;

  //============================================================================
  // Constructors
  //============================================================================

  //!  \brief Constructor for creating a serial, 3D Cartesian mesh.
  //!  \param[in] mesh The minimum coordinates of the domain.
  explicit portage_mesh_t(mesh_t & mesh) :
    cells_(mesh.cells()), faces_(mesh.faces()),
    vertices_(mesh.vertices()), mesh_(&mesh)
  {
    // base class (AuxMeshTopology) method that has to be called here
    // and not in the constructor of the base class because it needs
    // access to methods in this class which in turn need access to
    // its member variables. But these member vars don't get
    // initialized until the base class is constructed
    if ( mesh_t::num_dimensions == 3 )
      portage_mesh_aux_t::build_aux_entities(); 
  }

  //! Default constructor deleted
  portage_mesh_t() = default;

  //! Default copy constructor
  portage_mesh_t(const portage_mesh_t &) = default;

  //! Default assignment operator
  portage_mesh_t & operator=(const portage_mesh_t &) = default;

  //============================================================================
  // Public Members Required By Portage
  //============================================================================

  //! Dimension of space or mesh points
  constexpr auto space_dimension() const 
  {
    return mesh_t::num_dimensions;
  }

  //! Cell area/volume
  auto cell_volume(size_t cell_id) const 
  {
    return cells_[cell_id]->volume();
  }

  //! Dual cell area/volume
  //! \param [in] node_id the node index
  auto dual_cell_volume(size_t node_id) const 
  {
    //auto v = mesh_->vertices()[node_id];
    //real_t vol = 0.0;
    //for (auto corner : mesh_->corners(v))
    //  vol += corner->area();
    //return vol;

    raise_implemented_error( "dual_cell_volume not implemented yet!" );
    return 0.0;
  }

  //! Number of owned cells in the mesh
  size_t num_owned_cells() const 
  {
    return mesh_->num_cells();
  }

  //! Number of owned faces in the mesh
  size_t num_owned_faces() const 
  {
    return mesh_->num_faces();
  }

  //! Number of owned edges in the mesh
  size_t num_owned_edges() const 
  {
    return mesh_->num_edges();
  }

  //! Number of owned nodes in the mesh
  size_t num_owned_nodes() const 
  {
    return mesh_->num_vertices();
  }

  //! Number of ghost cells in the mesh
  size_t num_ghost_cells() const 
  {
    return 0;
  }

  //! Number of ghost faces in the mesh
  size_t num_ghost_faces() const 
  {
    return 0;
  }

  //! Number of ghost nodes in the mesh
  size_t num_ghost_nodes() const 
  {
    return 0;
  }

  //! Number of items of given entity
  //! \param [in] entity  The enumerated entity of interest
  //! \param [in] entity_type   The type of entity information (ghost, shared, 
  //!   all, etc...) 
  size_t num_entities(
    entity_kind_t entity, 
    entity_type_t entity_type = entity_type_t::ALL
  ) const 
  {
    switch(entity) {
      case entity_kind_t::NODE :
        return num_owned_nodes();
      case entity_kind_t::EDGE :
        return num_owned_edges();
      case entity_kind_t::FACE :
        return num_owned_faces();
      case entity_kind_t::CELL :
        return num_owned_cells();
      case entity_kind_t::WEDGE :
        return mesh_->num_wedges();
        break;
      case entity_kind_t::CORNER :
        return mesh_->num_corners();
      default :
        raise_runtime_error("Unknown entity type");
        return 0;
    }
  }

  //! The begin iterator for iterating over mesh entities.
  //! \param [in] entity  The enumerated entity of interest
  //! \param [in] entity_type   The type of entity information (ghost, shared, 
  //!   all, etc...) 
  auto begin(
    entity_kind_t entity,
    entity_type_t entity_type = entity_type_t::ALL
  ) const 
  {
    int start_index = 0;
    return Portage::make_counting_iterator(start_index);
  }

  //! The end Iterator for iterating over mesh entities.
  //! \param [in] entity  The enumerated entity of interest
  //! \param [in] entity_type   The type of entity information (ghost, shared, 
  //!   all, etc...) 
  auto end(
    entity_kind_t entity,
    entity_type_t entity_type = entity_type_t::ALL
  ) const 
  {
    return begin(entity, entity_type) + num_entities(entity, entity_type);
  }

  //! Get list of nodes for a cell
  //! \param [in] cell_id  The id of the cell
  //! \param [in,out] nodes  The list of nodes to populate
  template< typename T >
  void cell_get_nodes(size_t cell_id, std::vector<T> * nodes) const 
  {
    auto c = cells_[cell_id];
    nodes->clear();
    for (auto v : mesh_->vertices(c))
      nodes->emplace_back( v.id() );
  }



  //! Get cell faces and the directions in which they are used
  //! \param [in] cell_id  The id of the cell in question
  //! \param [in,out] faces  The list of cell faces to populate
  //! \param [in,out] face_dirs  The list of face directions to populate
  template< typename T >
  void cell_get_faces_and_dirs(
    size_t cell_id, 
    std::vector<T> *faces,
    std::vector<T> *face_dirs
  ) const 
  {
    faces->clear();
    face_dirs->clear();
    auto c = cells_[cell_id];
    for (auto f : mesh_->faces(c)) {
      faces->emplace_back( f.id() );
      face_dirs->emplace_back( mesh_->cells(f)[0] == c ? 1 : -1 );
    }
  }

  //! Get nodes of a face
  //! \param [in] face_id  The id of the face in question
  //! \param [in,out] nodes  The list of face nodes to populate
  template< typename T >
  void face_get_nodes(
    size_t face_id, 
    std::vector<T> *nodes
  ) const 
  {
    auto f = faces_[face_id];
    nodes->clear();
    for (auto v : mesh_->vertices(f))
      nodes->emplace_back( v.id() );
    //std::reverse( nodes->begin(), nodes->end() );
  }

  //! Get node connected neighbors of cell
  //! \param [in] cell_id  The id of the cell in question
  //! \param [in] type  The type of indexes to include (ghost, shared, 
  //!   all, etc...) 
  //! \param [in,out] adj_cells  The list of cell neighbors to populate
  template< typename T >
  void cell_get_node_adj_cells(
    size_t cell_id,
    entity_type_t const type,
    std::vector<T> *adj_cells
  ) const 
  {
    auto this_cell = cells_[cell_id];
    adj_cells->clear();
    // Loop over all nodes of this cell
    for (auto node : mesh_->vertices(this_cell)) {
      // Loop over all cells associated with this node
      for (auto cell : mesh_->cells(node)) {
        if (cell != this_cell)
          adj_cells->emplace_back(cell.global_id());
      }
    }
  }

  //! \brief Get "adjacent" nodes of given node
  //!
  //! Get "adjacent" nodes of given node - nodes that share a common
  //! cell with given node
  //!
  //! \param [in] node_id  The node index
  //! \param [in] type  The type of indexes to include (ghost, shared, 
  //!   all, etc...) 
  //! \param [in,out] adj_nodes  The list of node neighbors to populate
  template< typename T >
  void node_get_cell_adj_nodes(
    size_t node_id,
    entity_type_t type,
    std::vector<T> *adj_nodes
  ) const 
  {
    auto this_node = vertices_[node_id];
    adj_nodes->clear();
    // Loop over cells associated with this node
    for (auto cell : mesh_->cells(this_node)) {
      // Loop over nodes of this cell
      for (auto node : mesh_->vertices(cell)) {
        if (this_node != node)
          adj_nodes->emplace_back(node.global_id());
      }
    }
  }

  //! @brief Get adjacent "dual cells" of a given "dual cell"
  //!
  //! \param [in] node_id  The node index
  //! \param [in] type  The type of indexes to include (ghost, shared, 
  //!   all, etc...) 
  //! \param [in,out] adj_nodes  The list of node neighbors to populate
  template < typename T >
  void dual_cell_get_node_adj_cells(
    size_t node_id,
    entity_type_t const type,
    std::vector<T> *adj_nodes
  ) const 
  {
    auto this_node = vertices_[node_id];
    adj_nodes->clear();
    // Loop over cells associated with this node
    for (auto cell : mesh_->cells(this_node)) {
      // Loop over the nodes associated with this cell
      for (auto node : mesh_->vertices(cell)) {
        if (this_node != node)
          adj_nodes->emplace_back(node.global_id());
      }
    }
  }



  //!  \brief Get the coords of a node
  //!  \param[in] node_id The ID of the node.
  //!  \param[in,out] pp The Portage::Point object containing the coordinate
  //!    information.
  //!
  //!  \remark FleCSI Burton specialization doesn't currently fully support 1D.
  void node_get_coordinates(
    size_t node_id, 
    point_1d_t * pp
  ) const 
  {
    auto v = vertices_[node_id];
    *pp = portage::make_1d_point( v->coordinates() );
  }

  void node_get_coordinates(
    size_t node_id, 
    point_2d_t * pp
  ) const 
  {
    auto v = vertices_[node_id];
    *pp = portage::make_2d_point( v->coordinates() );
  }

  void node_get_coordinates(
    size_t node_id, 
    point_3d_t * pp
  ) const 
  {
    auto v = vertices_[node_id];
    *pp = portage::make_3d_point( v->coordinates() );
  }


  //! \brief Get the coodinates of the nodes of a cell.
  //! \param[in] cell_id The ID of the cell.
  //! \param[in,out] point_list The vector of Portage::Point objects containing
  //!   the coordinates of a node.  The length of the vector is equal to the 
  //!   number of nodes in the cell with ID @c cellid.
  //! \remark FleCSI Burton specialization doesn't currently fully support 1D.
  void cell_get_coordinates(
    size_t const cell_id,
    std::vector<point_t> *point_list
  )  const 
  {
    // Get this cell object
    auto cell = cells_[cell_id];

    // Loop over vertices of this cell to get their coordinates
    point_list->clear();
    auto verts = mesh_->vertices(cell);
    for (auto v : verts)
      point_list->emplace_back( portage::make_point(v->coordinates()) );
  }

  //! \brief 2D version of coords of nodes of a dual cell
  //! \param[in] node_id The ID of the node or dual cell in the dual mesh.
  //! \param[in,out] pplist The vector of Portage::Point objects containing the
  //!   coordinates of a node in the dual mesh / cell in the regular mesh.  The
  //!   length of the vector is equal to the number of nodes in the dual mesh
  //!   cell with ID @c nodeid.
  //!  
  //! The vertices are ordered CCW. For node @c nodeid not on a
  //! boundary, the vector @c pplist starts with a random vertex, but it is 
  //! still ordered CCW. Use the dual_cell_coordinates_canonical_rotation() 
  //! function to rotate the @c pplist into a canonical (unique) form.
  //! 
  //! \todo worry about boundary cases
  void dual_cell_get_coordinates(
    size_t node_id,
    std::vector<point_t> *point_list
  ) const
  {
    raise_implemented_error("dual_cell_get_coordinates not implemented yet!");
  }

  //! \brief Centroid of a cell.
  //! \param[in]  cell_id The ID of the cell.
  //! \param[in,out] centroid The vector of coordinates of the cell @c cellid's
  //!   centroid.  The length of the vector is equal to the dimension of the 
  //!    mesh.
  void cell_centroid(
    size_t cell_id,
    point_t * centroid
  ) const 
  {
    auto this_cell = cells_[cell_id];
    *centroid = portage::make_point(this_cell->centroid());
  }

  //! \brief Centroid of a dual cell.
  //! \param[in] node_id The ID of the node in the normal mesh / cell in the 
  //!   dual mesh.
  //! \param[in,out] centroid The vector of coordinates of the node in the 
  //!   normal mesh / the cell in the dual mesh with ID @c nodeid.  The length 
  //!   of the vector is equal to the dimension of the mesh.
  //!
  //! \todo NOTE: THIS IS ASSUMED TO BE THE NODE COORDINATE BECAUSE
  //!   THE NODAL VARIABLES LIVE THERE, BUT FOR DISTORTED GRIDS, THE
  //!   NODE COORDINATED MAY NOT BE THE CENTROID OF THE DUAL CELL
  void dual_cell_centroid(
    size_t node_id,
    point_t *centroid) const 
  {
    auto this_node = vertices_[node_id];
    *centroid = portage::make_point(this_node->coordinates());
  }


  //============================================================================
  // Public Members Required By Portage
  // These are only needed in 3D
  //============================================================================

  //! \brief Get the type of the cell - PARALLEL_OWNED or PARALLEL_GHOST
  //!
  //! Assumes a 1-1 correspondence between integer values of the
  //! enum types to avoid switch statements
  //! 
  //! \param [in] cell_id  The index of the cell in question.
  auto cell_get_type(size_t cell_id) const 
  {
    return entity_type_t::PARALLEL_OWNED;
  }

  //! Get the element type of a cell - TRI, QUAD, POLYGON, TET, HEX,
  //! PRISM OR POLYHEDRON

  //! \brief Get the element type of a cell
  //!
  //! Can be one of: TRI, QUAD, POLYGON, TET, HEX,
  //!   PRISM OR POLYHEDRON.
  //! 
  //! \param [in] cell_id  The index of the cell in question.
  auto cell_get_element_type(size_t cell_id) const 
  {
    // search for the type in the map
    auto tp = cells_[cell_id]->type();
    auto it = shapes_to_portage.find( tp );
    // if it was found, return the mapped type
    if ( it != shapes_to_portage.end() )
      return it->second;
    // otherwise we dont know what it is
    return element_type_t::UNKNOWN_TOPOLOGY;
  }

  //============================================================================
  // Private Members
  //============================================================================

private:

  //! \brief a pointer to the mesh
  const mesh_t * mesh_ = nullptr;
  //! \brief the list of cells
  decltype( mesh_->cells() ) cells_;
  //! \brief the list of faces
  decltype( mesh_->faces() ) faces_;
  //! \brief the list of veritices
  decltype( mesh_->vertices() ) vertices_;

};

////////////////////////////////////////////////////////////////////////////////
// Static Initialization
////////////////////////////////////////////////////////////////////////////////

template< typename M >
const typename portage_mesh_t<M>::shape_map_t 
  portage_mesh_t<M>::shapes_to_portage = 
{ 
  { 
    shape_t::triangle, 
    element_type_t::TRI 
  },
  { 
    shape_t::quadrilateral, 
    element_type_t::QUAD 
  },
  { 
    shape_t::polygon, 
    element_type_t::POLYGON 
  },
  { 
    shape_t::tetrahedron, 
    element_type_t::TET
  },
  { 
    shape_t::hexahedron, 
    element_type_t::HEX 
  },
  { 
    shape_t::polyhedron, 
    element_type_t::POLYHEDRON 
  }
};


} // namespace
} // namespace

