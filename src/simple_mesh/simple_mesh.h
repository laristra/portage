/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_SIMPLE_MESH_SIMPLE_MESH_H_
#define SRC_SIMPLE_MESH_SIMPLE_MESH_H_

#include <vector>
#include <memory>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

/*!
  @file simple_mesh.h
  @brief A very light-weight, simple mesh infrastructure.
 */

namespace Portage {

/*!
  @class Simple_Mesh "simple_mesh.h"
  @brief A very light-weight, serial, 2D/3D Cartesian mesh.

  A Simple_Mesh is a non-distributed (i.e. serial), 2D/3D, regular Cartesian mesh.
  The user need only specify the domain extents and the number of cells per
  direction, and the mesh class will create all of the connectivity information.

  Simple_Mesh only knows about cells, faces, and nodes.  

  In a 3D mesh, each cell is a rectangular prism with 8 nodes and 6 faces. 
  Likewise, each non-domain-boundary node is connected to 8 cells and each 
  non-domain-boundary face is connected to two cells.  

  In a 2D mesh, each cell is a quadrilateral with 4 nodes and 4 faces. 
  Likewise, each non-domain-boundary node is connected to 4 cells and each 
  non-domain-boundary face is connected to two cells.  
  
  Nodes and faces on the domain boundary have less connections as there are no
  ghost mesh entities in Simple_Mesh.
 */
class Simple_Mesh {
  /// Convenience for denoting a mesh entity ID.
  typedef int ID;

 public:
  /*!
    @brief Constructor for creating a serial, 3D Cartesian mesh.
    @param[in] x0,y0,z0 The minimum coordinates of the domain.
    @param[in] x1,y1,z1 The maximum coordinates of the domain.
    @param[in] nx,ny,nz The number of _cells_ in each direction.

    By specifying the spatial extents and number of cells in each
    direction, we create a Cartesian mesh in three dimensions.
    Connectivity information is automatically built from global IDs.
    This mesh class has _zero_ ghost mesh entities.
  */
Simple_Mesh(double x0, double y0, double z0,
            double x1, double y1, double z1,
            int nx, int ny, int nz) :
            nx_(nx), ny_(ny), nz_(nz),
            x0_(x0), y0_(y0), z0_(z0),
            x1_(x1), y1_(y1), z1_(z1) {
    spacedim_ = 3;
    num_cells_ = nx*ny*nz;
    num_nodes_ = (nx+1)*(ny+1)*(nz+1);
    num_faces_ = (nx_+1)*ny_*nz_ + nx_*(ny_+1)*nz_ + nx_*ny_*(nz_+1);

    nodes_per_face_     = 4;
    nodes_per_cell_     = 8;
    faces_per_cell_     = 6;
    cells_per_node_aug_ = 9;

    // Construct the nodal coordinates from extents and number of nodes
    build_node_coords_3d();

    // Build cell <--> node, cell <--> face, face --> node adjacencies
    build_cfn_adjacencies_3d();

    // Build ownership information - no ghosts in Simple Mesh
    nodeids_owned_.resize(num_nodes_);
    for (int i(0); i < num_nodes_; ++i)
      nodeids_owned_[i] = i;

    cellids_owned_.resize(num_cells_);
    for (int i(0); i < num_cells_; ++i)
      cellids_owned_[i] = i;

    faceids_owned_.resize(num_faces_);
    for (int i(0); i < num_faces_; ++i)
      faceids_owned_[i] = i;
  }

  /*!
    @brief Constructor for creating a serial, 2D Cartesian mesh.
    @param[in] x0,y0 The minimum coordinates of the domain.
    @param[in] x1,y1 The maximum coordinates of the domain.
    @param[in] nx,ny The number of _cells_ in each direction.

    By specifying the spatial extents and number of cells in each
    direction, we create a Cartesian mesh in two dimensions.
    Connectivity information is automatically built from global IDs.
    This mesh class has _zero_ ghost mesh entities.
  */
  Simple_Mesh(double x0, double y0,
              double x1, double y1, 
              int nx, int ny) :
              nx_(nx), ny_(ny), 
              x0_(x0), y0_(y0),
              x1_(x1), y1_(y1)  {
    spacedim_ = 2; 
    num_cells_ = nx*ny;
    num_nodes_ = (nx+1)*(ny+1);
    num_faces_ = (nx_+1)*ny_ + nx_*(ny_+1);

    nodes_per_face_     = 2;
    nodes_per_cell_     = 4;
    faces_per_cell_     = 4;
    cells_per_node_aug_ = 5;

    // Construct the nodal coordinates from extents and number of nodes
    build_node_coords_2d();

    // Build cell <--> node, cell <--> face, face --> node adjacencies
    build_cfn_adjacencies_2d();

    // Build ownership information - no ghosts in Simple Mesh
    nodeids_owned_.resize(num_nodes_);
    for (int i(0); i < num_nodes_; ++i)
      nodeids_owned_[i] = i;

    cellids_owned_.resize(num_cells_);
    for (int i(0); i < num_cells_; ++i)
      cellids_owned_[i] = i;

    faceids_owned_.resize(num_faces_);
    for (int i(0); i < num_faces_; ++i)
      faceids_owned_[i] = i;
  }

  /// Assignment operator (disabled).
  Simple_Mesh & operator=(const Simple_Mesh &) = delete;

  /// Destructor
  ~Simple_Mesh() {}

  /// Spatial dimension of the mesh
  inline int space_dimension() const {
    return spacedim_;
  }

  /*!
    @brief Determine the number of a specific mesh entity.
    @param[in] kind The type of entity, e.g. @c CELL.
    @param[in] type The type of the entity, e.g. @c PARALLEL_OWNED
    @returns The number of the specified mesh entity.
   */
  int num_entities(const Entity_kind kind,
                   const Entity_type type) const {
    switch (kind) {
      case Entity_kind::NODE:
        switch (type) {
          case Entity_type::PARALLEL_OWNED:
            return nodeids_owned_.size();
          case Entity_type::PARALLEL_GHOST:
            // Simple_Mesh has no ghosts.
            return 0;
          case Entity_type::ALL:
            // Simple_Mesh has no ghosts.
            return nodeids_owned_.size();
          default:
            return 0;
        }
      case Entity_kind::CELL:
        switch (type) {
          case Entity_type::PARALLEL_OWNED:
            return cellids_owned_.size();
          case Entity_type::PARALLEL_GHOST:
            // Simple_Mesh has no ghosts.
            return 0;
          case Entity_type::ALL:
            // Simple_Mesh has no ghosts.
            return cellids_owned_.size();
          default:
            return 0;
        }
      case Entity_kind::FACE:
        switch (type) {
          case Entity_type::PARALLEL_OWNED:
            return faceids_owned_.size();
          case Entity_type::PARALLEL_GHOST:
            // Simple_Mesh has no ghosts.
            return 0;
          case Entity_type::ALL:
            // Simple_Mesh has no ghosts.
            return faceids_owned_.size();
          default:
            return 0;
        }
      default:
        return 0;
    }
  }

  /*!
    @brief For a given cell, get the list of faces and the direction of their
    normals.
    @param[in] cellid The ID of the cell.
    @param[out] faces The vector of face IDs corresponding to cell @c cellid.
    @param[out] fdirs The vector of face directions corresponding to each face
    in @c faces.
   */
  void cell_get_faces_and_dirs(const ID cellid,
                               std::vector<ID> *faces,
                               std::vector<int> *fdirs) const {
    auto offset = faces_per_cell_*cellid;
    faces->clear();
    fdirs->clear();
    for (int i(0); i < faces_per_cell_; ++i) {
      faces->push_back(cell_to_face_[i+offset]);
      fdirs->push_back(cell_face_dirs_[i+offset]);
    }
  }
  // @TODO: replace with std::copy?
  /*!
    @brief For a given cell, get the list of nodes.
    @param[in] cellid The ID of the cell.
    @param[out] nodes The vector of node IDs corresponding to cell @c cellid.
   */
  void cell_get_nodes(const ID cellid,
                      std::vector<ID> *nodes) const {
    auto offset = nodes_per_cell_*cellid;
    nodes->clear();
    for (int i(0); i < nodes_per_cell_; ++i)
      nodes->push_back(cell_to_node_[i+offset]);
  }
  /*!
    @brief For a given face, get the list of nodes.
    @param[in] faceid The ID of the face.
    @param[out] nodes The vector of node IDs corresponding to face @c faceid.
   */
  void face_get_nodes(const ID faceid,
                      std::vector<ID> *nodes) const {
    auto offset = nodes_per_face_*faceid;
    nodes->clear();
    for (int i(0); i < nodes_per_face_; ++i)
      nodes->push_back(face_to_node_[i+offset]);
  }

  /*!
    @brief For a given node, get all the cells attached to this node.
    @param[in] nodeid The ID of the node.
    @param[out] cells The vector of cell IDs attached to node @c nodeid.
   */
  void node_get_cells(const ID nodeid,
                      std::vector<ID> *cells) const {
    auto offset = cells_per_node_aug_*nodeid;
    cells->clear();
    for (int i(0); i < node_to_cell_[offset]; ++i) {
      cells->push_back(node_to_cell_[i+offset+1]);
    }
  }

  // General specification - specialization follows at bottom of file
  /*!
    @brief Get the coordinates of a node.
    @tparam D Dimension of the node.
    @param[in] nodeid The ID of the node.
    @param[out] pp The @c Point object of dimension @c D containing the
    coordinates of node @nodeid.

    This is the general specification.  
   */
  template<long D>
  void node_get_coordinates(const ID nodeid,
                            Point<D> *pp) const {
    assert(D == space_dimension());
  }

 private:
  /*!
    @brief Constructs and stores the node coordinates from the extents and
    number of cells per direction passed to the 3D constructor.
   */
  void build_node_coords_3d() {
    coordinates3d_.clear();

    double hx = (x1_ - x0_)/nx_;
    double hy = (y1_ - y0_)/ny_;
    double hz = (z1_ - z0_)/nz_;

    for (int iz(0); iz <= nz_; ++iz)
      for (int iy(0); iy <= ny_; ++iy)
        for (int ix(0); ix <= nx_; ++ix) {
          coordinates3d_.emplace_back(x0_+ix*hx,
                                    y0_+iy*hy,
                                    z0_+iz*hz);
        }
  }

  /*
    @brief Builds the cell-face-node adjacency information for 3D mesh.
   */
  void build_cfn_adjacencies_3d() {
    // downward adjacencies
    cell_to_node_.resize(nodes_per_cell_*num_cells_);
    cell_to_face_.resize(faces_per_cell_*num_cells_);
    cell_face_dirs_.resize(faces_per_cell_*num_cells_);
    face_to_node_.resize(nodes_per_face_*num_faces_);
    // upward adjacencies
    node_to_cell_.resize(cells_per_node_aug_*num_nodes_);
    face_to_cell_.resize(2*num_faces_, -1);

    // cell adjacencies
    for (int iz(0); iz < nz_; ++iz)
      for (int iy(0); iy < ny_; ++iy)
        for (int ix(0); ix < nx_; ++ix) {
          auto thisCell = cell_index_(ix, iy, iz);
          auto cstart = nodes_per_cell_ * thisCell;
          auto fstart = faces_per_cell_ * thisCell;

          /*
            Node ordering in hex:
                                  z
               7-------6          ^   y
              /|      /|          |  /
             / |     / |          | /
            4--|----5  |          |/
            |  3----|--2          +-------> x
            | /     | /
            0-------1
           */
          cell_to_node_[cstart  ] = node_index_(ix, iy, iz);
          cell_to_node_[cstart+1] = node_index_(ix+1, iy, iz);
          cell_to_node_[cstart+2] = node_index_(ix+1, iy+1, iz);
          cell_to_node_[cstart+3] = node_index_(ix, iy+1, iz);
          cell_to_node_[cstart+4] = node_index_(ix, iy, iz+1);
          cell_to_node_[cstart+5] = node_index_(ix+1, iy, iz+1);
          cell_to_node_[cstart+6] = node_index_(ix+1, iy+1, iz+1);
          cell_to_node_[cstart+7] = node_index_(ix, iy+1, iz+1);

          // loop over the 8 nodes attached to this cell, and assign
          // this cell to its cell connectivity
          // order shouldn't matter here
          for (int iiz(iz); iiz <= iz+1; ++iiz)
            for (int iiy(iy); iiy <= iy+1; ++iiy)
              for (int iix(ix); iix <= ix+1; ++iix) {
                auto thisNode = node_index_(iix, iiy, iiz);
                auto cnstart = cells_per_node_aug_ * thisNode;
                auto & c_at_n = node_to_cell_[cnstart];
                node_to_cell_[cnstart+c_at_n+1] = thisCell;
                c_at_n++;
              }

          /*
            Face ordering in hex:
                                 z
               +-------+         ^   y
              /|      /|         |  /
             / | 5  2/ |         | /
            +--|----+ 1|         |/
            |3 +0---|--+         +-------> x
            | /  4  | /
            +-------+

            The "dirs" indicate whether or not the nodes are listed in
            a ccw (dir=+1) or cw (dir=-1) orientation when look DOWN
            the cell's OUTWARD normal at that face.  Images of each
            face are below when we construct the face adjacencies.
           */
          cell_to_face_[fstart  ] = xzface_index_(ix, iy, iz);    // front
          cell_face_dirs_[fstart  ] = 1;
          cell_to_face_[fstart+1] = yzface_index_(ix+1, iy, iz);  // right
          cell_face_dirs_[fstart+1] = 1;
          cell_to_face_[fstart+2] = xzface_index_(ix, iy+1, iz);  // back
          cell_face_dirs_[fstart+2] = -1;
          cell_to_face_[fstart+3] = yzface_index_(ix, iy, iz);    // left
          cell_face_dirs_[fstart+3] = -1;
          cell_to_face_[fstart+4] = xyface_index_(ix, iy, iz);    // bottom
          cell_face_dirs_[fstart+4] = -1;
          cell_to_face_[fstart+5] = xyface_index_(ix, iy, iz+1);  // top
          cell_face_dirs_[fstart+5] = 1;

          // go over the 6 faces attached to this cell, and assign
          // this cell to its cell connectivity
          auto cfstart = 2*xzface_index_(ix, iy, iz) + 1;  // +y side of front
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*yzface_index_(ix+1, iy, iz);          // -x side of right
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xzface_index_(ix, iy+1, iz);          // -y side of back
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*yzface_index_(ix, iy, iz) + 1;        // +x side of left
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xyface_index_(ix, iy, iz) + 1;        // +z side of bottom
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xyface_index_(ix, iy, iz+1);          // -z side of top
          face_to_cell_[cfstart] = thisCell;
        }

    // face adjacencies
    /* xy faces
         3-------2
        /       /
       /       /
      0-------1
    */
    for (int iz(0); iz <= nz_; ++iz)
      for (int iy(0); iy < ny_; ++iy)
        for (int ix(0); ix < nx_; ++ix) {
          auto thisFace = xyface_index_(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node_[nstart  ] = node_index_(ix, iy, iz);
          face_to_node_[nstart+1] = node_index_(ix+1, iy, iz);
          face_to_node_[nstart+2] = node_index_(ix+1, iy+1, iz);
          face_to_node_[nstart+3] = node_index_(ix, iy+1, iz);
        }
    /* xz faces
       3-------2
       |       |
       |       |
       |       |
       0-------1
     */
    for (int iz(0); iz < nz_; ++iz)
      for (int iy(0); iy <= ny_; ++iy)
        for (int ix(0); ix < nx_; ++ix) {
          auto thisFace = xzface_index_(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node_[nstart  ] = node_index_(ix, iy, iz);
          face_to_node_[nstart+1] = node_index_(ix+1, iy, iz);
          face_to_node_[nstart+2] = node_index_(ix+1, iy, iz+1);
          face_to_node_[nstart+3] = node_index_(ix, iy, iz+1);
        }
    /* yz faces
          2
         /|
        / |
       3  |
       |  1
       | /
       |/
       0
     */
    for (int iz(0); iz < nz_; ++iz)
      for (int iy(0); iy < ny_; ++iy)
        for (int ix(0); ix <= nx_; ++ix) {
          auto thisFace = yzface_index_(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node_[nstart  ] = node_index_(ix, iy, iz);
          face_to_node_[nstart+1] = node_index_(ix, iy+1, iz);
          face_to_node_[nstart+2] = node_index_(ix, iy+1, iz+1);
          face_to_node_[nstart+3] = node_index_(ix, iy, iz+1);
        }
  }


 /*!
    @brief Constructs and stores the node coordinates from the extents and
    number of cells per direction passed to the 2D constructor.
   */
  void build_node_coords_2d() {
    coordinates2d_.clear();

    double hx = (x1_ - x0_)/nx_;
    double hy = (y1_ - y0_)/ny_;

    for (int iy(0); iy <= ny_; ++iy)
      for (int ix(0); ix <= nx_; ++ix) {
        coordinates2d_.emplace_back(x0_+ix*hx,
                                  y0_+iy*hy);
      }
  }

  /*
    @brief Builds the cell-face-node adjacency information for 2D mesh. 
   */
  void build_cfn_adjacencies_2d() {
    // downward adjacencies
    cell_to_node_.resize(nodes_per_cell_*num_cells_);
    cell_to_face_.resize(faces_per_cell_*num_cells_);
    cell_face_dirs_.resize(faces_per_cell_*num_cells_);
    face_to_node_.resize(nodes_per_face_*num_faces_);
    // upward adjacencies
    node_to_cell_.resize(cells_per_node_aug_*num_nodes_);
    face_to_cell_.resize(2*num_faces_, -1);

    // cell adjacencies
    for (int iy(0); iy < ny_; ++iy)
      for (int ix(0); ix < nx_; ++ix) {
        auto thisCell = cell_index_(ix, iy, 0);
        auto cstart = nodes_per_cell_ * thisCell;
        auto fstart = faces_per_cell_ * thisCell;

          /*
            Node ordering in quad:
                                 y
                                 ^    
            3-------2            |
            |       |            |
            |       |            |
            0-------1            +-------> x
           */

          cell_to_node_[cstart  ] = node_index_(ix, iy, 0);
          cell_to_node_[cstart+1] = node_index_(ix+1, iy, 0);
          cell_to_node_[cstart+2] = node_index_(ix+1, iy+1, 0);
          cell_to_node_[cstart+3] = node_index_(ix, iy+1, 0);


          // loop over the 4 nodes attached to this cell, and assign
          // this cell to its cell connectivity
          // order shouldn't matter here
          for (int iiy(iy); iiy <= iy+1; ++iiy)
            for (int iix(ix); iix <= ix+1; ++iix) {
              auto thisNode = node_index_(iix, iiy, 0);
              auto cnstart = cells_per_node_aug_ * thisNode;
              auto & c_at_n = node_to_cell_[cnstart];
              node_to_cell_[cnstart+c_at_n+1] = thisCell;
              c_at_n++;
            }

          /*
            Face ordering in quad:        
                                 y
                2                ^    
            +-------+            |
            |       |            |
          3 |       | 1          |        
            |       |            |
            +-------+            +-------> x              
                0

            The "dirs" indicate whether or not the nodes are listed in
            a ccw (dir=+1) or cw (dir=-1) orientation when look DOWN
            the cell's OUTWARD normal at that face.  Images of each
            face are below when we construct the face adjacencies.
           */
          cell_to_face_[fstart  ] = xface_index_(ix, iy);    // bottom
          cell_face_dirs_[fstart  ] = 1;
          cell_to_face_[fstart+1] = yface_index_(ix+1, iy);  // right
          cell_face_dirs_[fstart+1] = 1;
          cell_to_face_[fstart+2] = xface_index_(ix, iy+1);  // up
          cell_face_dirs_[fstart+2] = -1;
          cell_to_face_[fstart+3] = yface_index_(ix, iy);    // left
          cell_face_dirs_[fstart+3] = -1;


          // go over the 4 faces attached to this cell, and assign
          // this cell to its cell connectivity
          auto cfstart = 2*xface_index_(ix, iy);    // +y side of bottom
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xface_index_(ix, iy+1)+1;     // -y side of up
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*yface_index_(ix+1, iy)+1;     // -x side of right
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*yface_index_(ix, iy) ;        // +x side of left
          face_to_cell_[cfstart] = thisCell;
         
        }

    // face adjacencies

    for (int iy(0); iy <= ny_; ++iy)
      for (int ix(0); ix < nx_; ++ix) {
        auto thisFace = xface_index_(ix, iy);
        auto nstart = nodes_per_face_ * thisFace;
        face_to_node_[nstart  ] = node_index_(ix, iy, 0);
        face_to_node_[nstart+1] = node_index_(ix+1, iy, 0);
      }

    for (int iy(0); iy < ny_; ++iy)
      for (int ix(0); ix <= nx_; ++ix) {
        auto thisFace = yface_index_(ix, iy);
        auto nstart = nodes_per_face_ * thisFace;
        face_to_node_[nstart  ] = node_index_(ix, iy, 0);
        face_to_node_[nstart+1] = node_index_(ix, iy+1, 0);
      }
  
  }


  /// @c Simple_Mesh can be 2D or 3D 
  int spacedim_ ;

  /// Number of cells in the three coordinate directions.
  int nx_, ny_, nz_;

  /// Coordinates of lower left front and upper right back of domain.
  double x0_, x1_, y0_, y1_, z0_, z1_;

  /// Node positions.
  std::vector<Point<3>> coordinates3d_;
  std::vector<Point<2>> coordinates2d_;

  /// Set in constructor for 2D quads and 3D hexes for now.
  int nodes_per_face_ ;
  int nodes_per_cell_ ;
  int faces_per_cell_ ;
  int cells_per_node_aug_ ;   // 1 entry for the num cells actually attached

  /// Cache of stored sizes.
  int num_cells_;
  int num_nodes_;
  int num_faces_;

  /// Storage for connectivity information.
  std::vector<ID> cell_to_face_;
  std::vector<int> cell_face_dirs_;
  std::vector<ID> cell_to_node_;
  std::vector<ID> face_to_node_;
  std::vector<ID> face_to_cell_;
  std::vector<ID> node_to_cell_;

  /// Cache of entity ID lists.
  std::vector<ID> nodeids_owned_;
  std::vector<ID> faceids_owned_;
  std::vector<ID> cellids_owned_;

  /// Helper functions for looking up indices.
  ID node_index_(int i, int j, int k) const {
    return i + j*(nx_+1) + k*(nx_+1)*(ny_+1);
  }
  ID cell_index_(int i, int j, int k) const {
    return i + j*nx_ + k*nx_*ny_;
  }

  ID xface_index_(int i, int j) const {
    return i + j*nx_ ;
  }
  ID yface_index_(int i, int j) const {
    return i + j*(nx_+1) + xface_index_(0, ny_+1);
  }

  ID xyface_index_(int i, int j, int k) const {
    return i + j*nx_ + k*nx_*ny_;
  }
  ID xzface_index_(int i, int j, int k) const {
    return i + j*nx_ + k*nx_*(ny_+1) + xyface_index_(0, 0, nz_+1);
  }
  ID yzface_index_(int i, int j, int k) const {
    return i + j*(nx_+1) + k*(nx_+1)*ny_ + xzface_index_(0, 0, nz_);
  }
};  // class Simple_Mesh

// Specializations
/*!
  @brief Get the 2D/3D coordinates of a specific node as @c Portage::Point object.
  @param[in] nodeid The ID of the node.
  @param[out] pp The @c Portage::Point containing the coordinates for node
  @c nodeid.
 */
template<>
void Simple_Mesh::node_get_coordinates<3>(const ID nodeid,
                                          Point<3> *pp) const {
  assert(spacedim_ == 3);
  *pp = coordinates3d_[nodeid];
}

template<>
void Simple_Mesh::node_get_coordinates<2>(const ID nodeid,
                                          Point<2> *pp) const {
  assert(spacedim_ == 2);
  *pp = coordinates2d_[nodeid];
}


}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_MESH_H_
