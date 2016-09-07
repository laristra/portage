/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SRC_INTERPOLATE_GRADIENT_H_
#define SRC_INTERPOLATE_GRADIENT_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/Matrix.h"

namespace Portage {


/*!
  @brief Compute least squares gradient from set of values
  @param[in] coords Vector of coordinates at which values are given
  @param[in] vals   Vector of values at said coordinates

  Compute a least squares gradient from a set of values. The first
  point is assumed to be the point where the gradient must be computed
  and the first value is assumed to the value at this reference point

  This operator does not know anything about a mesh.

*/

template<long D>
Vector<D> ls_gradient(std::vector<Point<D>> const & coords,
                      std::vector<double> const & vals) {

  Point<D> coord0 = coords[0];

  double val0 = vals[0];

  // There are nvals but the first is the reference point where we
  // are trying to compute the gradient; so the matrix sizes etc
  // will only be nvals-1

  int nvals = vals.size();

  // Each row of A contains the components of the vector from
  // coord0 to the candidate point being used in the Least Squares
  // approximation (X_i-X_0).

  Matrix A(nvals-1, D);
  for (int i = 0; i < nvals-1; ++i) {
    for (int j = 0; j < D; ++j)
      A[i][j] = coords[i+1][j]-coord0[j];
  }


  // A is a matrix of size nvals-1 by D (where D is the space
  // dimension). So transpose(A)*A is D by D

  Matrix AT = A.transpose();

  Matrix ATA = AT*A;

  // Each entry/row of F contains the difference between the
  // function value at the candidate point and the function value
  // at the point where we are computing (f-f_0)

  std::vector<double> F(nvals-1);
  for (int i = 0; i < nvals-1; ++i)
    F[i] = vals[i+1]-val0;

  // F is a vector of nvals. So transpose(A)*F is vector of D
  // (where D is the space dimension)

  Vector<D> ATF = Vector<D>(AT*F);

  // Inverse of ATA

  Matrix ATAinv = ATA.inverse();

  // Gradient of length D

  return ATAinv*ATF;
}


/*! @class Limited_Gradient gradient.h
    @brief Compute limited gradient of a field or components of a field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
    @tparam on_what An enum type which indicates different entity types


*/

template<typename MeshType, typename StateType, Entity_kind on_what,
         long D>
class Limited_Gradient {
 public:
  /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] on_what An enum that indicates what type of entity the field is on
      @param[in] var_name Name of field for which the gradient is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)

      @todo must remove assumption that field is scalar
   */

  Limited_Gradient(MeshType const & mesh, StateType const & state,
                   std::string const var_name,
                   LimiterType limiter_type) :
      mesh_(mesh), state_(state),
      var_name_(var_name), limtype_(limiter_type) {

    // Extract the field data from the statemanager

    state.get_data(on_what, var_name, &vals_);
  }

  /// @todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Gradient(const Limited_Gradient &) = delete;

  /// Assignment operator (disabled)

  Limited_Gradient & operator = (const Limited_Gradient &) = delete;

  /// Destructor

  ~Limited_Gradient() {}

  /// Functor - not implemented for all types - see specialization for
  /// cells, nodes

  Vector<D> operator()(int entity_id) {
    std::cerr << "Limited gradient not implementd for this entity kind\n";
  }

 private:
  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  std::string var_name_;
  double const *vals_;
};




///////////////////////////////////////////////////////////////////////////////

/*! @class Limited_Gradient<MeshType,StateType,CELL> gradient.h
    @brief Specialization of limited gradient class for @c cell-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
*/

template<typename MeshType, typename StateType, long D>
class Limited_Gradient<MeshType, StateType, CELL, D> {
 public:
  /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] var_name Name of field for which the gradient is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)

      @todo must remove assumption that field is scalar
   */

  Limited_Gradient(MeshType const & mesh, StateType const & state,
                   std::string const var_name,
                   LimiterType limiter_type) :
      mesh_(mesh), state_(state), var_name_(var_name),
      limtype_(limiter_type) {

    // Extract the field data from the statemanager
    state.get_data(CELL, var_name, &vals_);

    // Collect and keep the list of neighbors for each NODE as it may
    // be expensive to go to the mesh layer and collect this data for
    // each cell during the actual gradient calculation

    int ncells = mesh_.num_entities(CELL);
    cell_neighbors_.resize(ncells);

    for (int c = 0; c < ncells; ++c)
      mesh_.cell_get_node_adj_cells(c, ALL, &(cell_neighbors_[c]));
  }

  /// @todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Gradient(const Limited_Gradient &) = delete;

  /// Assignment operator (disabled)

  Limited_Gradient & operator = (const Limited_Gradient &) = delete;

  /// Destructor

  ~Limited_Gradient() {}

  /// Functor

  Vector<D> operator()(int cellid);

 private:
  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  std::string var_name_;
  double const *vals_;
  std::vector<std::vector<int>> cell_neighbors_;
};

// @brief Implementation of Limited_Gradient functor for CELLs

template<typename MeshType, typename StateType, long D>
Vector<D>
Limited_Gradient<MeshType, StateType, CELL, D>::operator() (int const cellid) {

  assert(D == mesh_.space_dimension());
  assert(D == 2 || D == 3);
  double phi = 1.0;
  Vector<D> grad;

  std::vector<int> const & nbrids = cell_neighbors_[cellid];

  std::vector<Point<D>> cellcenters(nbrids.size()+1);
  std::vector<double> cellvalues(nbrids.size()+1);
  
  mesh_.cell_centroid(cellid, &(cellcenters[0]));
  cellvalues[0] = vals_[cellid];
  
  int i = 1;
  for (auto nbrcell : nbrids) {
    mesh_.cell_centroid(nbrcell, &(cellcenters[i]));
    cellvalues[i] = vals_[nbrcell];
    i++;
  }

  grad = ls_gradient(cellcenters, cellvalues);


  // Limit the gradient to enforce monotonicity preservation
  
  if (limtype_ == BARTH_JESPERSEN) {
    
    phi = 1.0;
    
    // Min and max vals of function (cell centered vals) among neighbors
    /// @todo: must remove assumption the field is scalar
    
    double minval = vals_[cellid];
    double maxval = vals_[cellid];
    
    int nnbr = nbrids.size();
    for (int i = 0; i < nnbr; ++i) {
      minval = std::min(cellvalues[i], minval);
      maxval = std::max(cellvalues[i], maxval);
    }
    
    // Find the min and max of the reconstructed function in the cell
    // Since the reconstruction is linear, this will occur at one of
    // the nodes of the cell. So find the values of the reconstructed
    // function at the nodes of the cell
    
    int dim = mesh_.space_dimension();
    
    double cellcenval = vals_[cellid];
    std::vector<Point<D>> cellcoords;
    mesh_.cell_get_coordinates(cellid, &cellcoords);
    
    for (auto coord : cellcoords) {
      Vector<D> vec = coord-cellcenters[0];
      double diff = dot(grad, vec);
      double extremeval = (diff > 0.0) ? maxval : minval;
      double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff;
      phi = std::min(phi_new, phi);
    }
  }

  // Limited gradient is phi*grad

  return phi*grad;
}




///////////////////////////////////////////////////////////////////////////////

/*! @class Limited_Gradient<MeshType,StateType,NODE> gradient.h
    @brief Specialization of limited gradient class for @c node-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
*/


template<typename MeshType, typename StateType, long D>
class Limited_Gradient<MeshType, StateType, NODE, D> {
 public:

  /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] var_name Name of field for which the gradient is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)

      @todo must remove assumption that field is scalar
   */

  Limited_Gradient(MeshType const & mesh, StateType const & state,
                   std::string const var_name,
                   LimiterType limiter_type) :
      mesh_(mesh), state_(state), var_name_(var_name),
      limtype_(limiter_type) {

    // Extract the field data from the statemanager
    state.get_data(NODE, var_name, &vals_);

    // Collect and keep the list of neighbors for each NODE as it may
    // be expensive to go to the mesh layer and collect this data for
    // each cell during the actual gradient calculation

    int nnodes = mesh_.num_entities(NODE);
    node_neighbors_.resize(nnodes);

    for (int n = 0; n < nnodes; ++n)
      mesh_.dual_cell_get_node_adj_cells(n, ALL, &(node_neighbors_[n]));
  }

  /// \todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Gradient(const Limited_Gradient &) = delete;

  /// Assignment operator (disabled)

  Limited_Gradient & operator = (const Limited_Gradient &) = delete;

  /// Destructor

  ~Limited_Gradient() {}

  /// Functor

  Vector<D> operator()(int cellid);

 private:

  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  std::string const & var_name_;
  double const *vals_;
  std::vector<std::vector<int>> node_neighbors_;
};

// @brief Limited gradient functor implementation for NODE

template<typename MeshType, typename StateType, long D>
Vector<D>
Limited_Gradient<MeshType, StateType, NODE, D>::operator() (int const nodeid) {

  assert(D == mesh_.space_dimension());
  assert(D == 2 || D == 3);
  double phi = 1.0;
  Vector<D> grad;

  std::vector<int> const & nbrids = node_neighbors_[nodeid];
  
  std::vector<Point<D>> nodecoords(nbrids.size()+1);
  std::vector<double> nodevalues(nbrids.size()+1);
  
  Point<D> point;
  
  mesh_.node_get_coordinates(nodeid, &(nodecoords[0]));
  nodevalues[0] = vals_[nodeid];
  
  int i = 1;
  for (auto const & nbrnode : nbrids) {
    mesh_.node_get_coordinates(nbrnode, &nodecoords[i]);
    nodevalues[i] = vals_[nbrnode];
    i++;
  }
  
  grad = ls_gradient(nodecoords, nodevalues);
  
  if (limtype_ == BARTH_JESPERSEN) {
    
    // Min and max vals of function (cell centered vals) among neighbors
    
    double minval = vals_[nodeid];
    double maxval = vals_[nodeid];
    
    for (auto const & val : nodevalues) {
      minval = std::min(val, minval);
      maxval = std::max(val, maxval);
    }
    
    // Find the min and max of the reconstructed function in the cell
    // Since the reconstruction is linear, this will occur at one of
    // the nodes of the cell. So find the values of the reconstructed
    // function at the nodes of the cell
    
    double nodeval = vals_[nodeid];
    
    std::vector<Point<D>> dualcellcoords;
    mesh_.dual_cell_get_coordinates(nodeid, &dualcellcoords);
    
    for (auto const & coord : dualcellcoords) {
      Vector<D> vec = coord-nodecoords[0];
      double diff = dot(grad, vec);
      double extremeval = (diff > 0.0) ? maxval : minval;
      double phi_new = (diff == 0.0) ? 1 : (extremeval-nodeval)/diff;
      phi = std::min(phi_new, phi);
    }
  }

  // Limited gradient is phi*grad

  return phi*grad;
}

}  // namespace Portage

#endif  // SRC_INTERPOLATE_GRADIENT_H_
