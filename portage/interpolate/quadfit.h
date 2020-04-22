/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERPOLATE_QUADFIT_H_
#define PORTAGE_INTERPOLATE_QUADFIT_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <iostream>

// portage includes
#include "portage/support/portage.h"

// wonton includes
#include "wonton/support/lsfits.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

namespace Portage {

using Wonton::Point;
using Wonton::Vector;

///////////////////////////////////////////////////////////////////////////////


/*! @class Limited_Quadfit quadfit.h
    @brief Compute limited quadfit of a field or components of a field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
    @tparam on_what An enum type which indicates different entity types


*/

template<int D, Entity_kind on_what, typename MeshType, typename StateType>
class Limited_Quadfit {
 public:
  /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] on_what An enum that indicates what type of entity the field is on
      @param[in] var_name Name of field for which the quadfit is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)
      @param[in] Boundary_Limiter_type An enum indicating the limiter type on the boundary

      @todo must remove assumption that field is scalar
   */

  Limited_Quadfit(MeshType const & mesh, StateType const & state,
                   std::string const var_name,
                   Limiter_type limiter_type,
                   Boundary_Limiter_type Boundary_Limiter_type)
    : mesh_(mesh),
      state_(state),
      var_name_(var_name),
      limtype_(limiter_type),
      bnd_limtype_(Boundary_Limiter_type) {

    // Extract the field data from the statemanager

    state.mesh_get_data(on_what, var_name, &vals_);
  }

  /// @todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Quadfit(const Limited_Quadfit &) = delete;

  /// Assignment operator (disabled)

  Limited_Quadfit & operator = (const Limited_Quadfit &) = delete;

  /// Destructor

  ~Limited_Quadfit() = default;

  /// Functor - not implemented for all types - see specialization for
  /// cells, nodes

  Vector<D*(D+3)/2> operator()(int entity_id) {
    std::cerr << "Limited quadfit not implemented for this entity kind\n";
  }

 private:
  MeshType const & mesh_;
  StateType const & state_;
  std::string var_name_;
  double const* vals_;
  Limiter_type limtype_;
  Boundary_Limiter_type bnd_limtype_;
};


///////////////////////////////////////////////////////////////////////////////

/*! @class Limited_Quadfit<MeshType,StateType,CELL> quadfit.h
    @brief Specialization of limited quadfit class for @c cell-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
*/

template<int D, typename MeshType, typename StateType>
class Limited_Quadfit<D, Entity_kind::CELL, MeshType, StateType> {
 public:
  /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] var_name Name of field for which the quadfit is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)
      @param[in] Boundary_Limiter_type An enum indicating the limiter type on the boundary

      @todo must remove assumption that field is scalar
   */

  Limited_Quadfit(MeshType const & mesh, StateType const & state,
                   std::string const var_name,
                   Limiter_type limiter_type,
                   Boundary_Limiter_type Boundary_Limiter_type)
    : mesh_(mesh),
      state_(state),
      var_name_(var_name),
      limtype_(limiter_type),
      bnd_limtype_(Boundary_Limiter_type) {

    // Extract the field data from the statemanager
    state.mesh_get_data(Entity_kind::CELL, var_name, &vals_);

    // Collect and keep the list of neighbors for each NODE as it may
    // be expensive to go to the mesh layer and collect this data for
    // each cell during the actual quadfit calculation

    int ncells = mesh_.num_entities(Entity_kind::CELL);
    cell_neighbors_.resize(ncells);

    Wonton::for_each(mesh_.begin(Entity_kind::CELL), mesh_.end(Entity_kind::CELL),
                      [this](int c) { mesh_.cell_get_node_adj_cells(
                             c, Entity_type::ALL, &(cell_neighbors_[c])); } );
  }

  /// @todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Quadfit(const Limited_Quadfit &) = delete;

  /// Assignment operator (disabled)

  Limited_Quadfit & operator = (const Limited_Quadfit &) = delete;

  /// Destructor

  ~Limited_Quadfit() = default;

  /// Functor

  Vector<D*(D+3)/2> operator()(int cellid);

 private:
  MeshType const & mesh_;
  StateType const & state_;
  std::string var_name_;
  double const *vals_;
  Limiter_type limtype_;
  Boundary_Limiter_type bnd_limtype_;
  std::vector<std::vector<int>> cell_neighbors_;
};

  /*! @brief Implementation of Limited_Quadfit functor for CELLs
      Limited _Quadfit - Fit to a field to a quadratic
      multinomial using a Least-Squared fit. Returns an
      array of parameters.  If the CELL is on a boundary
      the stencil is too small, so it drops to linear order.
      Uses an SVD decomposition for the LS regression.

  */
template<int D, typename MeshType, typename StateType>
  Vector<D*(D+3)/2>
Limited_Quadfit<D, Entity_kind::CELL, MeshType, StateType>::operator() (int const cellid) {

  assert(D == mesh_.space_dimension());
  assert(D == 2 || D == 3);
  double phi = 1.0;
  Vector<D*(D+3)/2> qfit;
  Vector<D*(D+3)/2> dvec;

  bool boundary_cell =  mesh_.on_exterior_boundary(Entity_kind::CELL, cellid);
  // Limit the boundary gradient to enforce monotonicity preservation
  if (bnd_limtype_ == BND_ZERO_GRADIENT && boundary_cell) {
    qfit.zero();
    return qfit;
  }

  std::vector<int> const & nbrids = cell_neighbors_[cellid];

  std::vector<Point<D>> cellcenters(nbrids.size()+1);
  std::vector<double> cellvalues(nbrids.size()+1);

  // get centroid and value for cellid at center of point cloud
  mesh_.cell_centroid(cellid, &(cellcenters[0]));

  cellvalues[0] = vals_[cellid];

  int i = 1;
  for (auto nbrcell : nbrids) {
    mesh_.cell_centroid(nbrcell, &(cellcenters[i]));
    cellvalues[i] = vals_[nbrcell];
    i++;
  }

  qfit = Wonton::ls_quadfit(cellcenters, cellvalues, boundary_cell);
  // Limit the gradient to enforce monotonicity preservation

  if (limtype_ == BARTH_JESPERSEN && 
      (!boundary_cell || bnd_limtype_ == BND_BARTH_JESPERSEN)) {

    // Min and max vals of function (cell centered vals) among neighbors
    /// @todo: must remove assumption the field is scalar

    double minval = vals_[cellid];
    double maxval = vals_[cellid];

    int nnbr = nbrids.size();
    for (int ic = 0; ic < nnbr; ++ic) {
      minval = std::min(cellvalues[ic], minval);
      maxval = std::max(cellvalues[ic], maxval);
    }

    // Find the min and max of the reconstructed function in the cell
    // Since the reconstruction is linear, this will occur at one of
    // the nodes of the cell. So find the values of the reconstructed
    // function at the nodes of the cell

    double cellcenval = vals_[cellid];
    std::vector<Point<D>> cellcoords;
    mesh_.cell_get_coordinates(cellid, &cellcoords);

    for (auto coord : cellcoords) {
      Vector<D> vec = coord-cellcenters[0];
      //Vector<D*(D+3)/2> dvec;
      for (int j = 0; j < D; ++j) {
        dvec[j] = vec[j];
        // Add the quadratic terms
        for (int k = 0; k < j; ++k) {
          dvec[j+k+D-1] = dvec[k]*dvec[j];
        }
      }

      double diff = dot(qfit, dvec);
      double extremeval = (diff > 0.0) ? maxval : minval;
      double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff;
      phi = std::min(phi_new, phi);
    }

  }

  // Limited quadfit is phi*fit
  return phi*qfit;
}


///////////////////////////////////////////////////////////////////////////////

/*! @class Limited_Quadfit<MeshType,StateType,NODE> quadfit.h
    @brief Specialization of limited quadfit class for @c node-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
*/


template<int D, typename MeshType, typename StateType>
class Limited_Quadfit<D, Entity_kind::NODE, MeshType, StateType> {
 public:

  /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] var_name Name of field for which the quadfit is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)
      @param[in] Boundary_Limiter_type An enum indicating the limiter type on the boundary

      @todo must remove assumption that field is scalar
   */

  Limited_Quadfit(MeshType const & mesh, StateType const & state,
                   std::string const var_name,
                   Limiter_type limiter_type, 
                   Boundary_Limiter_type Boundary_Limiter_type)
    : mesh_(mesh),
      state_(state),
      var_name_(var_name),
      limtype_(limiter_type),
      bnd_limtype_(Boundary_Limiter_type) {

    // Extract the field data from the statemanager
    state.mesh_get_data(Entity_kind::NODE, var_name, &vals_);

    // Collect and keep the list of neighbors for each NODE as it may
    // be expensive to go to the mesh layer and collect this data for
    // each cell during the actual quadfit calculation

    int nnodes = mesh_.num_entities(Entity_kind::NODE);
    node_neighbors_.resize(nnodes);

    Wonton::for_each(mesh_.begin(Entity_kind::NODE), mesh_.end(Entity_kind::NODE),
                      [this](int n) { mesh_.dual_cell_get_node_adj_cells(
                             n, Entity_type::ALL, &(node_neighbors_[n])); } );
  }

  /// \todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Quadfit(const Limited_Quadfit &) = delete;

  /// Assignment operator (disabled)

  Limited_Quadfit & operator = (const Limited_Quadfit &) = delete;

  /// Destructor

  ~Limited_Quadfit() = default;

  /// Functor

  Vector<D*(D+3)/2> operator()(int nodeid);

 private:
  MeshType const & mesh_;
  StateType const & state_;
  std::string var_name_;
  double const *vals_;
  Limiter_type limtype_;
  Boundary_Limiter_type bnd_limtype_;
  std::vector<std::vector<int>> node_neighbors_;
};

  /*! @brief Implementation of Limited_Quadfit functor for NODEs
   *  Limited _Quadfit - Fit to a field to a quadratic
   *  multinomial using a Least-Squared fit. Returns an
   *  array of parameters.  If the MODE is on a boundary,
   *  the stencil is too small, so it drops to linear order.
   *  Uses an SVD decomposition for the LS regression.
   */


template<int D, typename MeshType, typename StateType>
  Vector<D*(D+3)/2>
Limited_Quadfit<D, Entity_kind::NODE, MeshType, StateType>::operator() (int const nodeid) {

  assert(D == mesh_.space_dimension());
  assert(D == 2 || D == 3);
  double phi = 1.0;
  Vector<D*(D+3)/2> qfit;
  Vector<D*(D+3)/2> dvec;

  bool boundary_node =  mesh_.on_exterior_boundary(Entity_kind::NODE, nodeid);
  if (bnd_limtype_ == BND_ZERO_GRADIENT && boundary_node) {
    qfit.zero();
    return qfit;
  }

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

  qfit = Wonton::ls_quadfit(nodecoords, nodevalues, boundary_node);

  if (limtype_ == BARTH_JESPERSEN && 
      (!boundary_node || bnd_limtype_ == BND_BARTH_JESPERSEN)) {

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
      // Vector<D*(D+3)/2> dvec;
	for (int j = 0; j < D; ++j) {
	  dvec[j] = vec[j];
	  // Add the quadratic terms
	  for (int k = 0; k < j; ++k) {
	    dvec[j+k+D-1] = dvec[k]*dvec[j];
	  }
	}

      double diff = dot(qfit, dvec);
      double extremeval = (diff > 0.0) ? maxval : minval;
      double phi_new = (diff == 0.0) ? 1 : (extremeval-nodeval)/diff;
      phi = std::min(phi_new, phi);
    }
  }

  // Limited gradient is phi*qfit

  return phi*qfit;
}

}  // namespace Portage

#endif  // PORTAGE_INTERPOLATE_QUADFIT_H_
