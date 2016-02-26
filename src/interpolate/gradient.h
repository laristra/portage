/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SRC_INTERPOLATE_GRADIENT_H_
#define SRC_INTERPOLATE_GRADIENT_H_

#include <algorithm>
#include <tuple>
#include <stdexcept>
#include <string>

#include "portage/support/matrix.h"

namespace Portage {

/// Limiter type

typedef enum {NOLIMITER, VAN_LEER, BARTH_JESPERSEN, MINMOD, SUPERBEE}
  LimiterType;



/*!
  @brief Compute least squares gradient from set of values
  @param[in] coords Vector of coordinates at which values are given
  @param[in] vals   Vector of values at said coordinates

  Compute a least squares gradient from a set of values. The first
  point is assumed to be the point where the gradient must be computed
  and the first value is assumed to the value at this reference point

  This operator does not know anything about a mesh.

*/

std::vector<double>
ls_gradient(std::vector<std::vector<double>> const & coords,
            std::vector<double> const & vals) {

  std::vector<double> coord0 = coords[0];
  int dim = coord0.size();

  double val0 = vals[0];

  // There are nvals but the first is the reference point where we
  // are trying to compute the gradient; so the matrix sizes etc
  // will only be nvals-1

  int nvals = vals.size();

  // Each row of A contains the components of the vector from
  // coord0 to the candidate point being used in the Least Squares
  // approximation (X_i-X_0).

  std::vector<std::vector<double>> A(nvals-1);
  for (int i = 0; i < nvals-1; ++i) {
    A[i].resize(dim);
    for (int j = 0; j < dim; ++j)
      A[i][j] = coords[i+1][j]-coord0[j];
  }


  // A is a matrix of size nvals-1 by D (where D is the space
  // dimension). So transpose(A)*A is D by D

  std::vector<std::vector<double>> AT;  // transpose(A)
  mat_transpose(A, &AT);

  std::vector<std::vector<double>> ATA;  // transpose(A)*A
  matmat_multiply(AT, A, &ATA);

  // Each entry/row of F contains the difference between the
  // function value at the candidate point and the function value
  // at the point where we are computing (f-f_0)

  std::vector<double> F(nvals-1);
  for (int i = 0; i < nvals-1; ++i)
    F[i] = vals[i+1]-val0;

  // F is a vector of nvals. So transpose(A)*F is vector of D
  // (where D is the space dimension)

  std::vector<double> ATF;  // transpose(A)*F
  matvec_multiply(AT, F, &ATF);

  // Inverse of ATA

  std::vector<std::vector<double>> ATAinv;
  mat_invert(ATA, &ATAinv);

  // Gradient of length D

  std::vector<double> grad;
  matvec_multiply(ATAinv, ATF, &grad);

  return grad;
}


/*! @class Limited_Gradient gradient.h
    @brief Compute limited gradient of a field or components of a field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
    @tparam on_what An enum type which indicates different entity types


*/

template<typename MeshType, typename StateType, Entity_kind on_what>
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

  std::vector<double> operator()(int entity_id) {
    std::cerr << "Limited gradient not implementd for this entity kind\n";
  }

 private:
  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  std::string const & var_name_;
  double * vals_;
};




///////////////////////////////////////////////////////////////////////////////

/*! @class Limited_Gradient<MeshType,StateType,CELL> gradient.h
    @brief Specialization of limited gradient class for @c cell-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
*/

template<typename MeshType, typename StateType>
class Limited_Gradient<MeshType, StateType, CELL> {
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

  std::vector<double> operator()(int cellid);

 private:
  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  std::string const & var_name_;
  double * vals_;
};

// @brief Implementation of Limited_Gradient functor for CELLs

template<typename MeshType, typename StateType>
std::vector<double>
Limited_Gradient<MeshType, StateType, CELL> :: operator() (int const cellid) {

  int dim = mesh_.space_dimension();
  double phi = 1.0;
  std::vector<double> grad(dim, 0);

  if (dim == 2 || dim == 3) {
    std::vector<int> nbrids;
    std::vector<std::vector<double>> cellcenters;
    std::vector<double> cellvalues;
    std::vector<double> cen, nbrcen;

    mesh_.cell_get_node_adj_cells(cellid, ALL, &nbrids);

    mesh_.cell_centroid(cellid, &cen);
    cellcenters.emplace_back(cen);
    cellvalues.emplace_back(vals_[cellid]);

    for (auto nbrcell : nbrids) {
      mesh_.cell_centroid(nbrcell, &nbrcen);
      cellcenters.emplace_back(nbrcen);
      cellvalues.emplace_back(vals_[nbrcell]);
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

      /// @todo: The following code is exactly why we need Point class
      double cellcenval = vals_[cellid];
      if (dim == 2) {

        std::vector<std::pair<double, double>> cellcoords;
        mesh_.cell_get_coordinates(cellid, &cellcoords);

        for (auto coord : cellcoords) {
  
          // At any coord:
          // val = cellcenval + grad*(coord-cellcencoord)
          // diff = val-cellcenval = grad DOT (coord-cellcencoord);
  
          double diff = grad[0]*(std::get<0>(coord)-cen[0]) +
              grad[1]*(std::get<1>(coord)-cen[1]);
          double extremeval = (diff > 0.0) ? maxval : minval;
          double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff;
          phi = std::min(phi_new, phi);
        }
      }
      else if (dim == 3) {

        std::vector<std::tuple<double, double, double>> cellcoords;
        mesh_.cell_get_coordinates(cellid, &cellcoords);

        for (auto coord : cellcoords) {
  
          // At any coord:
          // val = cellcenval + grad*(coord-cellcencoord)
          // diff = val-cellcenval = grad DOT (coord-cellcencoord);
  
          double diff = grad[0]*(std::get<0>(coord)-cen[0]) +
              grad[1]*(std::get<1>(coord)-cen[1]) +
            grad[2]*(std::get<2>(coord)-cen[2]);
          double extremeval = (diff > 0.0) ? maxval : minval;
          double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff;
          phi = std::min(phi_new, phi);
        }
      }
    }
  }
  else {
    std::cerr << "Gradient implemented only for 2D or 3D\n";
  }

  // Limited gradient is phi*grad

  for (int i = 0; i < dim; ++i)
    grad[i] *= phi;

  return grad;
}




///////////////////////////////////////////////////////////////////////////////

/*! @class Limited_Gradient<MeshType,StateType,NODE> gradient.h
    @brief Specialization of limited gradient class for @c node-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
*/


template<typename MeshType, typename StateType>
class Limited_Gradient<MeshType, StateType, NODE> {
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

  std::vector<double> operator()(int cellid);

 private:

  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  std::string const & var_name_;
  double * vals_;
};

// @brief Limited gradient functor implementation for NODE

template<typename MeshType, typename StateType>
std::vector<double>
Limited_Gradient<MeshType, StateType, NODE> :: operator() (int const nodeid) {

  int dim = mesh_.space_dimension();
  double phi = 1.0;
  std::vector<double> grad(dim, 0);

  std::vector<int> nbrids;
  mesh_.dual_cell_get_node_adj_cells(nodeid, ALL, &nbrids);

  if (dim ==2) {
    std::vector<std::vector<double>> nodecoords;
    std::vector<double> nodevalues;
    std::vector<double> ndcoord(3), coord(3);

    std::pair<double, double> coord_pair;

    ndcoord.resize(2);
    coord.resize(2);

    mesh_.node_get_coordinates(nodeid, &coord_pair);
    ndcoord[0] = coord_pair.first;
    ndcoord[1] = coord_pair.second;
    nodecoords.emplace_back(ndcoord);
    nodevalues.emplace_back(vals_[nodeid]);

    for (auto const & nbrnode : nbrids) {
      mesh_.node_get_coordinates(nbrnode, &coord_pair);
      coord[0] = coord_pair.first;
      coord[1] = coord_pair.second;
      nodecoords.emplace_back(coord);
      nodevalues.emplace_back(vals_[nbrnode]);
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

      std::vector<std::pair<double,double>> dualcellcoords;
      mesh_.dual_cell_get_coordinates(nodeid, &dualcellcoords);

      for (auto const & coord : dualcellcoords) {
        // val = nodeval + grad*(coord-nodecoord)
        // double diff = val-nodeval = grad DOT (coord-nodecoord);

        double diff = grad[0]*(std::get<0>(coord)-ndcoord[0]) +
            grad[1]*(std::get<1>(coord)-ndcoord[1]);
        double extremeval = (diff > 0.0) ? maxval : minval;
        double phi_new = (diff == 0.0) ? 1 : (extremeval-nodeval)/diff;
        phi = std::min(phi_new, phi);
      }
    }

  }
  else if (dim == 3) {
    std::vector<std::vector<double>> nodecoords;
    std::vector<double> nodevalues;
    std::vector<double> ndcoord(3), coord(3);

    std::tuple<double, double, double> coord_tuple;

    mesh_.node_get_coordinates(nodeid, &coord_tuple);
    ndcoord[0] = std::get<0>(coord_tuple);
    ndcoord[1] = std::get<1>(coord_tuple);
    ndcoord[2] = std::get<2>(coord_tuple);
    nodecoords.emplace_back(ndcoord);
    nodevalues.emplace_back(vals_[nodeid]);

    for (auto const & nbrnode : nbrids) {
      mesh_.node_get_coordinates(nbrnode, &coord_tuple);
      coord[0] = std::get<0>(coord_tuple);
      coord[1] = std::get<1>(coord_tuple);
      coord[2] = std::get<2>(coord_tuple);
      nodecoords.emplace_back(coord);
      nodevalues.emplace_back(vals_[nbrnode]);
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

      double nodeval = vals_[nodeid];

      std::vector<std::tuple<double, double, double>> dualcellcoords;
      mesh_.dual_cell_get_coordinates(nodeid, &dualcellcoords);

      for (auto const & coord : dualcellcoords) {
        // val = nodeval + grad*(coord-nodecoord)
        // diff = val-nodeval = grad DOT (coord-nodecoord);

        double diff = grad[0]*(std::get<0>(coord)-ndcoord[0]) +
            grad[1]*(std::get<1>(coord)-ndcoord[1]) +
            grad[2]*(std::get<2>(coord)-ndcoord[2]);
        double extremeval = (diff > 0.0) ? maxval : minval;
        double phi_new = (diff == 0.0) ? 1 : (extremeval-nodeval)/diff;
        phi = std::min(phi_new, phi);
      }
    }
  }
  else {
    std::cerr << "Gradient only implemented for 2D and 3D\n";
  }

  // Limited gradient is phi*grad

  for (int i = 0; i < dim; ++i)
    grad[i] *= phi;

  return grad;
}

}  // namespace Portage

#endif  // SRC_INTERPOLATE_GRADIENT_H_
