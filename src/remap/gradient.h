/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef PORTAGE_GRADIENT_H
#define PORTAGE_GRADIENT_H

#include "portage/support/matrix.h"

#include <algorithm>
#include <tuple>

namespace Portage {

/// Limiter type

typedef enum {NOLIMITER, VAN_LEER, BARTH_JESPERSEN, MINMOD, SUPERBEE} 
  LimiterType;
  


/*!
  \brief Compute least squares gradient from set of values

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
  
  std::vector<std::vector<double>> AT; // transpose(A)
  mat_transpose(A,&AT);
  
  std::vector<std::vector<double>> ATA;  // transpose(A)*A
  matmat_multiply(AT,A,&ATA);
  
  // Each entry/row of F contains the difference between the
  // function value at the candidate point and the function value
  // at the point where we are computing (f-f_0)
  
  std::vector<double> F(nvals-1);
  for (int i = 0; i < nvals-1; ++i)
    F[i] = vals[i+1]-val0;
  
  // F is a vector of nvals. So transpose(A)*F is vector of D
  // (where D is the space dimension)
  
  std::vector<double> ATF;  // transpose(A)*F
  matvec_multiply(AT,F,&ATF);
  
  // Inverse of ATA
  
  std::vector<std::vector<double>> ATAinv;
  mat_invert(ATA,&ATAinv);
  
  // Gradient of length D
  
  std::vector<double> grad;
  matvec_multiply(ATAinv,ATF,&grad);
  
  return grad;
}


/// Class for computing limited gradient of a field or components of a field

template<typename MeshType, typename StateType, typename OnWhatType>
class Limited_Gradient {
 public:

  //! Constructor

  Limited_Gradient(MeshType const & mesh, StateType const & state,
                   OnWhatType const on_what, std::string const remap_var_name,
                   LimiterType limiter_type) :
      mesh_(mesh), state_(state), on_what_(on_what), 
      remap_var_name_(remap_var_name), limtype_(limiter_type) {

    // Extract the field data from the statemanager
    
    state.get_data(on_what, remap_var_name, &vals_);
    
  }

  //! \todo Seems to be needed when using this in a Thrust transform call?
  //
  //  //! Copy constructor (deleted)
  //
  //  Limited_Gradient(const Limited_Gradient &) = delete;

  //! Assignment operator (disabled)

  Limited_Gradient & operator = (const Limited_Gradient &) = delete;

  //! Destructor

  ~Limited_Gradient() {}

  //! Functor

  std::vector<double> operator()(int cellid);

 private:

  LimiterType limtype_;
  MeshType const & mesh_;
  StateType const & state_;
  OnWhatType const on_what_;
  std::string const & remap_var_name_;
  double * vals_;  //! \todo must remove assumption that field is scalar
};
      


// Limited gradient functor of the Remap_2ndOrder class for each cell
//! To-do: May want to specialize on OnWhatType to avoid having to use
//! the "if on_what_ == " check

template<typename MeshType, typename StateType, typename OnWhatType>
std::vector<double>                          
Limited_Gradient<MeshType,StateType,OnWhatType> :: operator() (int const cellid) {

  std::vector<int> nbrids;
  std::vector<std::vector<double>> cellcenters;
  std::vector<double> cellvalues;
  std::vector<double> cen, nbrcen;

  //  if ((Entity_kind) on_what_ == CELL) {
    mesh_.cell_get_node_adj_cells(cellid,ALL,&nbrids);
    
    mesh_.cell_centroid(cellid,&cen);
    cellcenters.emplace_back(cen);
    cellvalues.emplace_back(vals_[cellid]);
    
    for (auto nbrcell : nbrids) {
      mesh_.cell_centroid(nbrcell,&nbrcen);
      cellcenters.emplace_back(nbrcen);
      cellvalues.emplace_back(vals_[nbrcell]);
    }
    //  }
  // else if ((Entity_kind) on_what_ == NODE) {
  //   // cellid is the ID of the dual cell which is the same as the node

  //   mesh_.dual_cell_get_node_adj_cells(cellid,ALL,&nbrids);

  //   mesh_.dual_cell_centroid(cellid,&cen);
  //   cellcenters.emplace_back(cen);
  //   cellvalues.emplace_back(vals_[cellid]);
    
  //   for (auto nbrcell : nbrids) {
  //     mesh_.dual_cell_centroid(nbrcell,&cen);
  //     cellcenters.emplace_back(cen);
  //     cellvalues.emplace_back(vals_[nbrcell]);
  //   }
  // }
  
  std::vector<double> grad = ls_gradient(cellcenters, cellvalues);

  
  // Limit the gradient to enforce monotonicity preservation

  double phi;
  if (limtype_ == NOLIMITER)
    phi = 1.0;
  else if (limtype_ == BARTH_JESPERSEN) {

    phi = 1.0;

    // Min and max vals of function (cell centered vals) among neighbors
    //! To-do: must remove assumption the field is scalar

    double minval = vals_[cellid];
    double maxval = vals_[cellid];

    int nnbr = nbrids.size();
    for (int i = 0; i < nnbr; ++i) {
      minval = std::min(cellvalues[i],minval);
      maxval = std::max(cellvalues[i],maxval);
    }

    // Find the min and max of the reconstructed function in the cell
    // Since the reconstruction is linear, this will occur at one of
    // the nodes of the cell. So find the values of the reconstructed
    // function at the nodes of the cell

    int dim = mesh_.space_dimension();

    //! To-do: The following code is exactly why we need Point class
    double cellcenval = vals_[cellid];
    if (dim == 1) {
    }
    else if (dim == 2) { 

      std::vector<std::pair<double,double>> cellcoords;
      //      if ((Entity_kind) on_what_ == CELL)
        mesh_.cell_get_coordinates(cellid,&cellcoords);
        //      else if ((Entity_kind) on_what_ == NODE)
        //        mesh_.dual_cell_get_coordinates(cellid,&cellcoords);

      for (auto coord : cellcoords) {
        // nodeval = cellcenval + grad*(nodecoord-cellcencoord)
        // diff = nodeval-cellcenval = grad DOT (nodecoord-cellcencoord);

        double diff = grad[0]*(std::get<0>(coord)-cen[0]) + 
            grad[1]*(std::get<1>(coord)-cen[1]);
        double extremeval = (diff > 0.0) ? maxval : minval;
        double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff; 
        phi = std::min(phi_new,phi);
      }      
    }
    else if (dim == 3) {

      std::vector<std::tuple<double,double,double>> cellcoords;
      //      if ((Entity_kind) on_what_ == CELL)
        mesh_.cell_get_coordinates(cellid,&cellcoords);
        //      else if ((Entity_kind) on_what_ == NODE)
        //        mesh_.dual_cell_get_coordinates(cellid,&cellcoords);

      for (auto coord : cellcoords) {
        // nodeval = cellcenval + grad*(nodecoord-cellcencoord)
        // diff = nodeval-cellcenval = grad DOT (nodecoord-cellcencoord);

        double diff = grad[0]*(std::get<0>(coord)-cen[0]) +
            grad[1]*(std::get<1>(coord)-cen[1]) +
            grad[2]*(std::get<2>(coord)-cen[2]);
        double extremeval = (diff > 0.0) ? maxval : minval;
        double phi_new = (diff == 0.0) ? 1 : (extremeval-cellcenval)/diff; 
        phi = std::min(phi_new,phi);
      }      
    }

  }
  else { // This will fail
    assert(limtype_ == NOLIMITER || limtype_ == BARTH_JESPERSEN);
  }



  // Limited gradient is phi*grad

  int dim = grad.size();
  for (int i = 0; i < dim; ++i)
    grad[i] *= phi;
    
  return grad;
}

} // namespace Portage

#endif
