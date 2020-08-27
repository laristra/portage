/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <cmath>
#include <cassert>
#include <list>
#include <limits>
#include <vector>
#include "pairs.hh"

namespace Portage { namespace Meshfree { namespace Pairs {

/**
 * @brief Build structure for neighbor finding.
 *
 * Good old reliable method of cells in arbitrary dimensions. fastest.
 * Pairs arranged by x-point: all pairs for given x are contiguous.
 * Perfect linear scaling with number of points, but not efficient
 * for large variation in h and numerous isolated bunches in x or y
 * results of a typical optimized timing run for this method: (Mathematica input).
 *
 * (* 2D: equal size x and y boxes *)
 * (* finding containments for 64 per side, 4096 points.*)
 * (* finding containments for 128 per side, 16384 points.*)
 * (* finding containments for 256 per side, 65536 points.*)
 * (* finding containments for 512 per side, 262144 points.*)
 * times={{4096,0.05},{16384,0.17},{65536,0.68},{262144,2.96}};
 * ltimes = N@Log/@times;
 * ListPlot[ltimes, PlotJoined->True];
 * Print[Fit[ltimes, {1, x}, {x}]];
 * ptimes={{15896,0.05},{65414,0.17},{260597,0.68},{1.04724e+06,2.96}};
 * pltimes = N@Log/@ptimes;
 * ListPlot[pltimes, PlotJoined->True];
 * Print[Fit[pltimes, {1, x}, {x}]];
 * (* 3D: equal size x and y boxes *)
 * (* finding containments for 8 per side, 512 points.*)
 * (* finding containments for 16 per side, 4096 points.*)
 * (* finding containments for 32 per side, 32768 points.*)
 * times={{512,0.02},{4096,0.12},{32768,1}};
 * ltimes = N@Log/@times;
 * ListPlot[ltimes, PlotJoined->True];
 * Print[Fit[ltimes, {1, x}, {x}]];
 * ptimes={{4952,0.02},{34954,0.12},{277243,1}};
 * pltimes = N@Log/@ptimes;
 * ListPlot[pltimes, PlotJoined->True];
 * Print[Fit[pltimes, {1, x}, {x}]];
 *
 * @param x: contained points.
 * @param y: containing points with attached box.
 * @param h: box sizes, centered on y.
 * @param do_scatter: specify whether scatter or gather form.
 */
void CellPairFinder::init(vpile const& in_x, vpile const& in_y,
                          vpile const& in_h, bool in_do_scatter) {
  x = in_x;
  y = in_y;
  h = in_h;
  dim = in_x.size()[0];
  do_scatter = in_do_scatter;

  // get sizes and check
  ulong nx = x.size()[1];
  ulong ny = y.size()[1];
  assert(y.size()[0] == unsigned(dim));
  assert(h.size()[0] >= unsigned(dim));
  assert(h.size()[1] == (do_scatter ? nx : ny));

  // find min and max of enclosing boxes to get bounding box for cells
  cminmax = (do_scatter ? PairsMinMax(x, h) : PairsMinMax(y, h));
  delta = cminmax[1] - cminmax[0];

  // decide on size of grid - this is a maximum value
  ulong maxmemory = (do_scatter ? nx / 10 : ny / 10);
  size_t nsidemax = std::max<size_t>(static_cast<size_t>(ceil(pow(maxmemory * 1., (1. / dim)))), 1);

  // find avg of h
  pile havg(dim);
  havg = 0.;
  for (int m = 0; m < dim; m++)
    havg[m] += h[m].apply(fabs).sum();

  double const epsilon = std::numeric_limits<double>::epsilon();
  havg = (do_scatter ? 4. * havg / (1. * nx) + epsilon
                     : 4. * havg / (1. * ny) + epsilon);

  vulong nsideh(dim);
  for (int m = 0; m < dim; m++)
    nsideh[m] = static_cast<size_t>(ceil(delta[m] / havg[m]));

  // set number of cells each coordinate direction
  vulong nsides(dim);
  for (int m = 0; m < dim; m++)
    nsides[m] = std::min<ulong>(nsidemax, nsideh[m]);

  size_t ncells = 1;
  for (int m = 0; m < dim; m++)
    ncells *= nsides[m];

  nsidesm = nsides - static_cast<ulong>(1);

  // allocate grid of lists
  cells.resize(ncells);

  // strides
  strides.resize(dim);
  strides[dim - 1] = 1;
  for (int m = dim - 2; m >= 0; m--) {
    strides[m] = strides[m + 1] * nsides[m + 1];
  }

  // hash x-corners for scatter form
  if (do_scatter) {
    for (ulong i = 0; i < nx; i++) {
      // get lower left and upper right bounds and indices for source box
      pile xll(dim), xur(dim);
      for (int m = 0; m < dim; m++) {
        xll[m] = x[m][i] - 2. * h[m][i];
        xur[m] = x[m][i] + 2. * h[m][i];
      }
      vulong ixl(dim), ixu(dim);
      PairsIntegize(dim, xll, cminmax, delta, nsidesm, ixl);
      PairsIntegize(dim, xur, cminmax, delta, nsidesm, ixu);
      size_t ndxll = cellindex(dim, strides, ixl);
      size_t ndxur = cellindex(dim, strides, ixu) + 1;

      // add x to all cells covered or intersected by this source box
      for (size_t ndx = ndxll; ndx < ndxur; ndx++) {
        cells[ndx].push_back(i);
      }
    }
  } else {
    // hash x-points for gather form
    for (ulong i = 0; i < nx; i++) {
      // ignore x values outside bounding box
      bool outside = false;
      pile xi(dim);
      for (int m = 0; m < dim; m++) {
        if (x[m][i] <= cminmax[0][m] or x[m][i] >= cminmax[1][m]) { outside = true; }
        xi[m] = x[m][i];
      }
      if (outside) continue;

      // get cell indices this x
      vulong ix(dim);
      PairsIntegize(dim, xi, cminmax, delta, nsidesm, ix);

      // add this x to cell list
      ulong ndx = cellindex(dim, strides, ix);
      cells[ndx].push_back(i);
    }
  }
}

/**
 * @brief Get pairs for target point j, gather case.
 *
 * @param j: current point index.
 * @return a lst of neighbors of the j-th point.
 */
std::list<ulong> CellPairFinder::find_gather(const ulong j) const {
  std::list<ulong> pairlist;

  // get cell indices lower left and upper right corners of box
  pile yll(dim), yur(dim);
  for (int m = 0; m < dim; m++) {
    yll[m] = y[m][j] - 2. * h[m][j];
    yur[m] = y[m][j] + 2. * h[m][j];
  }
  vulong iyl(dim), iyu(dim);
  PairsIntegize(dim, yll, cminmax, delta, nsidesm, iyl);
  PairsIntegize(dim, yur, cminmax, delta, nsidesm, iyu);

  // total number of cells for this y
  ulong ncellsy = iyu[dim - 1] - iyl[dim - 1] + 1;
  vulong ystrides(dim);
  ystrides[dim - 1] = 1;
  for (int m = dim - 2; m >= 0; m--) {
    ystrides[m] = ncellsy;
    ncellsy *= iyu[m] - iyl[m] + 1;
  }

  // scan cells for this y
  for (ulong cell = 0; cell < ncellsy; cell++) {
    // convert local y-cell indices to global cell index
    vulong yndx(dim);
    cellindices(dim, ystrides, cell, yndx);
    vulong cellis(dim);
    for (int m = 0; m < dim; m++) cellis[m] = iyl[m] + yndx[m];
    size_t celli = cellindex(dim, strides, cellis);

    // determine if in interior or boundary of y-cell
    bool ybndry = false;
    for (int m = 0; m < dim; m++)
      if (cellis[m] == iyl[m] || cellis[m] == iyu[m]) ybndry = true;

    // loop over all x's in this cell's list
    for (auto&& i : cells[celli]) {
      // if on y-cell boundary, check that x's are contained
      bool inside = true;
      if (ybndry) {
        for (int m = 0; m < dim; m++) {
          if (x[m][i] <= yll[m]) inside = false;
          if (x[m][i] >= yur[m]) inside = false;
        }
      }

      // add pair: put x's in this y-cell onto neighbor list, if inside
      if (inside) {
        pairlist.push_back(i);
      }
    }  // for i
  }  // for cell

  return pairlist;
}  // CellPairFinder::find_gather

/**
 * @brief Get pairs for target point j, scatter case.
 *
 * @param j: current point index.
 * @return a lst of neighbors of the j-th point.
 */
std::list<ulong> CellPairFinder::find_scatter(const ulong j) const {
  std::list<ulong> pairlist;

  // get a compact representation of this point
  pile ypt(dim);
  for (int m = 0; m < dim; m++)
    ypt[m] = y[m][j];

  // check for completely outside source boxes
  bool outside = false;
  for (int m = 0; m < dim; m++) {
    if (y[m][j] <= cminmax[0][m]) {
      outside = true;
      break;
    }
    if (y[m][j] >= cminmax[1][m]) {
      outside = true;
      break;
    }
  }
  if (outside)
    return pairlist;

  // get cell indices of input y-point
  vulong iy(dim);
  PairsIntegize(dim, ypt, cminmax, delta, nsidesm, iy);
  size_t ndx = cellindex(dim, strides, iy);

  // loop over all x's in this y-cell's list
  for (auto&& i : cells[ndx]) {
    // get lower left and upper right coords for this x
    pile xll(dim), xur(dim);
    for (int m = 0; m < dim; m++) {
      xll[m] = x[m][i] - 2. * h[m][i];
      xur[m] = x[m][i] + 2. * h[m][i];
    }

    // check that y is contained
    bool inside = true;
    for (int m = 0; m < dim; m++) {
      if (y[m][j] <= xll[m]) {
        inside = false;
        break;
      }
      if (y[m][j] >= xur[m]) {
        inside = false;
        break;
      }
    }

    // add pair: put x's in this y-cell onto neighbor list, if inside
    if (inside) {
      pairlist.push_back(i);
    }
  }  // for i

  return pairlist;
}  // CellPairFinder::find_scatter


/**
 * @brief Find bounding box of all the y-boxes, variable h
 *
 * @param in_y
 * @param in_h
 * @return
 */
vpile CellPairFinder::PairsMinMax(const vpile &in_y, const vpile &in_h) const {
  std::vector<size_t> ny(in_y.size());  // sizes
  int const dimension = ny[0];
  vpile r(2, dimension);

  // find mins and maxes
  double shift = 2.0 * std::numeric_limits<double>::epsilon();  // smallest increment
  for (int m = 0; m < dimension; m++) {
    r[0][m] = (in_y[m] - 2. * in_h[m]).min();
    r[1][m] = (in_y[m] + 2. * in_h[m]).max();
  }
  r[0] -= shift;
  r[1] += shift;
  return r;
}

/**
 * @brief Find bounding box of all the y-boxes
 *
 * @param in_c
 * @param in_h
 * @return
 */
vpile CellPairFinder::PairsMinMax(const vpile &in_c, const pile &in_h) const {
  std::vector<size_t> nc(in_c.size());  // sizes
  int const dimension = nc[0];
  vpile r(2, dimension);

  // find mins and maxes
  double shift = 2.0 * std::numeric_limits<double>::epsilon();  // smallest increment
  for (int m = 0; m < dimension; m++) {
    r[0][m] = in_c[m].min() - 2. * in_h[m];
    r[1][m] = in_c[m].max() + 2. * in_h[m];
  }
  r[0] -= shift;
  r[1] += shift;
  return r;
}

/**
 * @brief convert real coordinates scaled to unit box to integers in [0,ulong_max]
 *
 * @param in_dim
 * @param in_value
 * @param in_minmax
 * @param in_delta
 * @param in_sizes
 * @param in_ivalue
 */
void CellPairFinder::PairsIntegize(int in_dim, const pile &in_value,
                                   const vpile &in_minmax, const pile &in_delta,
                                   const vulong &in_sizes, vulong &in_ivalue) const {
  for (int m = 0; m < in_dim; m++) {
    auto value = in_sizes[m] * (in_value[m] - in_minmax[0][m]) / in_delta[m];
    value = static_cast<ulong>(std::floor(value));
    in_ivalue[m] = std::max<ulong>(std::min<ulong>(value, in_sizes[m] - 1), 0);
  }
}

/**
 * @brief Get index from indices.
 *
 * @param in_dim
 * @param in_strides
 * @param in_indices
 * @return
 */
ulong CellPairFinder::cellindex(int in_dim, const vulong &in_strides,
                                const vulong &in_indices) const {
  ulong result = 0;
  for (int m = 0; m < in_dim; m++) {
    result += in_indices[m] * in_strides[m];
  }
  return result;
}

/**
 * @brief Update indices from base index
 *
 * @param in_dim
 * @param in_strides
 * @param in_index
 * @param in_indices
 */
void CellPairFinder::cellindices(int in_dim, const vulong &in_strides,
                                 const ulong &in_index, vulong &in_indices) const {
  size_t offset = 0;
  for (int m = 0; m < in_dim; m++) {
    in_indices[m] = (in_index - offset) / in_strides[m];
    offset += in_indices[m] * in_strides[m];
  }
}

}}} // namespace Portage::Meshfree::Pairs

