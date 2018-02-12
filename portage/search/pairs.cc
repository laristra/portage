/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "pairs.hh"

#include <climits>

#ifdef HAVE_CMATH
#include <cmath>
using std::ceil;
using std::floor;
using std::abs;
#else
#include <math.h>
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>
#endif

#include <list>
using std::list;

#include <algorithm>
using std::min;
using std::max;

#include <limits>
using std::numeric_limits;

#include <vector>
using std::vector;

namespace Portage {
namespace Meshfree {
namespace Pairs {

/// find bounding box of all the y-boxes, variable h
vpile PairsMinMax(const vpile& y, const vpile &h) {
  vector<size_t> ny(y.size());  // sizes
  size_t dim = ny[0];
  vpile r(2, dim);

  // find mins and maxes
  double shift = 2.0 * numeric_limits<double>::epsilon();  // smallest increment
  for (size_t m = 0; m < dim; m++) {
    r[0][m] = (y[m] - 2. * h[m]).min();
    r[1][m] = (y[m] + 2. * h[m]).max();
  }
  r[0] -= shift;
  r[1] += shift;
  return r;
}

/// find bounding box of all the y-boxes, constant h
vpile PairsMinMax(const vpile& y, const pile &h) {
  vector<size_t> ny(y.size());  // sizes
  size_t dim = ny[0];
  vpile r(2, dim);

  // find mins and maxes
  double shift = 2.0 * numeric_limits<double>::epsilon();  // smallest increment
  for (size_t m = 0; m < dim; m++) {
    r[0][m] = y[m].min() - 2. * h[m];
    r[1][m] = y[m].max() + 2. * h[m];
  }
  r[0] -= shift;
  r[1] += shift;
  return r;
}

/// convert real coordinates scaled to unit box to integers in [0,ulong_max]
inline void PairsIntegize(const size_t &dim, const pile &value,
                          const vpile &minmax, const pile &delta,
                          const vulong &sizes, vulong &ivalue) {
  for (size_t m = 0; m < dim; m++) {
    ivalue[m] = max<ulong>(
        min<ulong>(
            static_cast<ulong>(floor(
                sizes[m] * (value[m] - minmax[0][m]) / delta[m])),
            sizes[m] - 1),
        0);
  }
}

// index from indices
inline ulong cellindex(const size_t &dim, const vulong &strides,
                       const vulong &indices) {
  ulong result = 0;
  for (size_t m = 0; m < dim; m++) {
    result += indices[m] * strides[m];
  }
  return result;
}

// index from indices
inline void cellindices(const size_t &dim, const vulong &strides,
                        const ulong &index, vulong &indices) {
  size_t offset = 0;
  for (size_t m = 0; m < dim; m++) {
    indices[m] = (index - offset) / strides[m];
    offset += indices[m] * strides[m];
  }
}

/// good old reliable method of cells in arbitrary dimensions. fastest.
/** pairs arranged by x-point: all pairs for given x are contiguous.
 perfect linear scaling with number of points, but not efficient
 for large variation in h and numerous isolated bunches in x or y
 results of a typical optimized timing run for this method: (Mathematica input)

 (* 2D: equal size x and y boxes *)
 (* finding containments for 64 per side, 4096 points.*)
 (* finding containments for 128 per side, 16384 points.*)
 (* finding containments for 256 per side, 65536 points.*)
 (* finding containments for 512 per side, 262144 points.*)
 times={{4096,0.05},{16384,0.17},{65536,0.68},{262144,2.96}};
 ltimes = N@Log/@times;
 ListPlot[ltimes, PlotJoined->True];
 Print[Fit[ltimes, {1, x}, {x}]];
 ptimes={{15896,0.05},{65414,0.17},{260597,0.68},{1.04724e+06,2.96}};
 pltimes = N@Log/@ptimes;
 ListPlot[pltimes, PlotJoined->True];
 Print[Fit[pltimes, {1, x}, {x}]];
 (* 3D: equal size x and y boxes *)
 (* finding containments for 8 per side, 512 points.*)
 (* finding containments for 16 per side, 4096 points.*)
 (* finding containments for 32 per side, 32768 points.*)
 times={{512,0.02},{4096,0.12},{32768,1}};
 ltimes = N@Log/@times;
 ListPlot[ltimes, PlotJoined->True];
 Print[Fit[ltimes, {1, x}, {x}]];
 ptimes={{4952,0.02},{34954,0.12},{277243,1}};
 pltimes = N@Log/@ptimes;
 ListPlot[pltimes, PlotJoined->True];
 Print[Fit[pltimes, {1, x}, {x}]];

 \param x contained points
 \param y containing points with attached box
 \param h box sizes, centered on y
 \param half_pairs if true, checks that &x==&y and only returns pairs (i,j) with i>=j

 */
CellPairFinder::CellPairFinder(
    const vpile &xin, const vpile &yin, const vpile &hin,
    const bool do_scatter_in)
    : x(xin), y(yin), h(hin), dim(x.size()[0]),
      do_scatter(do_scatter_in)
{
  // get sizes and check
  ulong nx=x.size()[1];
  ulong ny=y.size()[1];
  assert (y.size()[0] == dim);
  assert (h.size()[0] >= dim);
  if (do_scatter)
    assert (h.size()[1] == nx);
  else
    assert (h.size()[1] == ny);

  // find max h
  hmax.resize(dim);
  for (size_t m = 0; m < dim; m++) {
    hmax[m] = h[m].max();
  }

  // find min and max of enclosing box of y points
  if (do_scatter)
    yminmax = PairsMinMax(y,hmax);
  else
    yminmax = PairsMinMax(y,h);
  delta = yminmax[1]-yminmax[0];

  // decide on size of grid - this is a maximum value
  ulong maxmemory=nx/10;
  size_t nsidemax=max<size_t>(static_cast<size_t>(ceil(pow(maxmemory*1., (1./dim)))),1);

  // find avg of min and max h
  pile havg(dim);
  havg = 0.;
  for (size_t m=0; m<dim; m++) for (size_t i=0; i<ny; i++) {
    havg[m] += fabs(h[m][i]);
  }
  havg = 4.*havg/(1.*ny) + numeric_limits<double>::epsilon();
  vulong nsideh(dim);
  for (size_t m=0; m<dim; m++) nsideh[m] = static_cast<size_t>(ceil(delta[m]/havg[m]));

  // set number of cells each coordinate direction
  vulong nsides(dim);
  for (size_t m=0;m<dim;m++) nsides[m] = min<ulong>(nsidemax, nsideh[m]);
  size_t ncells=1;
  for (size_t m=0;m<dim;m++) ncells*=nsides[m];
  nsidesm = nsides-static_cast<ulong>(1);

  // allocate grid of lists
  cells.resize(ncells);

  // strides
  strides.resize(dim);
  strides[dim-1] = 1;
  for(int m=dim-2; m>=0; m--) {
    strides[m] = strides[m+1]*nsides[m+1];
  }

  // hash x-points
  for (ulong i=0; i<nx; i++) {
    // ignore x values outside bounding box
    bool outside=false;
    pile xi(dim);
    for (size_t m=0;m<dim;m++) {
      if (x[m][i]<=yminmax[0][m]) outside=true;
      if (x[m][i]>=yminmax[1][m]) outside=true;
      xi[m]=x[m][i];
    }
    if (outside) continue;

    // get cell indices this x
    vulong ix(dim);
    PairsIntegize(dim, xi, yminmax, delta, nsidesm, ix);

    // add this x to cell list
    ulong index=cellindex(dim, strides, ix);
    cells[index].push_back(i);
  }

}  // CellPairFinder::CellPairFinder


/// get pairs for target point j, gather case
list<ulong> CellPairFinder::find_gather(const ulong j) const
{
  list<ulong> pairlist;

  // get cell indices lower left and upper right corners of box
  pile yll(dim), yur(dim);
  for (size_t m=0;m<dim;m++) {
    yll[m]=y[m][j]-2.*h[m][j];
    yur[m]=y[m][j]+2.*h[m][j];
  }
  vulong iyl(dim), iyu(dim);
  PairsIntegize(dim, yll, yminmax, delta, nsidesm, iyl);
  PairsIntegize(dim, yur, yminmax, delta, nsidesm, iyu);

  // total number of cells for this y
  ulong ncellsy=iyu[dim-1]-iyl[dim-1]+1;
  vulong ystrides(dim);
  ystrides[dim-1] = 1;
  for(int m=dim-2; m>=0; m--) {
    ystrides[m] = ncellsy;
    ncellsy*=iyu[m]-iyl[m]+1;
  }

  // scan cells for this y
  for (ulong cell=0; cell<ncellsy; cell++) {
    // convert local y-cell indices to global cell index
    vulong yndx(dim);
    cellindices(dim, ystrides, cell, yndx);
    vulong cellis(dim);
    for(size_t m=0; m<dim; m++) cellis[m] = iyl[m]+yndx[m];
    size_t celli = cellindex(dim, strides, cellis);

    // determine if in interior or boundary of y-cell
    bool ybndry=false;
    for (size_t m=0; m<dim; m++)
      if (cellis[m]==iyl[m] || cellis[m]==iyu[m]) ybndry=true;

    // loop over all x's in this cell's list
    for (vector<ulong>::const_iterator i=cells[celli].begin();
        i != cells[celli].end(); i++) {
      // if on y-cell boundary, check that x's are contained
      bool inside = true;
      if (ybndry) {
        for(size_t m=0; m<dim; m++) {
          if (x[m][*i] <= yll[m]) inside = false;
          if (x[m][*i] >= yur[m]) inside = false;
        }
      }

      // add pair: put x's in this y-cell onto neighbor list, if inside
      if (inside) {
        pairlist.push_back(*i);
      }
    }  // for i
  }  // for cell

  return pairlist;
}  // CellPairFinder::find_gather


/// get pairs for target point j, scatter case
list<ulong> CellPairFinder::find_scatter(const ulong j) const
{
  list<ulong> pairlist;

  // get cell indices lower left and upper right corners of box
  pile yllmax(dim), yurmax(dim);
  for (size_t m=0;m<dim;m++) {
    yllmax[m]=y[m][j]-2.*hmax[m];
    yurmax[m]=y[m][j]+2.*hmax[m];
  }
  vulong iyl(dim), iyu(dim);
  PairsIntegize(dim, yllmax, yminmax, delta, nsidesm, iyl);
  PairsIntegize(dim, yurmax, yminmax, delta, nsidesm, iyu);

  // total number of cells for this y
  ulong ncellsy=iyu[dim-1]-iyl[dim-1]+1;
  vulong ystrides(dim);
  ystrides[dim-1] = 1;
  for(int m=dim-2; m>=0; m--) {
    ystrides[m] = ncellsy;
    ncellsy*=iyu[m]-iyl[m]+1;
  }

  // scan cells for this y
  for (ulong cell=0; cell<ncellsy; cell++) {
    // convert local y-cell indices to global cell index
    vulong yndx(dim);
    cellindices(dim, ystrides, cell, yndx);
    vulong cellis(dim);
    for(size_t m=0; m<dim; m++) cellis[m] = iyl[m]+yndx[m];
    size_t celli = cellindex(dim, strides, cellis);

    // TODO:  Try adding a boundary test, similar to what's in
    // find_gather, based on hmin instead of individual h's.
    // Do profiling to see whether that speeds up the search
    // in typical larger cases.

    // loop over all x's in this cell's list
    for (vector<ulong>::const_iterator i=cells[celli].begin();
        i != cells[celli].end(); i++) {
      pile yll(dim), yur(dim);
      for (size_t m=0;m<dim;m++) {
        yll[m]=y[m][j]-2.*h[m][*i];
        yur[m]=y[m][j]+2.*h[m][*i];
      }
      // check that x's are contained
      bool inside = true;
      for(size_t m=0; m<dim; m++) {
        if (x[m][*i] <= yll[m]) inside = false;
        if (x[m][*i] >= yur[m]) inside = false;
      }

      // add pair: put x's in this y-cell onto neighbor list, if inside
      if (inside) {
        pairlist.push_back(*i);
      }
    }  // for i
  }  // for cell

  return pairlist;
}  // CellPairFinder::find_scatter


/// get pairs for target point j
list<ulong> CellPairFinder::find(const ulong j) const
{
  if (do_scatter)
    return(find_scatter(j));
  else  // gather
    return(find_gather(j));
}  // CellPairFinder::find


#if 0
/// neighbor finding using sort
/** order is approximately n^1.6 in 2D and n^1.8 in 3D, n=number of points
 number of open boxes grows as n^.5, time for sort as n^1.2
 this method feels like it should be faster - but isn't because of "straddlers",
 small boxes that straddle large partitions.
 more appropriate for wide ranges of h sizes and numerous isolated blobs
 results of a typical optimized timing run for this method: (Mathematica input)

 (* 2D: equal size x and y boxes *)
 (* finding containments for 64 per side, 4096 points.*)
 (* finding containments for 128 per side, 16384 points.*)
 (* finding containments for 256 per side, 65536 points.*)
 (* finding containments for 512 per side, 262144 points.*)
 times={{4096,0.09},{16384,0.72},{65536,6.47},{262144,73.88}};
 ltimes = N@Log/@times;
 ListPlot[ltimes, PlotJoined->True];
 Print[Fit[ltimes, {1, x}, {x}]];
 ptimes={{15896,0.09},{65414,0.72},{260597,6.47},{1.04724e+06,73.88}};
 pltimes = N@Log/@ptimes;
 ListPlot[pltimes, PlotJoined->True];
 Print[Fit[pltimes, {1, x}, {x}]];

 (* 3D: equal size x and y boxes *)
 (* finding containments for 8 per side, 512 points.*)
 (* finding containments for 16 per side, 4096 points.*)
 (* finding containments for 32 per side, 32768 points.*)
 times={{512,0.02},{4096,0.46},{32768,32.65}};
 ltimes = N@Log/@times;
 ListPlot[ltimes, PlotJoined->True];
 Print[Fit[ltimes, {1, x}, {x}]];
 ptimes={{4952,0.02},{34954,0.46},{277243,32.65}};
 pltimes = N@Log/@ptimes;
 ListPlot[pltimes, PlotJoined->True];
 Print[Fit[pltimes, {1, x}, {x}]];
 */

shared_ptr<vector<list<ulong>>> PairsContainSort(const vpile &x, const vpile &y, const vpile &h)
{
  // get sizes
  size_t dim=x.size()[0];
  ulong nx=x.size()[1];
  ulong ny=y.size()[1];

  // upper right and lower left corners of y-boxes
  vpile yl(dim,ny), yu(dim,ny);
  for (size_t m=0;m<dim;m++) {
    yl[m] = y[m]-2.*h[m];
    yu[m] = y[m]+2.*h[m];
  }

  // find min and max of enclosing box of y points
  vpile yminmax(PairsMinMax(y,h));
  pile delta(yminmax[1]-yminmax[0]);

  // get deepest possible level
  size_t nlevels = nbits_ulong/dim;

  //declare vector of data to be sorted
  vector<vulong> points(nx+2*ny, vulong(3));

#ifdef TIME_PAIRS
  clock_t t1, t2;
  double t;
  t1 = clock();
#endif

  // for each y box truncate bits of l.l. at deepest level that contains l.l. and u.r.
  for (size_t j=0;j<ny;j++) {
    // useful declarations
    pile yj(dim);
    vulong iy(dim);

    // get interleaved coordinate bits this yl
    for (size_t m=0;m<dim;m++) yj[m]=yl[m][j];
    PairsIntegize(dim, yj, yminmax, delta, iy);
    bitset<nbits_ulong> yllbits(0);
    for (size_t m=0;m<dim;m++) {
      bitset<nbits_ulong> by(iy[m]);
      for (size_t k=0; k<nlevels; k++) {
        size_t index = nbits_ulong - (k+1)*dim + m;
        yllbits[index] = by[nbits_ulong-k-1];
      }
    }

    // get interleaved coordinate bits this yu
    for (size_t m=0;m<dim;m++) yj[m]=yu[m][j];
    PairsIntegize(dim, yj, yminmax, delta, iy);
    bitset<nbits_ulong> yurbits(0);
    for (size_t m=0;m<dim;m++) {
      bitset<nbits_ulong> by(iy[m]);
      for (size_t k=0; k<nlevels; k++) {
        size_t index = nbits_ulong - (k+1)*dim + m;
        yurbits[index] = by[nbits_ulong-k-1];
      }
    }

    // put the data into the sort vector
    vulong pointval(3);
#define MINIMAL_OPEN
#ifdef MINIMAL_OPEN
    pointval[0] = yllbits.to_ulong();
    pointval[1] = j;
    pointval[2] = 0;
    points[2*j] = pointval;

    pointval[0] = yurbits.to_ulong();
    pointval[1] = j;
    pointval[2] = 1;
    points[2*j+1] = pointval;
#else 
    // scan down bits looking for first difference between l.l. and u.r.
    size_t klevel=nlevels;
    bitset<nbits_ulong> ylltrunc(0), yurtrunc(0);
    yurtrunc.set();
    for (size_t k=0;k<nlevels;k++) {
      bool bail=false;
      for (size_t m=0;m<dim;m++) {
        size_t index = nbits_ulong - (k+1)*dim + m;
        if (yllbits[index] != yurbits[index]) {
          klevel = k;
          bail=true;
          break;
        } else {
          ylltrunc[index] = yllbits[index];
          yurtrunc[index] = yurbits[index];
        }
      }
      if (bail) break;
    }

    // make sure all bits for last level are set correctly
    // l.l. has all zeros and u.r. has all ones on right end of bitset
    for (size_t m=0;m<dim;m++) {
      size_t index = nbits_ulong - (klevel+1)*dim + m;
      ylltrunc[index] = 0;
      yurtrunc.set(index);
    }

    pointval[0] = ylltrunc.to_ulong();
    pointval[1] = j;
    pointval[2] = 0;
    points[2*j] = pointval;

    pointval[0] = yurtrunc.to_ulong();
    pointval[1] = j;
    pointval[2] = 1;
    points[2*j+1] = pointval;
#endif

#ifdef DUMP_EM_Y
    cout << j << " y="; for(size_t m=0;m<dim;m++) cout << y[m][j] << "+-"<< 2.*h[m][j] << "\n";
    for(size_t k=0;k<nbits_ulong;k++) cout << yllbits [nbits_ulong-k-1]; cout << "\n";
    for(size_t k=0;k<nbits_ulong;k++) cout << yurbits [nbits_ulong-k-1]; cout << "\n";
#ifndef MINIMAL_OPEN
    for(size_t k=0;k<nbits_ulong;k++) cout << ylltrunc[nbits_ulong-k-1]; cout << "\n";
    for(size_t k=0;k<nbits_ulong;k++) cout << yurtrunc[nbits_ulong-k-1]; cout << "\n";
#endif
    cout << "\n";
#endif
  }

#ifdef TIME_PAIRS
  t2 = clock();
  t=(double(t2-t1)/CLOCKS_PER_SEC);
  cout << "PairsContain<"<<dim<<">:"<< t <<" sec for y-bits\n";
#endif

  // number of y's and x's that are in bounding box
  size_t pointcount=2*ny;

#ifdef TIME_PAIRS
  t1 = clock();
#endif

  // get x bits
  for (ulong i=0; i<nx; i++) {
    bool outside=false;
    pile xi(dim);
    for (size_t m=0;m<dim;m++) {
      if (x[m][i]<=yminmax[0][m]) outside=true;
      if (x[m][i]>=yminmax[1][m]) outside=true;
      xi[m]=x[m][i];
    }
    if (outside) continue;

    // get coordinate bits this x
    vulong ix(dim);
    PairsIntegize(dim, xi, yminmax, delta, ix);
    bitset<nbits_ulong> xbits(0);
    for (size_t m=0;m<dim;m++) {
      bitset<nbits_ulong> bx(ix[m]);
      for (size_t k=0; k<nlevels; k++) {
        size_t index = nbits_ulong - (k+1)*dim + m;
        xbits[index] = bx[nbits_ulong-k-1];
      }
    }

    // put the data into the sort vector
    vulong pointval(3);
    pointval[0] = xbits.to_ulong();
    pointval[1] = i;
    pointval[2] = 2;
    points[pointcount] = pointval;

    // bump counter of total number of points
    pointcount++;

#ifdef DUMP_EM_X
    cout << i << " x="; for(size_t m=0;m<dim;m++) cout << x[m][i] << "\n";
    for(size_t k=0;k<nbits_ulong;k++) cout << xbits [nbits_ulong-k-1]; cout << "\n";
    cout << "\n";
#endif
  }

#ifdef TIME_PAIRS
  t2 = clock();
  t=(double(t2-t1)/CLOCKS_PER_SEC);
  cout << "PairsContain<"<<dim<<">:"<< t <<" sec for x-bits\n";
#endif

  // remove unneeded points to shorten sort
  points.resize(pointcount);

#ifdef TIME_PAIRS
  t1 = clock();
#endif

  // sort list of xy points
  less_vulong0 less;
  vector<vulong>::iterator first=points.begin(), last=points.end();
  sort<vector<vulong>::iterator, less_vulong0>(first, last, less);

#ifdef TIME_PAIRS
  t2 = clock();
  t=(double(t2-t1)/CLOCKS_PER_SEC);
  cout << "PairsContain<"<<dim<<">:"<< t <<" sec for sort\n";
#endif

#ifdef DEBUG_PAIRS
  cout << "Sorted points:\n";
  for (size_t i=0; i<points.size(); i++) {
    cout << points[i][0] <<" "<< points[i][1] << " " << points[i][2] << "\n";
  }
#endif    

  // list of open y-boxes
  list<ulong> openy;

  // list of pairs
  auto pairlist_p = make_shared<vector<list<ulong>>>(nx);
  vector<list<ulong>> &pairlist(*pairlist_p);

  // vector of iterators for boxes on open list
  vector< list<ulong>::iterator > openyiter(ny,openy.begin());

#ifdef TIME_PAIRS
  t1 = clock();
  ulong nopenmax=0, nopenavg=0, nopen;
  ulong nreject, nrejmax=0, nrejavg=0;
#endif

  // scan sorted list, adding pairs as needed
  for (size_t i =0; i<points.size(); i++) {
    // add a new y-box to the open list
    if (points[i][2]==0) {
      size_t index = points[i][1];
      list<ulong>::iterator end=openy.end();
      list<ulong>::iterator back;
      back = openy.insert(end, index);
      openyiter[index] = back;
    }

    // remove an old y-box from the open list
    if (points[i][2]==1) {
      openy.erase( openyiter[ points[i][1] ] );
    }

    // add pair for this x and all open y's
    nreject=0;
    if (points[i][2]==2 && openy.size()>0) {
      size_t ix = points[i][1];
      for (list<ulong>::iterator ybox=openy.begin(); ybox!=openy.end(); ybox++) {
        ulong jy = (*ybox);
        bool inside = true;
        for (size_t m=0;m<dim;m++) {
          if (x[m][ix]<yl[m][jy]) inside = false;
          if (x[m][ix]>yu[m][jy]) inside = false;
        }
        if (inside) {
          pairlist[ix].push_back(jy);
#ifdef TIME_PAIRS
        } else {
          nreject++;
#endif
        }
      }

#ifdef TIME_PAIRS
      nopen = openy.size();
      nopenmax = max<ulong>(nopenmax, nopen);
      nopenavg += nopen;
      nrejmax = max<ulong>(nrejmax, nreject);
      nrejavg += nreject;
#endif
    }

#ifdef DEBUG_PAIRS
    if (openy.size()>0) {
      cout << "openy(" << i << ")=";
      int openycount=0;
      for (list<ulong>::iterator j=openy.begin();j!=openy.end();j++) {
        cout << *j << " ";
        openycount++;
        if (openycount == 30) {
          cout << "\n";
          openycount = 0;
        }
      }
      cout << "\n";
    }
#endif
  }

#ifdef TIME_PAIRS
  cout << "PairsContain<"<<dim<<">:" << " open_max = " << nopenmax << "\n";
  cout << "PairsContain<"<<dim<<">:" << " open_average = " << (1.*nopenavg)/nx << "\n";
  cout << "PairsContain<"<<dim<<">:" << " reject_max = " << nrejmax << "\n";
  cout << "PairsContain<"<<dim<<">:" << " reject_average = " << (1.*nrejavg)/nx << "\n";
  t2 = clock();
  t=(double(t2-t1)/CLOCKS_PER_SEC);
  cout << "PairsContain<"<<dim<<">:"<< t <<" sec for scan\n";
#endif

#ifdef DEBUG_PAIRS
  cout <<"Pairs: \n";
  int paircount=0;
  for (size_t i=0; i<nx; i++) {
    for (auto j=pairlist[i].begin(); j!=pairslist[i].end(); j++) {
      cout << "(" << i << "," << *j << ") ";
      paircount++;
      if (paircount >= 20) {
        cout <<"\n";
        paircount=0;
      }
    }
    cout << "\n";
#endif

    return pairlist_p;
  }

  /// Neighbor finding using hash tables for y-boxes. Slow.
  /** Loosely based on Warren & Salmon's 1993 paper. 
   Equivalent to method of cells, only fancier.
   This method shows perfect linear scaling with number of points.
   It is most suitable for x and y densely packed and mostly uniformly distributed.
   Results of a typical optimized timing run for this method: (Mathematica input)

   (* 2D: equal size x and y boxes *)
   (* finding containments for 64 per side, 4096 points.*)
   (* finding containments for 128 per side, 16384 points.*)
   (* finding containments for 256 per side, 65536 points.*)
   (* finding containments for 512 per side, 262144 points.*)
   times={{4096,0.13},{16384,0.52},{65536,2.02},{262144,8.4}};
   ltimes = N@Log/@times;
   ListPlot[ltimes, PlotJoined->True];
   Print[Fit[ltimes, {1, x}, {x}]];
   ptimes={{15896,0.13},{65414,0.52},{260597,2.02},{1.04724e+06,8.4}};
   pltimes = N@Log/@ptimes;
   ListPlot[pltimes, PlotJoined->True];
   Print[Fit[pltimes, {1, x}, {x}]];

   (* 3D: equal size x and y boxes *)
   (* finding containments for 8 per side, 512 points.*)
   (* finding containments for 16 per side, 4096 points.*)
   (* finding containments for 32 per side, 32768 points.*)
   times={{512,0.08},{4096,0.63},{32768,4.45}};
   ltimes = N@Log/@times;
   ListPlot[ltimes, PlotJoined->True];
   Print[Fit[ltimes, {1, x}, {x}]];
   ptimes={{4952,0.08},{34954,0.63},{277243,4.45}};
   pltimes = N@Log/@ptimes;
   ListPlot[pltimes, PlotJoined->True];
   Print[Fit[pltimes, {1, x}, {x}]];
   */
template<size_t ndim>
shared_ptr<vector<list<ulong>>> PairsContainHashY(const vpile &x, const vpile &y, const vpile &h)
{
  // get sizes
  size_t dim=x.size()[0], nboxes=y.size()[1];
  assert(ndim==dim);
  size_t nverts=1; for(size_t m=0;m<dim;m++) nverts*=2;
  ulong nx=x.size()[1];

  // total number of bits in coordinates
  const size_t nbits=ndim*nbits_ulong;

  // find min and max of enclosing box of points
  vpile minmax(PairsMinMax(y,h));
  pile delta(minmax[1]-minmax[0]);

  // find levels for y boxes and optimal number of hash levels
  vulong ylevels(nboxes);
  for (ulong j=0;j<nboxes;j++) {
    // integize h
    pile hj(dim); for (size_t m=0;m<dim;m++) hj[m]=4.*h[m][j];
    vulong ih(dim);
    PairsIntegize(dim, hj, minmax, delta, ih);

    // find level at which this box lives
    valarray<int> k1(0,dim), count(0,dim);
    for (size_t m=0;m<dim;m++) {
      bitset<nbits_ulong> hbits(ih[m]);
      for (int k=nbits_ulong-1;k>=0;k--) {
        if (hbits[k] && k1[m]==0) k1[m]=k+1;
        if (hbits[k]) count[m]++;
      }
      if (count[m]==0 & k1[m]>0) k1[m]--;
      if (k1[m]<0) k1[m]=0;
    }
    size_t k=k1.max();
    ylevels[j]=min(max(static_cast<int>(nbits_ulong-1-k), 0),
        static_cast<int>(nbits_ulong-1));
  }

  // maximum level for all y
  ulong maxylevel = ylevels.max();

  // log2(max # vertices) - adjust to fit into memory
  // max size of hash tables is 2^(maxlog2ny)
  ulong maxlog2ny = 21;

  // max level allowed to fit in memory
  ulong maxlevel = maxlog2ny/dim-1;

  // select optimum number of levels
  ulong nlevels = min<ulong>(maxylevel, maxlevel)+1;
  assert(nlevels*dim<=nbits_ulong);

  // number of bits for hash table for each level
  vector<ulong> hashbits(nlevels);
  for (size_t i=0;i<nlevels;i++) hashbits[i]=(i+1)*dim;

  // size of hash table for each level
  ulong nvertsdim=1; for(size_t m=0;m<dim;m++) nvertsdim*=2;
  vector<ulong> hashsize(nlevels);
  for(size_t level=0;level<nlevels;level++) {
    hashsize[level]=1;
    for(size_t i=0;i<=level;i++) hashsize[level]*=nvertsdim;
  }

  // hash mask for each level
  vector< bitset<nbits> > hashmask(nlevels, 0);
  for (size_t level=0;level<nlevels;level++) {
    for (size_t i=0;i<hashbits[level];i++) hashmask[level].set(nbits-i-1,true);
  }

  // hash bits to shift for each level
  vector<ulong> hashshift(nlevels);
  for (ulong level=0; level<nlevels; level++) {
    hashshift[level]=nbits-hashbits[level];
  }

  // allocate hash table for each level
  typedef list<ulong> hashlist_t;
  vector< vector< hashlist_t > > hashtable(nlevels);
  for (size_t level=0;level<nlevels;level++) {
    hashtable[level]=vector<hashlist_t>(hashsize[level]);
  }

  // fill hash table with partitions for each vertex of each y-box
  for (ulong j=0;j<nboxes;j++) {
    ulong level = ylevels[j];
    for (size_t i=0;i<nverts;i++) {
      // find vertex coordinates
      bitset<nbits_ulong> mbits(i);
      pile vert(dim);
      for(size_t m=0;m<dim;m++) {
        vert[m] = y[m][j] + 2.*(mbits[m]-.5)*2.*h[m][j];
      }
      // integize vertex coordinates
      vulong ivert(dim);
      PairsIntegize(dim, vert, minmax, delta, ivert);
      // interleave bits of coordinate integers
      bitset<nbits> vbits(0);
      for (size_t m=0;m<dim;m++) {
        bitset<nbits_ulong> bvert(ivert[m]);
        for (size_t k=0; k<nbits_ulong; k++) {
          vbits[k*dim+m] = bvert[k];
        }
#ifdef DUMP_EM_Y
        cout << "y=" <<j<<" "<<i<<" "<<m<<" "<<vert[m]<<" ";
        for(size_t k=0;k<nbits_ulong;k++) cout <<vbits[nbits_ulong-1-k];
        cout <<"\n";
#endif	    
      }
#ifdef DUMP_EM_Y
      cout <<"\n";
#endif	    
      // push box onto hash table
      bitset<nbits> partition = hashmask[level] & vbits;
      ulong hash = (partition>>hashshift[level]).to_ulong();
      hashtable[level][hash].push_back(j);
    }
  }

#ifdef DUMP_EM_TSIZES
  for (ulong level=0;level<nlevels;level++) {
    cout <<"--------------------------------------\n";
    ulong size=0;
    for (ulong hash=0; hash<hashsize[level]; hash++) {
      if (!hashtable[level][hash].empty()) size += hashtable[level][hash].size();
    }
    cout <<"level "<<level<<" hash size: total="<<size<<" avg= "<<size/hashsize[level]<<"\n";
  }
#endif

#ifdef DUMP_EM_SIZES
  for (ulong level=0;level<nlevels;level++) {
    cout <<"--------------------------------------\n";
    cout <<"level "<< level << ":\n";

    cout <<"hash table sizes: \n";
    size_t size_count=0;
    for (ulong hash=0; hash<hashsize[level]; hash++) {
      size_count++;
      if (!hashtable[level][hash].empty()) cout << hashtable[level][hash].size() << " ";
      else cout << 0 << " ";
      if (size_count==50) {size_count=0; cout <<"\n";}
    }
    if (size_count!=0) cout <<"\n";
  }
#endif

#ifdef DUMP_EM_TABLE
  for (ulong level=0;level<nlevels;level++) {
    cout <<"--------------------------------------\n";
    cout <<"level "<< level << ":\n";
    cout <<"hash table: \n";
    for (ulong hash=0; hash<hashsize[level]; hash++) {
      if (!hashtable[level][hash].empty()) {
        for (hashlist_t::iterator hashtest=hashtable[level][hash].begin();
            hashtest!=hashtable[level][hash].end();
            hashtest++)
        {
          cout << hex << hash << " " << dec<< " " << *hashtest << "\n";
        }
      }
    }
  }
#endif

  // allocate box counting vector
  vulong box_count(nboxes);

  // allocate list of pairs
  auto pairlist_p = make_shared<vector<list<ulong>>>(nx);
  vector<list<ulong>> &pairlist(*pairlist_p);

  // compare partitions this x against hash table entries for pairs
  for (ulong i=0; i<nx; i++) {

    // get coordinates this x and check inside y-box
    bool outside=false;
    pile xi(dim);
    for (size_t m=0;m<dim;m++) {
      if (x[m][i]<=minmax[0][m]) outside=true;
      if (x[m][i]>=minmax[1][m]) outside=true;
      xi[m]=x[m][i];
    }
    if (outside) continue;

    // get coordinate bits this x
    vulong ix(dim);
    PairsIntegize(dim, xi, minmax, delta, ix);
    bitset<nbits> xbits(0);
    for (size_t m=0;m<dim;m++) {
      bitset<nbits_ulong> bx(ix[m]);
      for (size_t k=0; k<nbits_ulong; k++) {
        xbits[k*dim+m] = bx[k];
      }
    }

    // scan levels this x
    for (size_t level=0;level<nlevels;level++) {

      // get partition and hash
      bitset<nbits> partition = hashmask[level] & xbits;
      ulong hash = (partition>>hashshift[level]).to_ulong();

      // reset box counters
      for( hashlist_t::iterator hashbox=hashtable[level][hash].begin();
          hashbox!=hashtable[level][hash].end(); hashbox++)
      {
        box_count[*hashbox]=0;
      }

      // scan boxes for strict containments, store on list
      for( hashlist_t::iterator hashbox=hashtable[level][hash].begin();
          hashbox!=hashtable[level][hash].end(); hashbox++)
      {
        // mark visited, add pair if needed, skip if already done
        ulong box=*hashbox;
        if (box_count[box]>0) continue;
        box_count[box]++;
        bool contained=true;
        double dist, adist;
        for(size_t m=0;m<dim;m++) {
          dist = x[m][i]-y[m][box];
          adist = fabs(dist);
          if (adist >= 2.0*h[m][box]) contained=false;
        }
#ifdef DUMP_EM_X
        cout <<"x="; for(size_t m=0;m<dim;m++) cout<<x[m][i]<<" ";
        cout <<" "; for(size_t k=0;k<nbits;k++) cout << xbits[nbits-k-1]; cout <<" ";
        cout<<"level "<<level<<" " << hash <<" part "<<hex<< partition<<
        dec << " box "<<box;
#endif	  
        if (contained) {
          pairlist[i].push_back(box);
#ifdef DUMP_EM_X
          cout << " contained";
#endif	 
        }
#ifdef DUMP_EM_X
        cout <<"\n";
#endif	  
      }
    }
  }

  return pairlist_p;
}

/// neighbor finding driver routine
shared_ptr<pairs_data_t> PairsFind(
    const vpile &x, const vpile &y, const vpile &h,
    const bool do_scatter,
    const contain_type type, const bool half_pairs)
{
  // get and check sizes
  vector<size_t> dimx = x.size(), dimy=y.size(), dimh=h.size();
  assert(dimx[0]==dimy[0] && dimy[0]<=dimh[0]);
  assert(dimy[1]==dimh[1]);
  size_t dim=dimx[0];

  // check other arguments
  if (half_pairs) assert(type == CELLS);

  // compute containments
  switch (type) {
    case CELLS:
    return PairsContainCells(x,y,h,half_pairs);
    case SORT:
    return PairsContainSort(x,y,h);
    case HASHY:
    switch (dim) {
      case 1:
      return PairsContainHashY<1>(x,y,h);
      case 2:
      return PairsContainHashY<2>(x,y,h);
      case 3:
      return PairsContainHashY<3>(x,y,h);
      case 4:
      return PairsContainHashY<4>(x,y,h);
      case 5:
      return PairsContainHashY<5>(x,y,h);
      case 6:
      return PairsContainHashY<6>(x,y,h);
      case 7:
      return PairsContainHashY<7>(x,y,h);
      default:
      assert(false);  // nothing else allowed
    }
    default:
    assert(false);  // not allowed
  }

  return make_shared<pairs_data_t>();
}
#endif

} // namespace Pairs
} // namespace Meshfree
} // namespace Portage

