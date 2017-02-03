/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



#ifndef SRC_SUPPORT_PORTAGE_H_
#define SRC_SUPPORT_PORTAGE_H_

#ifdef THRUST

#include "thrust/device_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

#else  // no thrust

#include <boost/iterator/counting_iterator.hpp>

#include <vector>
#include <algorithm>

#endif

/*
  @file portage.h
  @brief Several utility types and functions within the Portage namespace.
 */

/*
  @namespace Portage
  The Portage namespace houses all of the code within Portage.

  Cells (aka zones/elements) are the highest dimension entities in a mesh
  Nodes (aka vertices) are lowest dimension entities in a mesh
  Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
  BOUNDARY_FACE is a special type of entity that is need so that process
  kernels can define composite vectors (see src/data_structures) on
  exterior boundary faces of the mesh only

  Wedges are special subcell entities that are a simplicial
  decomposition of cell. In 3D, a wedge is tetrahedron formed by one
  point of the edge, the midpoint of the edge, the "center" of the
  face and the "center" of the cell volume. In 2D, a wedge is a
  triangle formed by an end-point of the edge, the mid-point of the
  edge and the center of the cell. In 1D, wedges are lines, that are
  formed by the endpoint of the cell and the midpoint of the
  cell. There are two wedges associated with an edge of cell face in
  3D.

  Corners are also subcell entities that are associated uniquely with
  a node of a cell. Each corner is the union of all the wedges incident
  upon that node in the cell

  Facets are the boundary entity between two wedges in adjacent
  cells. In 3D, a facet is a triangular subface of the cell face
  shared by two wedges in adjacent cells. In 2D, a facet is half of
  an edge that is shared by two wedges in adjacent cells
 */
namespace Portage {
// TODO:  Right now we're relying on the fact that this enum is
//        identical to Jali::Entity_kind.  Need to fix this.
/// The type of mesh entity.
enum Entity_kind {
  ALL_KIND = -3,     /*!< All possible types */
  ANY_KIND = -2,     /*!< Any of the possible types */
  UNKNOWN_KIND = -1, /*!< Usually indicates an error */
  NODE = 0,
  EDGE,
  FACE,
  CELL,
  SIDE,
  WEDGE,
  CORNER,
  FACET,
  BOUNDARY_FACE
};

const int NUM_ENTITY_KINDS = 8;

// Parallel status of entity
/// The parallel type of a given entity.
enum Entity_type {
  TYPE_UNKNOWN = -1,
  DELETED = 0,
  PARALLEL_OWNED = 1,   /*!< Owned by this processor */
  PARALLEL_GHOST = 2,   /*!< Owned by another processor */
  BOUNDARY_GHOST = 3,   /*!< Ghost/Virtual entity on boundary */
  ALL  = 4              /*!< PARALLEL_OWNED + PARALLEL_GHOST + BOUNDARY_GHOST */
};

/// Element (cell topology) type
enum Element_type {
  UNKNOWN_TOPOLOGY = 0,
  TRI,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHEDRON
};

/// Limiter type
typedef enum {NOLIMITER, VAN_LEER, BARTH_JESPERSEN, MINMOD, SUPERBEE}
  LimiterType;


#ifdef THRUST

template<typename T>
    using vector = thrust::device_vector<T>;

template<typename T>
    using pointer = thrust::device_ptr<T>;

typedef thrust::counting_iterator<int> counting_iterator;
inline counting_iterator make_counting_iterator(int const i) {
  return thrust::make_counting_iterator(i);
}

template<typename InputIterator, typename OutputIterator,
         typename UnaryFunction>
inline OutputIterator transform(InputIterator first, InputIterator last,
                                OutputIterator result, UnaryFunction op) {
  return thrust::transform(first, last, result, op);
}

template<typename InputIterator1, typename InputIterator2,
         typename OutputIterator, typename BinaryFunction>
inline OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                                InputIterator2 first2, OutputIterator result,
                                BinaryFunction op) {
  return thrust::transform(first1, last1, first2, result, op);
}

template<typename InputIterator, typename UnaryFunction>
inline void for_each(InputIterator first, InputIterator last,
                              UnaryFunction f) {
  thrust::for_each(first, last, f);
}

#else  // no thrust

template<typename T>
    using vector = std::vector<T>;

template<typename T>
    using pointer = T*;

typedef boost::counting_iterator<int> counting_iterator;
inline counting_iterator make_counting_iterator(int const i) {
  return boost::make_counting_iterator<int>(i);
}

template<typename InputIterator, typename OutputIterator,
    typename UnaryFunction>
inline OutputIterator transform(InputIterator first, InputIterator last,
                                OutputIterator result, UnaryFunction op) {
  return std::transform(first, last, result, op);
}

template<typename InputIterator1, typename InputIterator2,
         typename OutputIterator, typename BinaryFunction>
inline OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                                InputIterator2 first2, OutputIterator result,
                                BinaryFunction op) {
  return std::transform(first1, last1, first2, result, op);
}

template<typename InputIterator, typename UnaryFunction>
inline void for_each(InputIterator first, InputIterator last,
                     UnaryFunction f) {
  std::for_each(first, last, f);
}


#endif

struct Weights_t {
  int entityID;
  std::vector<double> weights;
};

}  // namespace Portage

#endif  // SRC_SUPPORT_PORTAGE_H_
