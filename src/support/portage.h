/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

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

namespace Portage {

// Cells (aka zones/elements) are the highest dimension entities in a mesh
// Nodes (aka vertices) are lowest dimension entities in a mesh
// Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
// BOUNDARY_FACE is a special type of entity that is need so that process
// kernels can define composite vectors (see src/data_structures) on
// exterior boundary faces of the mesh only
//
// Wedges are special subcell entities that are a simplicial
// decomposition of cell. In 3D, a wedge is tetrahedron formed by one
// point of the edge, the midpoint of the edge, the "center" of the
// face and the "center" of the cell volume. In 2D, a wedge is a
// triangle formed by an end-point of the edge, the mid-point of the
// edge and the center of the cell. In 1D, wedges are lines, that are
// formed by the endpoint of the cell and the midpoint of the
// cell. There are two wedges associated with an edge of cell face in
// 3D.
//
// Corners are also subcell entities that are associated uniquely with
// a node of a cell. Each corner is the union of all the wedges incident
// upon that node in the cell
//
// Facets are the boundary entity between two wedges in adjacent
// cells. In 3D, a facet is a triangular subface of the cell face
// shared by two wedges in adjacent cells. In 2D, a facet is half of
// an edge that is shared by two wedges in adjacent cells
//


// TODO:  Right now we're relying on the fact that this enum is
//        identical to Jali::Entity_kind.  Need to fix this.
enum Entity_kind {
  ALL_KIND = -3,
  ANY_KIND = -2,
  UNKNOWN_KIND = -1,
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

enum Entity_type {
  TYPE_UNKNOWN = -1,  // Initializer
  DELETED = 0,
  PARALLEL_OWNED = 1,          // Owned by this processor
  PARALLEL_GHOST = 2,          // Owned by another processor
  BOUNDARY_GHOST = 3,          // Ghost/Virtual entity on boundary
  ALL  = 4     // PARALLEL_OWNED + PARALLEL_GHOST + BOUNDARY_GHOST
};

// Element (cell topology) type

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
