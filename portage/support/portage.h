/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_SUPPORT_PORTAGE_H_
#define PORTAGE_SUPPORT_PORTAGE_H_

// Autogenerated file that contains configuration specific defines
// like PORTAGE_ENABLE_MPI and PORTAGE_ENABLE_THRUST
#include "portage-config.h"


#ifdef PORTAGE_ENABLE_THRUST

#include "thrust/device_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

#else  // no thrust

#include <boost/iterator/counting_iterator.hpp>

#include <vector>
#include <algorithm>
#include <string>

#endif

#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"
#include "wonton/support/Matrix.h"
#include "wonton/support/wonton.h"

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

using Wonton::Matrix;

// useful enums
using Wonton::Entity_kind;
using Wonton::Entity_type;
using Wonton::Element_type;
using Wonton::Field_type;
using Wonton::Data_layout;
using Wonton::Weights_t;

/// Limiter type
typedef enum {NOLIMITER, BARTH_JESPERSEN} Limiter_type;
constexpr int NUM_LIMITER_TYPE = 2;

Limiter_type DEFAULT_LIMITER = Limiter_type::BARTH_JESPERSEN;

std::string to_string(Limiter_type limiter_type) {
  static const std::string type2string[NUM_LIMITER_TYPE] =
      {"Limiter_type::NOLIMITER",
       "Limiter_type::BARTH_JESPERSEN"};

  int itype = static_cast<int>(limiter_type);
  return (itype >= 0 && itype < NUM_LIMITER_TYPE) ?
      type2string[itype] : "INVALID LIMITER TYPE";
}


/// Fixup options for partially filled cells
typedef enum {CONSTANT, LOCALLY_CONSERVATIVE, SHIFTED_CONSERVATIVE}
  Partial_fixup_type;
constexpr int NUM_PARTIAL_FIXUP_TYPE = 3;

Partial_fixup_type DEFAULT_PARTIAL_FIXUP_TYPE =
    Partial_fixup_type::LOCALLY_CONSERVATIVE;

std::string to_string(Partial_fixup_type partial_fixup_type) {
  static const std::string type2string[NUM_PARTIAL_FIXUP_TYPE] =
      {"Partial_fixup_type::CONSTANT",
       "Partial_fixup_type::LOCALLY_CONSERVATIVE",
       "Partial_fixup_type::SHIFTED_CONSERVATIVE"};

  int itype = static_cast<int>(partial_fixup_type);
  return (itype >= 0 && itype < NUM_PARTIAL_FIXUP_TYPE) ? type2string[itype] :
      "INVALID PARTIAL FIXUP TYPE";
}


/// Fixup options for empty cells
typedef enum {LEAVE_EMPTY, EXTRAPOLATE, FILL}
  Empty_fixup_type;
constexpr int NUM_EMPTY_FIXUP_TYPE = 3;

Empty_fixup_type DEFAULT_EMPTY_FIXUP_TYPE = Empty_fixup_type::LEAVE_EMPTY;

std::string to_string(Empty_fixup_type empty_fixup_type) {
  static const std::string type2string[NUM_EMPTY_FIXUP_TYPE] =
      {"Empty_fixup_type::LEAVE_EMPTY",
       "Empty_fixup_type::EXTRAPOLATE",
       "Empty_fixup_type::FILL"};

  int itype = static_cast<int>(empty_fixup_type);
  return (itype >= 0 && itype < NUM_EMPTY_FIXUP_TYPE) ? type2string[itype] :
      "INVALID EMPTY FIXUP TYPE";
}


/// Intersection and other tolerances to handle tiny values
struct NumericTolerances_t {
    bool   tolerances_set                   = false;
    double polygon_convexity_eps            = error_value_;
    double minimal_intersection_volume      = error_value_;
    double intersect_bb_eps                 = error_value_;
    double min_relative_volume              = error_value_;
    double mmdriver_relative_min_mat_vol    = error_value_;
    double mismatch_fixup_constant_value    = error_value_;
    double mismatch_conservation_all_mat    = error_value_;

    void use_default()
    {
        tolerances_set                  =   true;
        polygon_convexity_eps           =  1e-14;
        minimal_intersection_volume     = -1e-14;
        intersect_bb_eps                =  1e-12;
        min_relative_volume             =  1e-12;
        mmdriver_relative_min_mat_vol   =  1e-10;
        mismatch_fixup_constant_value   =  100*std::numeric_limits<double>::epsilon();
        mismatch_conservation_all_mat   =  100*std::numeric_limits<double>::epsilon();
    }

    private:
        double error_value_ = 1e5;
};

// Iterators and transforms that depend on Thrust vs. std
#ifdef PORTAGE_ENABLE_THRUST

template<typename T>
    using vector = thrust::device_vector<T>;

template<typename T>
    using pointer = thrust::device_ptr<T>;

typedef thrust::counting_iterator<unsigned int> counting_iterator;
inline counting_iterator make_counting_iterator(unsigned int const i) {
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

typedef boost::counting_iterator<unsigned int> counting_iterator;
inline counting_iterator make_counting_iterator(unsigned int const i) {
  return boost::make_counting_iterator<unsigned int>(i);
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

}  // namespace Portage

#endif  // PORTAGE_SUPPORT_PORTAGE_H_
