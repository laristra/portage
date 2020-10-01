/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_SUPPORT_PORTAGE_H_
#define PORTAGE_SUPPORT_PORTAGE_H_

// Autogenerated file that contains configuration specific defines
// like WONTON_ENABLE_MPI etc
#include "portage-config.h"

#include <vector>
#include <string>
#include <limits>
#include <algorithm>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"
#include "wonton/support/Matrix.h"

#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/support/tangram.h"
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

constexpr Limiter_type DEFAULT_LIMITER = Limiter_type::BARTH_JESPERSEN;

/// Boundary limiter type
typedef enum {BND_NOLIMITER, BND_ZERO_GRADIENT, BND_BARTH_JESPERSEN} Boundary_Limiter_type;
constexpr int NUM_Boundary_Limiter_type = 2;

constexpr Boundary_Limiter_type DEFAULT_BND_LIMITER = Boundary_Limiter_type::BND_NOLIMITER;

inline std::string to_string(Limiter_type limiter_type) {
  switch(limiter_type) {
    case NOLIMITER: return std::string("Limiter_type::NOLIMITER");
    case BARTH_JESPERSEN: return std::string("Limiter_type::BARTH_JESPERSEN");
    default: return std::string("INVALID LIMITER TYPE");
  }
}

inline std::string to_string(Boundary_Limiter_type boundary_limiter_type) {
  switch(boundary_limiter_type) {
    case BND_NOLIMITER: return std::string("Boundary_Limiter_type::BND_NOLIMITER");
    case BND_ZERO_GRADIENT: return std::string("Boundary_Limiter_type::BND_ZERO_GRADIENT");
    case BND_BARTH_JESPERSEN: return std::string("Boundary_Limiter_type::BND_BARTH_JESPERSEN");
    default: return std::string("INVALID BOUNDARY LIMITER TYPE");
  }
}

/// Fixup options for partially filled cells
typedef enum {CONSTANT, LOCALLY_CONSERVATIVE, GLOBALLY_CONSERVATIVE}
  Partial_fixup_type;
constexpr int NUM_PARTIAL_FIXUP_TYPE = 3;

constexpr Partial_fixup_type DEFAULT_PARTIAL_FIXUP_TYPE =
    Partial_fixup_type::LOCALLY_CONSERVATIVE;

inline std::string to_string(Partial_fixup_type partial_fixup_type) {
  static const std::string type2string[NUM_PARTIAL_FIXUP_TYPE] =
      {"Partial_fixup_type::CONSTANT",
       "Partial_fixup_type::LOCALLY_CONSERVATIVE",
       "Partial_fixup_type::GLOBALLY_CONSERVATIVE"};

  int itype = static_cast<int>(partial_fixup_type);
  return (itype >= 0 && itype < NUM_PARTIAL_FIXUP_TYPE) ? type2string[itype] :
      "INVALID PARTIAL FIXUP TYPE";
}


/// Fixup options for empty cells
typedef enum {LEAVE_EMPTY, EXTRAPOLATE, FILL}
  Empty_fixup_type;
constexpr int NUM_EMPTY_FIXUP_TYPE = 3;

constexpr Empty_fixup_type DEFAULT_EMPTY_FIXUP_TYPE = Empty_fixup_type::LEAVE_EMPTY;

inline std::string to_string(Empty_fixup_type empty_fixup_type) {
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
    // Flag if custom tolerances were used. If user is setting his own
    // tolerances, this flaq need to be set to true.
    bool   user_tolerances;

    // Check wheather the volume returned by r2d reduce is positive
    // (or slightly negative). If the volume is smaller, we throw an
    // error.
    double minimal_intersection_volume;

    // Distance tolerance: two points within that distance from each
    // other are considered coincident. Used for bounding box check
    // in Portage intersect and passed to Tangram in multi-material runs
    double min_absolute_distance;

    // Volume tolerance: intersections and material polytopes with
    // the volume below this tolerance are ignored. In multi-material
    // runs this tolerance is passed to Tangram. Target multi-material 
    // cells will not contain any material with volume below tolerance,
    // interface reconstruction results on the source mesh will not
    // contain material polytopes for materials with volume below this
    // tolerance.
    double min_absolute_volume;

    // Default relative tolerance on aggregated field values to detect
    // mesh mismatch
    double relative_conservation_eps;

    // Default number of iterations for mismatch repair
    int max_num_fixup_iter;
};

// Default values for tolerances
template <int D>
const NumericTolerances_t DEFAULT_NUMERIC_TOLERANCES = {
  false,                                           //user_tolerances
  -1.0e-14,                                        //minimal_intersection_volume
  sqrt(D)*std::numeric_limits<double>::epsilon(),  //min_absolute_distance: for two points
  // to be distinct, we need at least one coordinate to differ by machine epsilon or more.
  // In the worst case, when difference along all coordinate axes is the same, it corresponds to
  // the distance of sqrt(D)*machine_epsilon or more. Because we use distance criterion in 
  // order to not introduce a direction bias, we have to impose the worst case tolerance even 
  // though some points that are closer can technically be distinguished.
  std::numeric_limits<double>::epsilon(),          //min_absolute_volume
  100*std::numeric_limits<double>::epsilon(),      //relative_conservation_eps
  5                                                //max_num_fixup_iter
};


template<typename T>
bool contains(std::vector<T> const& haystack, T const& needle) {
  return std::find(haystack.begin(), haystack.end(), needle) != haystack.end();
}

}  // namespace Portage

#endif  // PORTAGE_SUPPORT_PORTAGE_H_
