/*
This file is part of the Ristra Tangram project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

// All the include files we want to check to make sure they won't
// cause multiply defined symbols. If some function does cause multiply
// defined symbol errors, the options are to
//
// 0. Make sure you have include guards so that multiple inclusions of
// the include file in a source file is ok
// 1. Move it to a .cc file
// 2. Make it inline
// 3. Make it static
// 4. Enclose it in its own namespace
//
// Option 3 and 4 will cause each translation unit (compiled source file) to
// have it's own copy of the function

#include "portage/distributed/mpi_bounding_boxes.h"
#include "portage/distributed/mpi_particle_distribute.h"

#include "portage/driver/coredriver.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/driver/mmdriver.h"
#include "portage/driver/parts.h"
#include "portage/driver/uberdriver.h"
#ifdef PORTAGE_HAS_TANGRAM
#include "portage/driver/write_to_gmv.h"
#endif

#include "portage/interpolate/gradient.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/interpolate/interpolate_nth_order.h"

#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/intersect/intersect_clipper.h"
#include "portage/intersect/intersect_polys_r2d.h"
#include "portage/intersect/intersect_polys_r3d.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/intersect/intersect_rNd.h"

#include "portage/search/BoundBox.h"
#include "portage/search/kdtree.h"
#include "portage/search/search_direct_product.h"
#include "portage/search/search_kdtree.h"
#include "portage/search/search_simple.h"

#include "portage/support/mpi_collate.h"
#include "portage/support/portage.h"
#include "portage/support/timer.h"

#ifdef WONTON_ENABLE_KOKKOS
  #include "portage/search/search_simple_points.h"
  #include "portage/search/search_points_by_cells.h"
  #include "portage/accumulate/accumulate.h"
  #include "portage/estimate/estimate.h"
  #include "portage/driver/driver_mesh_swarm_mesh.h"
  #include "portage/driver/driver_swarm.h"
  #include "portage/support/basis.h"
  #include "portage/support/faceted_setup.h"
  #include "portage/support/operator.h"
  #include "portage/support/weight.h"
#endif
