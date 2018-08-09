/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

/*!
  @file DummyInterfaceReconstructor.h
  @brief A placeholder for a real interface reconstructor from Tangram
*/
#ifndef DUMMY_INTERFACE_RECONSTRUCTOR_H
#define DUMMY_INTERFACE_RECONSTRUCTOR_H


#ifdef HAVE_TANGRAM
#include "tangram/support/tangram.h"
#include "tangram/driver/CellMatPoly.h"
#endif

namespace Portage {

// Dummy interface reconstruction class, if external interface
// reconstruction methods are not found
template<class Mesh_Wrapper, 
         int Dim, 
         class MatPoly_Splitter=void, 
         class MatPoly_Clipper=void>
class DummyInterfaceReconstructor {
 public:
  DummyInterfaceReconstructor(Mesh_Wrapper const& mesh) {}

#ifdef HAVE_TANGRAM
  DummyInterfaceReconstructor(Mesh_Wrapper const& mesh,
                              const Tangram::IterativeMethodTolerances_t& im_tols,
                              const bool all_convex = false) {}

  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs,
                            std::vector<Tangram::Point<Dim>> const& cell_mat_centroids) {}
#endif

  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs) {}

  void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {}

#ifdef HAVE_TANGRAM
  std::shared_ptr<Tangram::CellMatPoly<Dim>> operator()(const int cell_op_ID) const {}
#endif

};

}

#endif  // DUMMY_INTERFACE_RECONSTRUCTOR_H
