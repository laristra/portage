/*!
  @file DummyInterfaceReconstructor.h
  @brief A placeholder for a real interface reconstructor from Tangram
*/
#ifdef HAVE_TANGRAM

#include "tangram/support/tangram.h"
#include "tangram/driver/CellMatPoly.h"


namespace Portage {

// Dummy interface reconstruction class, if external interface
// reconstruction methods are not found
template<class Mesh_Wrapper, int Dim> class DummyInterfaceReconstructor {
 public:
  DummyInterfaceReconstructor(Mesh_Wrapper const & mesh) {}
  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs,
                            std::vector<Tangram::Point<Dim>> const& cell_mat_centroids) {}
  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs) {}

  void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {}
  std::shared_ptr<Tangram::CellMatPoly<Dim>> operator()(const int cell_op_ID) const {}
};

}

#endif  // HAVE_TANGRAM
