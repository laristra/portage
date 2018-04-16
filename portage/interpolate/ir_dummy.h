#ifndef SRC_IR_DUMMY_H_
#define SRC_IR_DUMMY_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/Matrix.h"
#include "portage/support/lsfits.h"

#include "tangram/driver/driver.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"

namespace Portage {
  namespace detail{ 
    // TODO: move this to a common header that can be shared between 
    // mmdriver.h (which also uses DummyInterfaceReconstructor)

    // Dummy interface reconstruction class, if external interface
    // reconstruction methods are not found
    template<class Mesh_Wrapper, int Dim> class DummyInterfaceReconstructor {
    public:
      DummyInterfaceReconstructor(Mesh_Wrapper const & mesh) {};
      void set_volume_fractions(std::vector<int> const& cell_num_mats,
				std::vector<int> const& cell_mat_ids,
				std::vector<double> const& cell_mat_volfracs,
				std::vector<Tangram::Point<Dim>> const& cell_mat_centroids) {}; 
      void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
				std::vector<double> const& cell_mat_volfracs) {};
      
      void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {}
      std::shared_ptr<Tangram::CellMatPoly<Dim>> operator()(const int cell_op_ID) const {}
    };
  }
}
#endif // SRC_IR_DUMMY_H_
