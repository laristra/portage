/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef CORNER_GET_CENTROID_HH_
#define CORNER_GET_CENTROID_HH_

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/support/Point.h"

template<int D>
void corner_get_centroid(
    int cn, const Wonton::Jali_Mesh_Wrapper& mesh,
    Wonton::Point<D>& xcn);

#endif
