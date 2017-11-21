#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(intersect_HEADERS
    clipper.hpp
    intersectClipper.h
    r2d.h
    intersect_r2d.h
    r3d.h
    intersect_r3d.h
    PARENT_SCOPE
)

# Add source files
set(intersect_SOURCES
    clipper.cpp
    r2d.c
    r3d.c
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
  if (Jali_DIR)
    cinch_add_unit(intersectClipper
      SOURCES intersectClipperTest.cc clipper.cpp
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      )

    cinch_add_unit(intersect_r2d
      SOURCES intersect_r2d_test.cc r2d.c
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      )

    cinch_add_unit(intersect_r3d_jali
      SOURCES intersect_r3d_jali_test.cc r3d.c
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      # this is relative to thisCMakeLists.txt directory
      INPUTS ../../test_data/1x3d.exo ../../test_data/cube_0.exo ../../test_data/cube_0.5.exo ../../test_data/cube_1.exo ../../test_data/cube_2.exo
      POLICY MPI
      )
  endif (Jali_DIR)

  cinch_add_unit(intersect_r3d
    SOURCES intersect_r3d_test.cc r3d.c
    LIBRARIES portage 
    POLICY SERIAL
    )
endif(ENABLE_UNIT_TESTS)

