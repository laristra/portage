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
    intersect_polys_r2d.h
    intersect_r2d.h
    r3d.h
    intersect_polys_r3d.h
    intersect_r3d.h
    dummy_interface_reconstructor.h
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

  endif (Jali_DIR)

  cinch_add_unit(intersect_r2d
    SOURCES intersect_r2d_test.cc r2d.c
    LIBRARIES 
    POLICY SERIAL
    )
    
  cinch_add_unit(intersect_r3d
    SOURCES intersect_r3d_test.cc r3d.c
    LIBRARIES portage  
    POLICY SERIAL
    )

  #[[  
  if (TANGRAM_FOUND AND XMOF2D_FOUND)
       cinch_add_unit(test_tangram_intersect
       SOURCES test/test_tangram_intersect
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} ${XMOF2D_LIBRARIES}
       POLICY MPI
       THREADS 1)       
  endif (TANGRAM_FOUND AND XMOF2D_FOUND)
  ]]

  if (TANGRAM_FOUND AND XMOF2D_FOUND)
      include_directories(${TANGRAM_INCLUDE_DIRS})
      cinch_add_unit(test_tangram_intersect_2D
       SOURCES test/test_tangram_intersect_2D
       LIBRARIES portage 
       POLICY SERIAL
       )       
      cinch_add_unit(test_tangram_intersect_3D
       SOURCES test/test_tangram_intersect_3D
       LIBRARIES portage 
       POLICY SERIAL
       )       
  endif (TANGRAM_FOUND AND XMOF2D_FOUND)
     
 

endif(ENABLE_UNIT_TESTS)

