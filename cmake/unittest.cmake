#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#----------------------------------------------------------------------------
# Google Test
#----------------------------------------------------------------------------

find_package(GTest QUIET)  # This will catch externally installed GTest

if (NOT GTest_FOUND OR NOT TARGET GTest::gtest)  # build from submodule
  find_package(Threads QUIET)  # Find thread libraries for system

  option(INSTALL_GTEST OFF)
  add_subdirectory(googletest)

  if (NOT TARGET GTest::gtest)
    if (TARGET gtest)
      # Add aliases that will allow us to refer to the targets the same way
      # whether we use an extern googletest installation or the one we have
      # as a submodule
      add_library(GTest::gtest ALIAS gtest)
      add_library(GTest::gtest_main ALIAS gtest_main)
    else ()
      message(FATAL_ERROR "Added googletest subdirectory but cannot find target GTest::gtest or gtest")
    endif ()
  endif ()
endif ()

#[===========================================================================[
.. command:: portage_add_unittest

The ``portage_add_unittest`` function creates a custom unit test with
various runtime policies::

portage_add_unittest(<name> [<option>...])

General options are:

``SOURCES <sources>...``
The sources necessary to build the test executable
``INPUTS <inputs>...``
The input files used to run the test
``POLICY <policy>``
The runtime policy to use when executing the test (SERIAL, MPI)
``THREADS <threads>...``
The number of threads to run the test with
``LIBRARIES <libraries>...``
List of libraries to link target against
``DEFINES <defines>...``
Defines to set when building test target
``ARGUMENTS <test arguments>
Arguments supplied to the test command line
#]===========================================================================]

function(portage_add_unittest name)

  #--------------------------------------------------------------------------#
  # Setup argument options.
  #--------------------------------------------------------------------------#

  set(options)
  set(one_value_args POLICY)
  set(multi_value_args SOURCES INPUTS LIBRARIES DEFINES THREADS)
  cmake_parse_arguments(test "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})
  
  #--------------------------------------------------------------------------#
  # Make sure that the user specified sources and the list contains Main.cc
  #--------------------------------------------------------------------------#
  
  if (NOT test_SOURCES)
    message(FATAL_ERROR "You must specify unit test source files using SOURCES")
  endif ()

  # Include the correct main file
  if (test_POLICY STREQUAL "MPI")
    list(APPEND test_SOURCES ${CMAKE_SOURCE_DIR}/cmake/test_main_mpi.cc)
  endif ()

  
  # Set up the unit test executable and dependencies
  
  add_executable(${name} ${test_SOURCES})
  target_link_libraries(${name} PRIVATE GTest::gtest)

  if (NOT test_POLICY OR test_POLICY STREQUAL "SERIAL")
    # use the standard gtest main program
    target_link_libraries(${name} PRIVATE GTest::gtest_main)
  endif ()

  # Add any user supplied compiler definitions
  
  target_compile_definitions(${name} PRIVATE ${test_DEFINES})

  # Add in any extra libraries
  
  if (test_LIBRARIES)
    target_link_libraries(${name} PRIVATE ${test_LIBRARIES})
  endif ()

  # If test needs any input files, copy them to the testing directory
  if (INPUTS)
    set(_INPUT_FILES)
    foreach (input ${test_INPUTS})
      get_filename_component(_OUTPUT_NAME ${input} NAME)
      get_filename_component(_PATH ${input} ABSOLUTE)
      configure_file(${_PATH} ${PROJECT_BINARY_DIR}/${_OUTPUT_NAME})
      list(APPEND _INPUT_FILES ${PROJECT_BINARY_DIR}/${_OUTPUT_NAME})
    endforeach ()

    add_custom_target(${name}_inputs DEPENDS ${_INPUT_FILES})
    add_dependencies(${name} ${name}_inputs)
  endif ()

  
  if (test_POLICY STREQUAL "MPI")
    if (NOT test_THREADS)
      set(test_THREADS 1)
    endif ()
    
    foreach (instance ${test_THREADS})
      if (ENABLE_JENKINS_OUTPUT)
        set(_OUTPUT ${name}_${instance}.xml)
        set(_GTEST_FLAGS "--gtest_output=xml:${_OUTPUT}")
      endif()
      
      add_test(NAME ${name}
        COMMAND
	${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${test_THREADS}
        $<TARGET_FILE:${name}> ${test_ARGUMENTS}
        ${GTEST_FLAGS}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    endforeach ()

  else()

    if (ENABLE_JENKINS_OUTPUT)
      set(_OUTPUT ${PROJECT_BINARY_DIR}/${name}.xml)
      set(GTEST_FLAGS "--gtest_output=xml:${_OUTPUT}")
    endif()

    add_test(
      NAME ${name}
      COMMAND
      $<TARGET_FILE:${name}>
      ${GTEST_FLAGS}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  endif ()

endfunction (portage_add_unittest)
