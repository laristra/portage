# -----------------------------------------------------------------------------
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#     https://github.com/laristra/portage/blob/master/LICENSE
#
# It has been adapted from cinch:
#     https://gitlab.lanl.gov/laristra/cinch
# -----------------------------------------------------------------------------

function(add_doxygen)

  #----------------------------------------------------------------------#
  # Find Doxygen
  #----------------------------------------------------------------------#

  find_package(Doxygen REQUIRED)

  #----------------------------------------------------------------------#
  # Create the output directory if it doesn't exist. This is where
  # the documentation target will be written.
  #
  # NOTE: This differs depending on whether or not the project is
  # a top-level project or not.  Subprojects are put under their
  # respective project names.
  #----------------------------------------------------------------------#

  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(CONFIG_INFOTAG)
  else()
    set(CONFIG_INFOTAG "${PROJECT_NAME}.")
  endif()

  if(CONFIG_INFOTAG)
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/doc/${PROJECT_NAME})
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/${PROJECT_NAME})
    endif(NOT EXISTS ${CMAKE_BINARY_DIR}/doc/${PROJECT_NAME})

    set(_directory ${CMAKE_BINARY_DIR}/doc/${PROJECT_NAME})

    #------------------------------------------------------------------#
    # This variable is used in the doxygen configuration file.  It
    # will be used in the configure_file call below.
    #------------------------------------------------------------------#

    set(${PROJECT_NAME}_DOXYGEN_TARGET ${PROJECT_NAME}/doxygen)

    #------------------------------------------------------------------#
    # Install under the main project name in its own directory
    #------------------------------------------------------------------#

    set(_install ${CMAKE_PROJECT_NAME}/${PROJECT_NAME})

    #------------------------------------------------------------------#
    # Add dependency for sub doxygen targets
    #------------------------------------------------------------------#

    if ( TARGET doxygen )
      add_dependencies(doxygen ${CONFIG_INFOTAG}doxygen)
    endif()

    if ( TARGET install-doxygen )
      add_dependencies(install-doxygen
          ${CONFIG_INFOTAG}install-doxygen)
    endif()

  else()

    #------------------------------------------------------------------#
    # Output target is in 'doc'
    #------------------------------------------------------------------#

    set(_directory ${CMAKE_BINARY_DIR}/doc)

    #------------------------------------------------------------------#
    # This variable is used in the doxygen configuration file.  It
    # will be used in the configure_file call below.
    #------------------------------------------------------------------#

    set(${PROJECT_NAME}_DOXYGEN_TARGET doxygen)

    #------------------------------------------------------------------#
    # Install in its own directory
    #------------------------------------------------------------------#

    set(_install ${CMAKE_PROJECT_NAME})
  endif(CONFIG_INFOTAG)

  #----------------------------------------------------------------------#
  # Create directory for intermediate files
  #----------------------------------------------------------------------#

  if(NOT EXISTS ${_directory}/.doxygen)
    file(MAKE_DIRECTORY ${_directory}/.doxygen)
  endif(NOT EXISTS ${_directory}/.doxygen)

  #----------------------------------------------------------------------#
  # Generate doxygen configuration file
  #----------------------------------------------------------------------#

  if(ENABLE_DOXYGEN_WARN)
    set(DOXYGEN_WARN YES)
  else()
    set(DOXYGEN_WARN NO)
  endif(ENABLE_DOXYGEN_WARN)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.conf.in
      ${_directory}/.doxygen/doxygen.conf)

  #----------------------------------------------------------------------#
  # Add the doxygen target
  #----------------------------------------------------------------------#

  add_custom_target(${CONFIG_INFOTAG}doxygen
      ${DOXYGEN} -u ${_directory}/.doxygen/doxygen.conf
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.conf.in)

  #----------------------------------------------------------------------#
  # Add install target
  #----------------------------------------------------------------------#

  add_custom_target(${CONFIG_INFOTAG}install-doxygen
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${_directory}/doxygen
      $ENV{DESTDIR}/${CMAKE_INSTALL_PREFIX}/share/${_install})

endfunction(add_doxygen)