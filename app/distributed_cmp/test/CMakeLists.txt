#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
file(COPY distributed_cmp.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field.txt.0 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field.txt.1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field.txt.2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME distributed_cmp
         COMMAND ./distributed_cmp.sh)
set_property(TEST distributed_cmp
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_CURRENT_BINARY_DIR}/..)


