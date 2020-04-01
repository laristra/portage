#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(support_HEADERS
    portage.h
    weight.h
    basis.h
    operator.h
    operator_references.h
    faceted_setup.h
    PARENT_SCOPE
)


# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_weight
    SOURCES test/test_weight.cc
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_basis
    SOURCES test/test_basis.cc
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_operator
    SOURCES test/test_operator.cc
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_faceted_setup
    SOURCES test/test_faceted_setup.cc
    POLICY SERIAL
    )

endif(ENABLE_UNIT_TESTS)
