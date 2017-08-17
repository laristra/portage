



# Comparator for testing infrastructure
add_executable(apptest_cmp apptest_cmp.cc)

if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
