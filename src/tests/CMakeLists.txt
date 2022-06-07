# (c) Pawsey Supercomputing Research Centre, 2022
# Copyright by Pascal Jahan Elahi
# All rights reserved
#
# Contributed by Pascal Jahan Elahi
#
# This file is part of NBodylib.

set(tests
    test_kdtree
)

foreach(test ${tests})
  
  add_executable(${test} ${test}.cxx)
  target_compile_definitions(${test} PRIVATE ${NBODY_DEFINES})
  target_compile_options(${test} PRIVATE ${NBODY_CXX_FLAGS})
	target_include_directories(${test} PRIVATE ${NBODY_INCLUDE_DIRS})
  target_link_libraries(${test} ${GSL_LIBRARIES} nbodylib)
  if (NBODY_LINK_FLAGS)
    set_target_properties(${test} PROPERTIES LINK_FLAGS ${NBODY_LINK_FLAGS})
  endif()
endforeach()
