# Copyright (c) 2022 Chan Beom Park <cbpark@gmail.com>

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  foreach(exe m2 m2cons m2ccons m2vertex m2vertex_ineq m2consvertex m2consvertex_ineq)
    add_executable(${exe} ${exe}.cc)
    target_link_libraries(${exe} YAM2_shared ${nlopt_LIBRARY} ${ROOT_LIBRARIES})
  endforeach()
endif()

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION "share/YAM2/exmaples"
  FILES_MATCHING
  PATTERN "*.cc")
