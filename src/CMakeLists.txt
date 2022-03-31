# Copyright (c) 2022 Chan Beom Park <cbpark@gmail.com>

set(YAM2_SRC
  constraint.cc
  gradient.cc
  input.cc
  invisibles.cc
  momentum.cc
  yam2.cc)
add_library(YAM2_OBJ OBJECT ${YAM2_SRC})
set_property(TARGET YAM2_OBJ PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(YAM2_OBJ
  PUBLIC
  "${PROJECT_SOURCE_DIR}/include"
  PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(YAM2_shared SHARED $<TARGET_OBJECTS:YAM2_OBJ>)
set_target_properties(YAM2_shared PROPERTIES OUTPUT_NAME YAM2)
add_library(YAM2_static STATIC $<TARGET_OBJECTS:YAM2_OBJ>)
set_target_properties(YAM2_static PROPERTIES OUTPUT_NAME YAM2)
target_link_libraries(YAM2_shared YAM2_static)

install(TARGETS YAM2_shared YAM2_static DESTINATION lib)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/YAM2"
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.h")
