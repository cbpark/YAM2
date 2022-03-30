# Copyright (c) 2022 Chan Beom Park <cbpark@gmail.com>

cmake_minimum_required(VERSION 3.10)

project(YAM2
  VERSION 0.3.3
  DESCRIPTION "Yet another library for the M2 variables"
  HOMEPAGE_URL "https://github.com/cbpark/YAM2"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: Debug, Release." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra -pedantic")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
find_package(nlopt REQUIRED)
include_directories(${nlopt_INCLUDE_DIR})

find_package(ROOT QUIET COMPONENTS Physics)
if(${ROOT_FOUND})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${ROOT_INCLUDE_DIRS} -DHAS_ROOT")
endif()

add_subdirectory(src)
add_subdirectory(examples)
