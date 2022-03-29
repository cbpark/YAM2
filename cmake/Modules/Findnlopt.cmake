# This sets the following variables:
#   NLOPT_FOUND
#   NLOPT_INCLUDE_DIRS
#   NLOPT_LIBRARIES

set(nlopt_DIR "" CACHE PATH "Root of the path to NLopt.")

find_path(nlopt_INCLUDE_DIR nlopt.h ${nlopt_DIR}/include)

set(nlopt_NAMES nlopt nlopt_cxx)
find_library(nlopt_LIBRARY
  NAMES ${nlopt_NAMES}
  PATHS ${nlopt_DIR}/lib ${nlopt_DIR}/lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(nlopt DEFAULT_MSG nlopt_LIBRARY nlopt_INCLUDE_DIR)

if(nlopt_FOUND)
  set(nlopt_LIBRARIES ${nlopt_LIBRARY})
else()
  set(nlopt_LIBRARIES)
endif()
