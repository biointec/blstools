cmake_minimum_required(VERSION 2.6)
project(blstools)

# project version
set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 1)
set(${PROJECT_NAME}_PATCH_LEVEL 0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# check the Fortran / C name mangling scheme
include(cmake/FortranScheme.cmake)

# set some definitions
add_definitions("-DBLSTOOLS_MAJOR_VERSION=${${PROJECT_NAME}_MAJOR_VERSION}")
add_definitions("-DBLSTOOLS_MINOR_VERSION=${${PROJECT_NAME}_MINOR_VERSION}")
add_definitions("-DBLSTOOLS_PATCH_LEVEL=${${PROJECT_NAME}_PATCH_LEVEL}")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a recent C++ compiler.")
endif()

find_package(BLAS REQUIRED)

# write the config.h configuration file
configure_file(cmake/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)
add_definitions(-DHAVE_CONFIG_H)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/src")

add_executable(blstools src/blstools.cpp src/bls.cpp src/phylotree.cpp
        src/matrix.cpp src/pwmscan.cpp src/motif.cpp src/sequence.cpp)
target_link_libraries(blstools ${BLAS_LIBRARIES})

install(TARGETS blstools RUNTIME DESTINATION bin)
