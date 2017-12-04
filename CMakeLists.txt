cmake_minimum_required(VERSION 2.8)
project(blstools)

# project version
set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 3)
set(${PROJECT_NAME}_PATCH_LEVEL 0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a more recent C++ compiler.")
endif()

find_package(BLAS REQUIRED)
find_package(Threads REQUIRED)
find_package(CUDA)

# write the config.h configuration file
configure_file(cmake/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)
add_definitions(-DHAVE_CONFIG_H)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/src" ${CUDA_INCLUDE_DIRS})

add_executable(blstools src/dict.cpp src/hist.cpp src/blstools.cpp src/bls.cpp
        src/phylotree.cpp src/matrix.cpp src/pwmscan.cpp src/motif.cpp
        src/sequence.cpp  src/ortho.cpp src/species.cpp)
target_link_libraries(blstools ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS blstools RUNTIME DESTINATION bin)
