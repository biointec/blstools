cmake_minimum_required(VERSION 2.6)
project(blstools)

# project version
set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 1)
set(${PROJECT_NAME}_PATCH_LEVEL 0)

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

add_executable(blstools src/blstools.cpp src/bls.cpp src/phylotree.cpp)

install(TARGETS blstools RUNTIME DESTINATION bin)
