# supported minimum version of cmake
cmake_minimum_required(VERSION 2.6)

# name of the project
project(esf_multi)

# specify the default flags
set(CMAKE_CXX_FLAGS "-std=c++11" CACHE STRING "default compiler flags" FORCE)
include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_SOURCE_DIR}/src)

# version number
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION.PATCH}")

# location of dependency
set(DEP_DIR ${CMAKE_SOURCE_DIR}/ext)

add_subdirectory(src)

# build
enable_testing()

set(GTEST ${DEP_DIR}/gtest)
add_subdirectory(${GTEST})

add_subdirectory(test)
