# supported minimum version of cmake
cmake_minimum_required(VERSION 2.6.2)

# name of the project
project(esf)

# version number
set(ESF_VERSION_MAJOR 0)
set(ESF_VERSION_MINOR 0)
set(ESF_VERSION_PATCH 0)

set(ESF_VERSION "${ESF_VERSION_MAJOR}.${ESF_VERSION_MINOR}.${ESF_VERSION.PATCH}")

# specify the default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


################################################################################
# build dependencies
################################################################################
set(DEP_DIR ${CMAKE_SOURCE_DIR}/ext)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Lua51 REQUIRED)

find_package(Eigen3 3.1.0 REQUIRED)
if (3.0.6 VERSION_GREATER ${EIGEN3_VERSION})
   message(FATAL_ERROR "Eigen3 is too old (required at least version 3.0.6)")
endif()

include_directories(${CMAKE_SOURCE_DIR}/src ${Lua51_INCLUDE_DIRS})

# Eigen3 generates lots of errors during build.  In order to suppress those errors,
# use -isystem instead of -I.
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

link_directories(${CMAKE_BINARY_DIR}/src ${Lua51_LIBRARIES})


################################################################################
# build main program
################################################################################
add_subdirectory(src)



################################################################################
# Test
################################################################################
# enable_testing()

# # google test (gtest) for unit tests
# set(GTEST ${CMAKE_SOURCE_DIR}/ext/gtest)
# add_subdirectory(${GTEST})

# add_subdirectory(test)
