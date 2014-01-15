# supported minimum version of cmake
cmake_minimum_required(VERSION 2.6)

# name of the project
project(esfp)

# version number
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION.PATCH}")

# specify the default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_SOURCE_DIR}/src)


# compile static libraries
set(BUILD_SHARED_LIBS OFF)

################################################################################
# build dependencies
################################################################################
set(DEP_DIR ${CMAKE_SOURCE_DIR}/ext)

# Lua for configuration

# disable install_* commands, because only a static library is needed.
macro(install)
endmacro()

add_subdirectory(${DEP_DIR}/lua)
set_target_properties(lua luac PROPERTIES EXCLUDE_FROM_ALL ON)
macro(install)
_install(${ARGV})
endmacro(install)

################################################################################
# build main program
################################################################################
add_subdirectory(src)



################################################################################
# Test
################################################################################
enable_testing()

# google test (gtest) for unit tests
set(GTEST ${DEP_DIR}/gtest)
add_subdirectory(${GTEST})

add_subdirectory(test)
