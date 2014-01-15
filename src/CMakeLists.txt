project(${CMAKE_PROJECT_NAME})


# enable warnings if they are available and reasonable.

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-Wall HAS_W_ALL)
if(HAS_W_ALL)
add_definitions(-Wall)
endif()

check_cxx_compiler_flag(-Wextra HAS_W_EXTRA)
if(HAS_W_EXTRA)
add_definitions(-Wextra)
endif()

check_cxx_compiler_flag(-Weverything HAS_W_EVERYTHING)
if(HAS_W_EVERYTHING)
add_definitions(-Weverything)
endif()

check_cxx_compiler_flag(-Wno-c++98-compat-pedantic HAS_W_NO_CXX98_COMPAT_PEDANTIC)
if(HAS_W_NO_CXX98_COMPAT_PEDANTIC AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
add_definitions(-Wno-c++98-compat-pedantic)
endif()

check_cxx_compiler_flag(-Wno-c++98-compat HAS_W_NO_CXX98_COMPAT)
if(HAS_W_NO_CXX98_COMPAT AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
add_definitions(-Wno-c++98-compat)
endif()

check_cxx_compiler_flag(-pedantic HAS_PEDANTIC)
if(HAS_PEDANTIC)
add_definitions(-pedantic)
endif()

set(LIB_NAME ${PROJECT_NAME})
set(LIB_NAME ${LIB_NAME} PARENT_SCOPE)

# suppress those warnings from Eigen
include_directories(SYSTEM ${DEP_DIR}/eigen)

set(LIB_SRC
  afs.cc
  allele.cc
  esf_prob.cc
  hit_prob.cc
  init.cc
  param.cc
  state.cc
)

add_library(${LIB_NAME} STATIC ${LIB_SRC})
