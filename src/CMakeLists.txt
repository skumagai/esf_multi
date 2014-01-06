project(${CMAKE_PROJECT_NAME})

set(LIB_NAME ${PROJECT_NAME}_lib)
set(LIB_NAME ${LIB_NAME} PARENT_SCOPE)

include_directories(${DEP_DIR}/eigen)

set(LIB_SRC
  afs.cc
  allele.cc
  esf_prob.cc
  hit_prob.cc
  init.cc
  param.cc
  state.cc
  util.cc
)

add_library(${LIB_NAME} STATIC ${LIB_SRC})
