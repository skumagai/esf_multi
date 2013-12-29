project(${CMAKE_PROJECT_NAME})

set(LIB_NAME ${PROJECT_NAME}_lib)
set(LIB_NAME ${LIB_NAME} PARENT_SCOPE)

message(STATUS ${LIB_NAME})

set(LIB_SRC
  # init.cc
  # state.cc
  allele.cc
  afs.cc
  # hit_prob.cc
)

add_library(${LIB_NAME} STATIC ${LIB_SRC})
