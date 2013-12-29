project(${CMAKE_PROJECT_NAME}_test)

include_directories("${GTEST}/include")
link_directories(${GTEST})

set(test_SRC
  allele_test.cc
  afs_test.cc
)

add_executable(${PROJECT_NAME} ${test_SRC})

target_link_libraries(${PROJECT_NAME} gtest gtest_main ${LIB_NAME})

add_test(UnitTest ${PROJECT_NAME})
