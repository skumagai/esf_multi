project(${CMAKE_PROJECT_NAME}_test)

include_directories(SYSTEM "${GTEST}/include")
link_directories(${GTEST})

set(test_SRC
  allele_test.cc
  afs_test.cc
  esf_prob_test.cc
  hit_prob_test.cc
  init_test.cc
  param_test.cc
  state_test.cc
  util_test.cc
)

add_executable(${PROJECT_NAME} ${test_SRC})

target_link_libraries(${PROJECT_NAME} gtest gtest_main ${LIB_NAME})

add_test(AlleleTest ${PROJECT_NAME} --gtest_filter="AlleleTest.*")

add_test(AFSTets ${PROJECT_NAME} --gtest_filter="AFSTest.*")

add_test(ESFProbTest ${PROJECT_NAME} --gtest_filter="ESFProbTest.*")

add_test(HitProbTest ${PROJECT_NAME} --gtest_filter="HitProbTest.*")

add_test(InitTest ${PROJECT_NAME} --gtest_filter="InitTest.*")

add_test(ParamTest ${PROJECT_NAME} --gtest_filter="ParamTest.*")

add_test(StateTest ${PROJECT_NAME} --gtest_filter="StateTest.*")

add_test(UtilTest ${PROJECT_NAME} --gtest_filter="UtilTest.*")
