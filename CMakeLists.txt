cmake_minimum_required(VERSION 2.8)
project(simplicial CXX)


# requirements: CHANGE HERE!!!
#   * Gudhi
include_directories("./lib")
#   * Eigen3
set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
# end requirementes

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

set(STUPID_SOURCES "./src/main.cpp")
add_executable(my_stupid_test ${STUPID_SOURCES})

set(OTHER_STUPID_SOURCES "./src/libtest.cpp")
add_executable(my_other_stupid_test ${OTHER_STUPID_SOURCES})
