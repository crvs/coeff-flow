cmake_minimum_required(VERSION 2.8)
project(simplicial CXX)


# requirements: CHANGE LOCATIONS HERE!!!
#   * Gudhi
#   * Eigen3
set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
#   * boost
# end requirementes

# header files

set(CMAKE_BUILD_TYPE Debug)
include_directories("./lib")
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

set(STUPID_SOURCES "./src/main.cpp")
add_executable(my_stupid_test ${STUPID_SOURCES})

set(YOTHER_STUPID_SOURCES "./src/graphtest.cpp")
add_executable(yother_stupid_test ${YOTHER_STUPID_SOURCES})

# set(CMAKE_CXX_FLAGS "-lspatialindex_c -lspatialindex ${CMAKE_CXX_FLAGS}")
set(OTHER_STUPID_SOURCES "./src/libtest.cpp")
add_executable(my_other_stupid_test ${OTHER_STUPID_SOURCES})
