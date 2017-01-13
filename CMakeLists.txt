cmake_minimum_required(VERSION 2.8)
project(simplicial CXX)


# user libraries should include GUDHI
# requirements:
#   * Gudhi
#   * Boost
# include_directories($ENV{USER_LIBS})

include_directories("./lib")

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

set(STUPID_SOURCES "./src/main.cpp")
add_executable(my_stupid_test ${STUPID_SOURCES})

set(OTHER_STUPID_SOURCES "./src/libtest.cpp")
add_executable(my_other_stupid_test ${OTHER_STUPID_SOURCES})
