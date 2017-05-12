cmake_minimum_required(VERSION 2.8)
project(simplicial CXX)

# set(CMAKE_CXX_COMPILER_ID "Clang")
# message("${CMAKE_CXX_COMPILER_ID}")

# requirements: CHANGE LOCATIONS HERE!!!
#   * Gudhi under /usr/include/gudhi
#   * Eigen3 under /usr/include/eigen3
#   * boost under /usr/include/boost
#   * CGAL under /usr/include/CGAL
#
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")
set(EIGEN3_INCLUDE_DIR "/home/crvs/eigen-3.3.2")
include_directories(${EIGEN3_INCLUDE_DIR})
set(CGAL_INCLUDE_DIR "/usr/local/include")
include_directories(${CGAL_INCLUDE_DIR})


find_package(CGAL REQUIRED COMPONTENTS Core)
include(${CGAL_USE_FILE})
#
# end_requirements

#include_directories("/usr/lib/gcc/x86_65-linuxsgnu/4.8/include/omp")

# Using CLANG instead of gcc for reasons
set(CMAKE_C_COMPILER             "/usr/bin/clang-3.8")
set(CMAKE_C_FLAGS                "-Wall -std=c99")
set(CMAKE_C_FLAGS_DEBUG          "-g")
set(CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_CXX_COMPILER             "/usr/bin/clang++-3.8")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")# -fopenmp")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
# Linker flags for
# set(CMAKE_CXX_FLAGS "-lspatialindex_c -lspatialindex ${CMAKE_CXX_FLAGS}")

set(CMAKE_AR      "/usr/bin/llvm-ar")
set(CMAKE_LINKER  "/usr/bin/llvm-ld")
set(CMAKE_NM      "/usr/bin/llvm-nm")
set(CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
set(CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

# Build type
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

include_directories("./lib") # where the header files are

add_library(scomplex SHARED "./lib/scomplex/simplicial_complex.cpp")
add_library(pathsnap SHARED "./lib/scomplex/path_snapper.cpp")
add_library(clusterer SHARED "./lib/scomplex/clusterer.cpp")
target_link_libraries(clusterer pathsnap)

add_executable(test0 "src/test0.cpp")
target_link_libraries(test0 scomplex)

add_executable(test1 "src/test1.cpp")
target_link_libraries(test1 scomplex pathsnap)

add_executable(test2 "src/test2.cpp")
target_link_libraries(test2 scomplex pathsnap)

add_executable(test3 "src/test3.cpp")
target_link_libraries(test3 scomplex pathsnap)

add_executable(test4 "src/test4.cpp")
target_link_libraries(test4 scomplex pathsnap)

add_executable(test5 "src/test5.cpp")
target_link_libraries(test5 scomplex pathsnap)

add_executable(test6 "src/test6.cpp")
target_link_libraries(test6 scomplex pathsnap)

add_executable(test7 "src/test7.cpp")
target_link_libraries(test7 scomplex pathsnap clusterer)

find_library(HAVE_YAML NAME yaml-cpp)
add_executable(yamltest "src/testing_facility.cpp")
target_link_libraries(yamltest yaml-cpp scomplex pathsnap)
