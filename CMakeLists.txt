cmake_minimum_required(VERSION 2.8)
project(coefficient-flow CXX)

# requirements:
#   * Gudhi under /usr/include/gudhi
#   * Eigen3 under /usr/include/eigen3
#   * boost under /usr/include/boost
#   * CGAL under /usr/include/CGAL
#   * libyaml under /usr/local/include/libyaml
#   * pcl 1.8 under /usr/include/pcl-1.8
#

# in case Eigen, and/or CGAL, and/or GUDHI are not found by find_package
# try setting the paths to them manually
#
# set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")
# include_directories(${EIGEN3_INCLUDE_DIR})
# set(CGAL_INCLUDE_DIR "/usr/local/include/eigen3")
# include_directories(${CGAL_INCLUDE_DIR})
#

# GUDHI version should be 1.3 or higher, provided they don't break the
# SimplexTree interface
#
# libyaml-cpp this should actually be found in "/usr/include" so there is no
# real need for the find_package, but it's always nice to make sure
#
# pcl-1.8 version 1.7 doesn't let one read ply files into meshes which is
# rather strange, but it is a thing.
#

find_package(Boost REQUIRED QUIET)
message(STATUS "Boost version: ${BOOST_VERSION}")
message(STATUS "Boost directory: ${BOOST_INCLUDE_DIRS}")

find_package(Eigen3 3.1.2 REQUIRED QUIET)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
message(STATUS "Eigen3 directory: ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CGAL REQUIRED QUIET)
message(STATUS "CGAL version: ${CGAL_VERSION}")
message(STATUS "CGAL directory: ${CGAL_INCLUDE_DIRS}")

find_package(yaml-cpp REQUIRED QUIET)
message(STATUS "yaml-cpp version: ${YAML_CPP_VERSION}")
message(STATUS "yaml-cpp directory: ${YAML_CPP_INCLUDE_DIR}")
include_directories(${YAML_CPP_INCLUDE_DIR})

find_package(PCL 1.8 REQUIRED QUIET)
message(STATUS "PCL version: ${PCL_VERSION}")
message(STATUS "PCL directory: ${PCL_INCLUDE_DIRS}")
include_directories(${PCL_INCLUDE_DIRS})

find_package(GUDHI 1.3 REQUIRED QUIET)
message(STATUS "Gudhi version: ${GUDHI_VERSION}")
message(STATUS "Gudhi directory: ${GUDHI_INCLUDE_DIRS}")

set(CMAKE_CXX_COMPILER             "/usr/bin/clang++-3.8")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_AR      "/usr/bin/llvm-ar")
set(CMAKE_LINKER  "/usr/bin/llvm-ld")
set(CMAKE_NM      "/usr/bin/llvm-nm")
set(CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
set(CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

set(CMAKE_BUILD_TYPE Release)

include_directories("./lib") # where all the library files are

add_library(scomplex SHARED "./lib/scomplex/simplicial_complex.cpp")
add_library(pathsnap SHARED "./lib/scomplex/path_snapper.cpp")

add_executable(qhull2ply "src/make_mesh.cpp")

add_executable(yamltest "src/testing_facility.cpp")
target_link_libraries(yamltest yaml-cpp scomplex pathsnap ${PCL_LIBRARIES})

message(INFO ${CMAKE_CURRENT_SOURCE_DIR})
message(INFO ${CMAKE_CURRENT_BINARY_DIR})
# TEST target
add_custom_target(timing_test
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/test/run_dumb_tests.sh" "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/test/dumbexample.yaml" "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/run_dumb_tests.sh" "1"
    DEPENDS yamltest qhull2ply "test/run_dumb_tests.sh" "test/dumbexample.yaml")

add_library(quotient SHARED "./lib/scomplex/quotient.cpp")
add_executable(chaintest "src/chaintest.cpp")
add_executable(exttest "src/exttest.cpp")

