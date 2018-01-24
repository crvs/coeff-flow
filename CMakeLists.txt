cmake_minimum_required(VERSION 2.8)
project(coefficient-flow CXX)

find_package(Boost REQUIRED QUIET)

find_package(Eigen3 3.1.2 REQUIRED QUIET)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CGAL REQUIRED QUIET)
include(${CGAL_USE_FILE})

find_package(yaml-cpp REQUIRED QUIET)
include_directories(${YAML_CPP_INCLUDE_DIR})

#find_package(PCL 1.8 REQUIRED QUIET)
#include_directories(${PCL_INCLUDE_DIRS})

find_package(GUDHI REQUIRED QUIET)
# definitions for GUDHI in case it is only locally compiled and not installed
if(NOT DEFINED(${GUDHI_INCLUDE_DIR}))
set(GUDHI_INCLUDE_DIR "${GUDHI_DIR}/include")
endif()
include_directories(${GUDHI_INCLUDE_DIR})

set(CMAKE_CXX_COMPILER             "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS                " -std=c++11 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_AR      "/usr/bin/llvm-ar")
set(CMAKE_LINKER  "/usr/bin/llvm-ld")
set(CMAKE_NM      "/usr/bin/llvm-nm")
set(CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
set(CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

set(CMAKE_BUILD_TYPE DEBUG)
# set(CMAKE_BUILD_TYPE RELEASE)

include_directories("./lib") # where all the library files are

# source available from
# https://github.com/ddiakopoulos/tinyply
add_library(tinyply SHARED "./lib/tinyply.cpp")

add_library(scomplex SHARED "./lib/SimpComplex.cpp")

add_library(pathsnap SHARED "./lib/path_snapper.cpp")
target_link_libraries(pathsnap scomplex CGAL)

add_executable(qhull2ply "src/make_mesh.cpp")

add_custom_target(timing_test
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/test/run_dumb_tests.sh" "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/test/dumbexample.yaml" "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/run_dumb_tests.sh" "1"
    DEPENDS yamltest qhull2ply "test/run_dumb_tests.sh" "test/dumbexample.yaml")

add_custom_target(tags
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/base_tags" "${CMAKE_CURRENT_SOURCE_DIR}/tags"
    COMMAND "ctags" "-R" "-a" "${CMAKE_CURRENT_SOURCE_DIR}/tags" "${CMAKE_CURRENT_SOURCE_DIR}")

add_custom_target(all_tags
    COMMAND "ctags" "-R" "-f" "${CMAKE_CURRENT_SOURCE_DIR}/base_tags" "/usr/include"
    COMMAND "ctags" "-R" "-a" "${CMAKE_CURRENT_SOURCE_DIR}/base_tags" "/usr/local/include"
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/base_tags" "${CMAKE_CURRENT_SOURCE_DIR}/tags"
    COMMAND "ctags" "-R" "-a" "${CMAKE_CURRENT_SOURCE_DIR}/tags" "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(quotient SHARED "./lib/quotient.cpp")
target_link_libraries(quotient scomplex)

add_library(boundchain SHARED "./lib/BoundingChain.cpp")
target_link_libraries(boundchain scomplex)

add_executable(chaintest "src/chaintest.cpp")
add_executable(exttest "src/exttest.cpp")

add_executable(yamltest "src/testing_facility.cpp")
target_link_libraries(yamltest yaml-cpp scomplex boundchain pathsnap tinyply)

include_directories(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} "/usr/include/python2.7/")
add_library(gsimp SHARED python/bindings.cpp)
target_link_libraries(gsimp "python2.7" "boost_python" ${PYTHON_LIBRARIES} scomplex pathsnap quotient ${CGAL_LIBRARIES} ${CGAL_CORE} ${Boost_LIBRARIES} )
