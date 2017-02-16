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
include_directories(${EIGEN3_INCLUDE_DIR})
set(CGAL_INCLUDE_DIR "/usr/local/include")
include_directories(${CGAL_INCLUDE_DIR})

find_package(CGAL REQUIRED COMPONTENTS Core)
include(${CGAL_USE_FILE})
#
# end_requirements


# Using CLANG instead of gcc for reasons
set(CMAKE_C_COMPILER             "/usr/bin/clang")
set(CMAKE_C_FLAGS                "-Wall -std=c99")
set(CMAKE_C_FLAGS_DEBUG          "-g")
set(CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_CXX_COMPILER             "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
# Linker flags for
# set(CMAKE_CXX_FLAGS "-lspatialindex_c -lspatialindex ${CMAKE_CXX_FLAGS}")

set(CMAKE_AR      "/usr/bin/llvm-ar")
set(CMAKE_LINKER  "/usr/bin/llvm-ld")
set(CMAKE_NM      "/usr/bin/llvm-nm")
set(CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
set(CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

set(CMAKE_BUILD_TYPE Debug)

include_directories("./lib") # where the header files are

add_library(scomplex SHARED "./lib/scomplex/simplicial_complex.cpp")
add_library(coeffflow SHARED "./lib/scomplex/coeff_flow.cpp")

add_executable(impl_test "src/impl_test.cpp")
target_link_libraries(impl_test scomplex)

# set(STUPID_SOURCES "./src/main.cpp")
# add_executable(my_stupid_test ${STUPID_SOURCES})

# set(YOTHER_STUPID_SOURCES "./src/graphtest.cpp")
# add_executable(yother_stupid_test ${YOTHER_STUPID_SOURCES})

# set(OTHER_STUPID_SOURCES "./src/libtest.cpp")
# add_executable(my_other_stupid_test ${OTHER_STUPID_SOURCES})

