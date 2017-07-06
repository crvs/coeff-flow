# Coefficient Flow algorithm

## Introduction
This repository imelements the `Coefficient_Flow` algorithm presented in [1].

## Compiling

In order to compile this repository you will need the following packages:

- `cmake (v >= 2.8)`
- `clang (v == 3.8)`
- `Eigen (v == 3.3.2)`
- `Boost (v >= 1.54.0)`
- `CGAL  (v >= 4.9.1)`
- `Gudhi (v >= 1.3.1)`
- `pcl   (v >= 1.8.0)`
- `libyaml`

**Note:**
- In case you are using a different version of `clang`, everything should still compile without problems, but you will need to change the appropriate lines in `CMakeLists.txt`.
- pcl needs to be version 1.8 or higher because 1.7 does not support reading ply files into `pcl::PolygonMesh` elements.

## Running

The timing tests presented in [1] can be run easily by simply performing

```{bash}
mkdir build && cd build
cmake .. && make
make timing_test
```
this will take a long time to run as it will run a test for a random mesh comprising (about) `x 1ey` points, with `x in [1..9]` and `y in [1..5]`. The results of the test are output to the file `results.csv`.


### References

[1]: J.F. Carvalho, M. Vejdemo-Johansson, D. Kragic, F.T. Pokorny; An algorithm for calculating top-dimensional bounding chains. 2017, Unpublished.
