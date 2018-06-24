[![DOI](https://zenodo.org/badge/78863905.svg)](https://zenodo.org/badge/latestdoi/78863905)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Build](https://travis-ci.org/crvs/coeff-flow.svg?branch=master)](https://travis-ci.org/crvs/coeff-flow)

# Coefficient Flow algorithm

## Introduction
This repository imelements the `Coefficient_Flow` algorithm presented in [1].

## Compiling

In order to compile this repository you will need the following packages:

- `cmake (v >= 2.8)`
- `clang (v == 3.8)`
- `Eigen (v == 3.3.2)`
- `Boost (v >= 1.54.0)`
- `Gudhi (v >= 1.3.1)`

**Note:**
- In case you are using a different version of `clang`, everything should still compile without problems, but you will need to change the appropriate lines in `CMakeLists.txt`.

## Running

The timing tests presented in [1] can be run easily by simply performing

```{bash}
mkdir build && cd build
cmake .. && make
make timing_test
```
this will take a long time to run as it will run a test for a random mesh comprising (about) `x 1ey` points, with `x in [1..9]` and `y in [1..5]`. The results of the test are output to the file `results.csv`.

For the other tests presented in [1] we only provide the `yaml` files needed to run them, as the meshes are provided by a third party and can be found in [here [2]](https://graphics.stanford.edu/data/3Dscanrep/#bunny) and [here [3]](https://3d.si.edu/explorer/eulaema-bee#downloads). The aforementioned `yaml` files are stored in the `/test` folder and can be run (after the project has been built, and starting from the build directory, and making sure that the proper mesh is located in the same folder)

```{bash}
cp ../test/bunny2.yaml .
./yamltest bunny2.yaml
```

This will output two `ply` files which contain a copy of the mesh `bunzipper.ply` and the bounding chain to the cycle specified by the path contained in the yaml file. One of them was computed by solving a large linear system using least squares conjugate gradient descent, whereas the other was computed by employing `coefficient_flow`.

### References/links

[1]: Carvalho JF, Vejdemo-Johansson M, Kragic D, Pokorny FT. _An algorithm for calculating top-dimensional bounding chains._ 2017 [PeerJ Preprints 5:e3151v1](https://doi.org/10.7287/peerj.preprints.3151v1)

[2]: The Stanford 3D Scanning Repository https://graphics.stanford.edu/data/3Dscanrep/#bunny

[3]: Smithsonian X3D https://3d.si.edu/explorer/eulaema-bee#downloads
