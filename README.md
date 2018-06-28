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

### Running with docker

A docker file is provided in the `Docker` folder, and can be built by simply running `docker build -t coeff-flow .` from the `Docker` folder. Alternative an image of the same file can be downloaded from docker hub via `docker pull crvsf/coeff-flow` and there you can run the tests by:

```{bash}
docker run -it crvsf/coeff-flow
cd /coeff-flow build
make timing-test
```

to run the timing tests, the files for the other tests, are still not provided.

### Using Python bindings

In the python bindings we provide bindings for the `simplicial_complex` which uses only cells to for construction, and provides also bindings for the functions `cell_to_index` and `index_to_cell`. Furthermore we provide bindings for `coeff_flow` and `coeff_flow_embedded`. Below is a toy example run:

```{python}
>>> import coeffflow
>>> a = coeffflow.simplicial_complex([[0,1,2]])
>>> a.cell_to_index([0,2])
1
>>> a.index_to_cell(0,1)
[1]
>>> a.index_to_cell(1,2)
[2, 1]
>>> coeffflow.coeff_flow_embedded(a ,(1, [1,-1,1]))
(2, [1.0])
```

### References/links

[1]: Carvalho JF, Vejdemo-Johansson M, Kragic D, Pokorny FT. _An algorithm for calculating top-dimensional bounding chains._ 2017 [PeerJ Preprints 5:e3151v1](https://doi.org/10.7287/peerj.preprints.3151v1)

[2]: The Stanford 3D Scanning Repository https://graphics.stanford.edu/data/3Dscanrep/#bunny

[3]: Smithsonian X3D https://3d.si.edu/explorer/eulaema-bee#downloads
