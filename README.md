/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

[![Build Status](https://travis-ci.org/laristra/portage.svg?branch=master)](https://travis-ci.org/laristra/portage)
[![codecov.io](https://codecov.io/github/laristra/portage/coverage.svg?branch=master)](https://codecov.io/github/laristra/portage/portage?branch=master)
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=portage%3A%2Fmaster)](https://sonarqube.com/dashboard?id=portage%3A%2Fmaster)

# portage

The portage library provides a framework for general purpose data
remapping -- between meshes, between particles, or between meshes and
particles --- in computational physics applications.  Remapping is
facilitated thorugh the use of user-supplied _wrappers_ around
meshes/particle swarms with their data, and is broken into three
phases that operate through the wrappers: _search_ for intersection
candidates, perform the _intersection_ with candidates, then
_interpolate_ the results onto the new mesh or particle swarm.
Algorithms for each of the phases can be customized (e.g. order of
accuracy of the interpolation) and, throught he wrappers, take
advantage of hybrid parallelism (MPI+X).

## Getting Started

To obtain a copy of portage and its submodules, clone recursively:

```sh
git clone --recursive https://github.com/laristra/portage
```

### Prerequisites

portage uses some standard C++14 features, so a fairly modern compiler
is needed.  We regularly test with intel 17+ or gcc 5.3+, and openmpi
1.10.3+, however, the code and most tests can be built without MPI
support.  The build system _requires_ CMake version 3.0+.

The following libraries are also _required_:

- LAPACKE
- Boost (currently just for counting iterators)

portage provides wrappers for a few existing mesh types.  Building
support for these is _optional_:

- [Jali](http://github.com/lanl/jali): regularly tested with verison
  0.9.8.  You will need to set the `Jali_Dir` CMake variable if you
  wish to build support for Jali and its tests (see examples below).
- [FleCSI Burton Specialization](http://github.com/laristra/flecsi-sp):
  the Burton specialization in the `flecsi-sp` repository is built on
  top of [FleCSI](http://github.com/laristra/flecsi).  You will need
  both projects to build support for the Burton mesh specialization
  and its tests.  You will need to set the `FLECSI_INSTALL_DIR` CMake
  variable, see examples below.


### Installation instructions

Below we list copy & paste instructions for several local machines; we
have a script that parses this README file to execute the examples
below to ensure they build. You can easily adapt them for other
machines.

## Darwin

Execute the following from the portage root directory:

```c++
# darwin-fe
module load openmpi/2.0.1-intel_17.0.0 cmake/3.7.1
JALI_INSTALL_PREFIX=/projects/ngc/private/jali/0.9.8-intel-17.0.0-openmpi-2.0.1
TPL_INSTALL_PREFIX=/projects/ngc/private/jali-tpl/1.0.9-intel-17.0.0-openmpi-2.0.1
mkdir build
cd build
cmake \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D CMAKE_CXX_COMPILER=`which mpiCC` \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D Boost_INCLUDE_DIR:PATH=$TPL_INSTALL_PREFIX/include \
    ..
make -j16
ctest -j16 --output-on-failure
```

## Moonlight

Execute the following from the portage root directory:

```c++
# ml-fey
module load intel/17.0.1 openmpi/1.10.5 cmake
JALI_INSTALL_PREFIX=/usr/projects/ngc/private/jali/0.9.8-intel-17.0.1-openmpi-1.10.5
mkdir build
cd build
cmake \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D CMAKE_CXX_COMPILER=`which mpiCC` \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    ..
make -j16
ctest -j16 --output-on-failure
```

## Varan

Execute the following from the portage root directory:

```c++
# varan
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/17.0.1 openmpi/1.10.5 cmake
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.8-intel-17.0.1-openmpi-1.10.5
LAPACKE_INCLUDE_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-intel-17.0.1/include
LAPACKE_LIBRARY_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-intel-17.0.1
mkdir build
cd build
cmake \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D CMAKE_CXX_COMPILER=`which mpiCC` \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D PC_LAPACKE_NCLUDE_DIRS=$LAPACKE_INCLUDE_DIR \
    -D PC_LAPACKE_LIBRARY_DIRS=$LAPACKE_LIBRARY_DIR \
    ..
make -j2
ctest -j2 --output-on-failure
```

---

If you want to build an app that uses
[FleCSI](https://github.com/losalamos/flecsi), you can link against a built
verison of FleCSI on Varan.  An example is below:

```c++
# varan::flecsi
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load gcc/5.3.0 openmpi/1.10.3 cmake
FLECSI_INSTALL_PREFIX=/usr/local/codes/ngc/private/flecsi/gcc5.3_openmpi1.10.3
FLECSISP_INSTALL_PREFIX=/usr/local/codes/ngc/private/flecsi-sp/gcc5.3_openmpi1.10.3
LAPACKE_INCLUDE_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-gcc-5.3.0/include
LAPACKE_LIBRARY_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-gcc-5.3.0
mkdir build-flecsi
cd build-flecsi
cmake \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D CMAKE_CXX_COMPILER=`which mpiCC` \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D ENABLE_FleCSI=True \
    -D CMAKE_PREFIX_PATH="$FLECSI_INSTALL_PREFIX;$FLECSISP_INSTALL_PREFIX" \
    -D PC_LAPACKE_NCLUDE_DIRS=$LAPACKE_INCLUDE_DIR \
    -D PC_LAPACKE_LIBRARY_DIRS=$LAPACKE_LIBRARY_DIR \
    ..
make -j2
ctest -j2 --output-on-failure
```
