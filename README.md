[![Build Status](https://travis-ci.org/laristra/portage.svg?branch=master)](https://travis-ci.org/laristra/portage)
[![codecov.io](https://codecov.io/github/laristra/portage/coverage.svg?branch=master)](https://codecov.io/github/laristra/portage/portage?branch=master)
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=portage%3A%2Fmaster)](https://sonarqube.com/dashboard?id=portage%3A%2Fmaster)

# portage

The portage library provides a framework for general purpose data
remapping - between meshes, between particles, or between meshes and
particles - in computational physics applications.  Remapping is
facilitated through the use of user-supplied _wrappers_ around
meshes/particle swarms with their data.  The remap algorithm is
organized in three phases operating on the wrappers corresponding to
the original mesh/particles: _search_ for intersection candidates,
calculate the _intersection_ with candidates, then _interpolate_ the
results onto the new mesh or particle swarm.  Algorithms for each of
the phases can be customized (e.g. order of accuracy of the
interpolation) and, through the wrappers, take advantage of hybrid
parallelism (MPI+X).

## Getting Started

To obtain a copy of portage and its submodules from GitHub, clone
recursively:

```sh
git clone --recursive https://github.com/laristra/portage
```

If you are familiar with Docker, take a look at
our
[Dockerfile](https://github.com/laristra/portage/blob/master/docker/Dockerfile) for
a working build environment.  In particular, the Dockerfile builds off
of
the [portage-buildenv](https://github.com/laristra/portage-buildenv)
Dockerfile, and uses
our
[travis.yml](https://github.com/laristra/portage/blob/master/.travis.yml) file
with Travis CI.

### Prerequisites

Portage uses standard C++11 features, so a fairly modern compiler is
needed.  We regularly test with Intel 17+ or GCC 5.3+.  Utilizing the
full capabilities of portage will require an MPI implementation; we
regularly test with OpenMPI 1.10.3+ The build system _requires_ CMake
version 3.0+.

The following libraries are also _required_ (see examples below):

- LAPACKE (3.7.1+)

- **__Either__** Boost (1.53.0+) **__or__** Thrust (1.6.0+):
  We wrap some features of either one of these packages.  If you would
  like to run with OpenMP or TBB threads, then you _must_ use Thrust.

Portage provides wrappers for a few third-party mesh types.  Building
support for these is _optional_:

- [Jali](http://github.com/lanl/jali):

  We regularly test with verison 0.9.8.  You will need to set the
  `Jali_Dir` CMake variable if you wish to build support for Jali and
  its tests (see examples below).

- [FleCSI Burton Specialization](http://github.com/laristra/flecsi-sp):

  The Burton specialization in the `flecsi-sp` repository is built on
  top of [FleCSI](http://github.com/laristra/flecsi).  You will need
  _both_ projects to build support for the Burton mesh specialization
  and its tests.  You will need to set `ENABLE_FleCSI=True` and add
  the FleCSI and FleCSI-sp install paths to the `CMAKE_PREFIX_PATH`;
  see examples below.  Both FleCSI packages are under constant
  development.  This version of portage is known to work with hash
  `374b56b` of the FleCSI _stable_ branch, and hash `e78c594` of the
  FleCSI-SP _stable_ branch.

The [documentation](http://portage.lanl.gov) is built using doxygen (1.8+).

For more details regarding CMake settings, see
the [documentation](http://portage.lanl.gov) page.

### Installing

In the simplest case where you have the appropriate versions mentioned
above and Boost and LAPACKE are in the usual locations that CMake
searches, then the build step is:

```sh
portage $ mkdir build
portage $ cd build
portage/build $ cmake -DENABLE_APP_TESTS=True ..
portage/build $ make
```

This compiles the serial code and about a dozen application tests.  To
run the tests, simply execute

```sh
portage/build $ make test
```

If you wish to install the code into the `CMAKE_INSTALL_PREFIX` then
simply execute
```sh
portage/build $ make install
```

To build the documentation, one would configure with the
`-DENABLE_DOXYGEN=True` flag, and then `make doxygen`.

See the examples below, or the
[documentation](http://portage.lanl.gov) for more build instructions.

# License

This project is licensed under a modified 3-clause BSD license - see
the [LICENSE](https://github.com/laristra/portage/blob/master/LICENSE)
file for details.

# Release

This software has been approved for open source release and has been
assigned **LA-CC-16-084**.

----
----

# Example builds

Below we list copy & paste instructions for several local machines; we
have a script that parses this README file to execute the examples
below to ensure they build.

## Darwin

Execute the following from the portage root directory:

```c++
# machine=darwin-fe
module load openmpi/2.0.1-intel_17.0.0 cmake/3.7.1
JALI_INSTALL_PREFIX=/projects/ngc/private/jali/0.9.8-intel-17.0.0-openmpi-2.0.1
TPL_INSTALL_PREFIX=/projects/ngc/private/jali-tpl/1.0.9-intel-17.0.0-openmpi-2.0.1
mkdir build
cd build
cmake \
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

## Snow

Execute the following from the portage root directory:

```c++
# machine=sn-fey
module load intel/17.0.4 openmpi/2.1.2 cmake
JALI_INSTALL_PREFIX=/usr/projects/ngc/private/jali/0.9.8-intel-17.0.4-openmpi-2.1.2
TANGRAM_INSTALL_PREFIX=/usr/projects/ngc/private/tangram/133c1db580f-intel-17.0.4-openmpi-2.1.2
XMOF2D_INSTALL_PREFIX=/usr/projects/ngc/private/xmof2d/0.9-intel-17.0.4-openmpi-2.1.2
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D TANGRAM_DIR:FILEPATH=$TANGRAM_INSTALL_PREFIX \
    -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
    ..
make -j16
ctest -j16 --output-on-failure
```

## Varan

__Note the typo in the version of the build system we are using: it is
indeed `PC_LAPACKE_NCLUDE_DIRS`.  This will be fixed in the CMake
files in a coming release.__

Execute the following from the portage root directory:

```c++
# machine=varan
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/17.0.1 openmpi/1.10.5 cmake
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.8-intel-17.0.1-openmpi-1.10.5
TANGRAM_INSTALL_PREFIX=/usr/local/codes/ngc/private/tangram/475b813919f-intel-17.0.1-openmpi-1.10.5
XMOF2D_INSTALL_PREFIX=/usr/local/codes/ngc/private/xmof2d/0.9-intel-17.0.1-openmpi-1.10.5
LAPACKE_INCLUDE_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-intel-17.0.1/include
LAPACKE_LIBRARY_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-intel-17.0.1
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D TANGRAM_DIR:FILEPATH=$TANGRAM_INSTALL_PREFIX \
    -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
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
# machine=varan::flecsi
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load gcc/5.3.0 openmpi/1.10.3 cmake
FLECSI_INSTALL_PREFIX=/usr/local/codes/ngc/private/flecsi/gcc5.3_openmpi1.10.3
FLECSISP_INSTALL_PREFIX=/usr/local/codes/ngc/private/flecsi-sp/gcc5.3_openmpi1.10.3
TANGRAM_INSTALL_PREFIX=/usr/projects/ngc/private/tangram/133c1db580f-gcc-5.3.0-openmpi-1.10.3
XMOF2D_INSTALL_PREFIX=/usr/projects/ngc/private/xmof2d/0.9-gcc-5.3.0-openmpi-1.10.3
LAPACKE_INCLUDE_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-gcc-5.3.0/include
LAPACKE_LIBRARY_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-gcc-5.3.0
mkdir build-flecsi
cd build-flecsi
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D ENABLE_FleCSI=True \
    -D CMAKE_PREFIX_PATH="$FLECSI_INSTALL_PREFIX;$FLECSISP_INSTALL_PREFIX" \
    -D TANGRAM_DIR:FILEPATH=$TANGRAM_INSTALL_PREFIX \
    -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
    -D PC_LAPACKE_NCLUDE_DIRS=$LAPACKE_INCLUDE_DIR \
    -D PC_LAPACKE_LIBRARY_DIRS=$LAPACKE_LIBRARY_DIR \
    ..
make -j2
ctest -j2 --output-on-failure
```

