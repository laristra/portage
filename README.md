[![Build Status](https://travis-ci.com/laristra/portage.svg?branch=master)](https://travis-ci.com/laristra/portage)
[![codecov.io](https://codecov.io/github/laristra/portage/coverage.svg?branch=master)](https://codecov.io/github/laristra/portage/portage?branch=master)

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
needed.  We regularly test with Intel 18.0.1, GCC 6.4.0, and GCC 7.3.0.  
Utilizing the full capabilities of portage will require an MPI implementation; 
we regularly test with OpenMPI 2.1.2 The build system _requires_ CMake
version 3.13+.

The following libraries are also _required_ (see examples below):

- [Wonton](https://github.com/laristra/wonton)

  Wonton is a utility library for Portage containing some commonly
  used classes like Point, Vector, etc., some abstractions for on-node
  parallelism and also some mesh and state wrappers. Wonton itself has
  dependencies - the highly recommended one is LAPACKE
  (3.8.0+). On-node parallelism in Portage requires that Wonton be
  built with NVidia Thrust or Kokkos. Distributed parallelism requires
  that Wonton be built with MPI enabled. In the absence of the Thrust
  library, the Boost library must be linked into Wonton. See the
  Wonton README for details.
  
  If you specify, `PORTAGE_ENABLE_TANGRAM` and `TANGRAM_ROOT`,
  Wonton will be picked up automatically as Wonton is a dependency of
  Tangram as well. If not, you must specify the path to Wonton as
  `WONTON_ROOT`
  
  
The following libraries are _required if multi-material remapping is to be enabled_:

- [Tangram](https://github.com/laristra/tangram)

  Tangram is a material interface reconstruction library that is used
  to correctly remap between meshes with fractional amounts of
  multiple materials in some cells. In the future Tangram may become a
  required component of Portage.

The [documentation](https://portage.lanl.gov) is built using doxygen (1.8+).

For more details regarding CMake settings, see
the [documentation](https://portage.lanl.gov) page.

### Installing

In the simplest case the build step is:

```sh
portage $ mkdir build
portage $ cd build
portage/build $ cmake -DENABLE_APP_TESTS=True -DWONTON_ROOT:/path/to/wonton/installation ..
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
below to ensure they build. NOTE: If TANGRAM is **DISABLED**, the path
to WONTON must be specified explicitly using WONTON_ROOT or indicated
in CMAKE_PREFIX_PATH.

## Darwin

Execute the following from the portage root directory:

```sh
# machine=darwin-fe

INTEL_VERSION=18.0.3
MPI_VERSION=3.1.3
TANGRAM_VERSION=1.0.1
WONTON_VERSION=1.2.2

TPL_PREFIX=/usr/projects/ngc/private

# load the correct boost, compiler, and openmpi
module purge
module load cmake/3.15.3 openmpi/${MPI_VERSION}-intel_${INTEL_VERSION} boost/${BOOST_VERSION}

cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D PORTAGE_ENABLE_TANGRAM=True \
    -D WONTON_ROOT:FILEPATH=${TPL_PREFIX}/wonton/${WONTON_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION} \
    -D TANGRAM_ROOT:FILEPATH=${TPL_PREFIX}/tangram/${TANGRAM_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION} \
    -D PORTAGE_ENABLE_MPI=True \
    -D PORTAGE_ENABLE_THRUST=False \
	-D PORTAGE_ENABLE_Jali=True \
	-D PORTAGE_ENABLE_FleCSI=False \
    ..

make -j16
ctest -j16 --output-on-failure
```


## Snow

Execute the following from the portage root directory:

```sh
# machine=sn-fey
. /usr/share/lmod/lmod/init/sh

INTEL_VERSION=18.0.5
MPI_VERSION=2.1.2
WONTON_VERSION=dev
TANGRAM_VERSION=dev

TPL_PREFIX=/usr/projects/ngc/private

module load intel/${INTEL_VERSION} openmpi/${MPI_VERSION} cmake/3.14.6


mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D PORTAGE_ENABLE_TANGRAM=True \
    -D WONTON_ROOT:FILEPATH=${TPL_PREFIX}/wonton/${WONTON_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION} \
	-D TANGRAM_ROOT:FILEPATH=${TPL_PREFIX}/tangram/${TANGRAM_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION} \
    -D PORTAGE_ENABLE_MPI=True \
	-D PORTAGE_ENABLE_Jali=True \
	-D PORTAGE_ENABLE_FleCSI=False \
    ..
make -j4
ctest -j4 --output-on-failure
```

---

If you want to build an app for performance testing, you should include
Thrust and TCMalloc in your build.  The cmake command for this is:

```sh
# machine=sn-fey::thrust
. /usr/share/lmod/lmod/init/sh

INTEL_VERSION=18.0.5
MPI_VERSION=2.1.2
WONTON_VERSION=1.2.8
TANGRAM_VERSION=1.0.4

TPL_PREFIX=/usr/projects/ngc/private

module load intel/${INTEL_VERSION} openmpi/${MPI_VERSION} cmake/3.14.6
mkdir build-thrust
cd build-thrust
cmake \
   -D CMAKE_BUILD_TYPE=Release \
   -D ENABLE_UNIT_TESTS=True \
   -D ENABLE_APP_TESTS=True \
   -D PORTAGE_ENABLE_TANGRAM=True \
   -D WONTON_ROOT:FILEPATH=${TPL_PREFIX}/wonton/${WONTON_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION}-thrust \
   -D TANGRAM_ROOT:FILEPATH=${TPL_PREFIX}/tangram/${TANGRAM_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION}-thrust \
   -D PORTAGE_ENABLE_MPI=True \
   -D PORTAGE_ENABLE_THRUST=True \
   -D PORTAGE_ENABLE_Jali=True \
   -D PORTAGE_ENABLE_FleCSI=False \
   ..
make -j4
ctest -j4 --output-on-failure
```


## Varan

Execute the following from the portage root directory:

```sh
# machine=varan
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh

INTEL_VERSION=18.0.1
MPI_VERSION=2.1.2
WONTON_VERSION=1.2.8
TANGRAM_VERSION=1.0.4

TPL_PREFIX=/usr/local/codes/ngc/private

module load intel/${INTEL_VERSION} openmpi/${MPI_VERSION} cmake/3.14.0
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D PORTAGE_ENABLE_TANGRAM=True \
    -D WONTON_ROOT:FILEPATH=${TPL_PREFIX}/wonton/${WONTON_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION} \
    -D TANGRAM_ROOT:FILEPATH=${TPL_PREFIX}/tangram/${TANGRAM_VERSION}-intel-${INTEL_VERSION}-openmpi-${MPI_VERSION} \
    -D PORTAGE_ENABLE_MPI=True \
	-D PORTAGE_ENABLE_Jali=True \
	-D PORTAGE_ENABLE_FleCSI=False \
    ..
make -j2
ctest -j2 --output-on-failure
```
---


