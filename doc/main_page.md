# Welcome to portage!   {#mainpage}

portage is a framework that computational physics applications can use
to build a highly customized, hybrid parallel (MPI+X) conservative
remapping library for transfer of field data between meshes, between
particles, and between meshes and particles.

We aim to provide:
- A modern, modular design - pick and choose your preferred search,
  intersect, and interpolate method
- A conservative, intersection-based remap capability on general
  polytopal meshes and general particle shape functions
- A pick-and-choose framework for search, intersect, and interpolate
  algorithms
- High flexibility for application customization
- Algorithms that take advantage of both distrubuted and on-node parallelism
- An _Open Source Community_ for these tools!
- Use of client application's native mesh/particle and state
  (data/field) data structures

See the [Concepts](@ref concepts) page for a high-level discussion of
the methods used within portage.

See the [Example Use](@ref example) page for a simple example of
hooking portage up to a mesh and state manager.

---

# Details and Requirements

At a minimum, portage requires:
- A C++-11 compatible compiler; regular testing is performed with GCC
  5.3+ and Intel 17+.
- CMake 3.0+
- LAPACKE (3.7.1+)
- Boost (1.53.0+) **or** Thrust (1.6.0+)

Distributed parallelism of portage is currently supported through MPI;
regular testing is performed with OpenMPI 1.10.3+ .  Most application
tests and all of the units tests are currently only built if MPI is
used.  MPI is enabled in portage by setting the CMake variable
`ENABLE_MPI=True`.  In addition, you'll likely need to tell CMake to
use the MPI-wrapped compiler by setting something like

~~~
CMAKE_C_COMPILER=`which mpicc`
CMAKE_CXX_COMPILER=`which mpiCC`
~~~

On-node parallelism is exposed through
the [Thrust](https://thrust.github.io) library.  Enabling Thrust
within portage requires setting at least two CMake variables:
`ENABLE_THRUST=True` and `THRUST_DIR=<path_to_thrust_directory>`.
Additionally, one can specify the Thrust backend to utilize, with the
default being the OpenMP backend
`THRUST_BACKEND=THRUST_DEVICE_SYSTEM_OMP`.  portage also supports the
`THRUST_DEVICE_SYSTEM_TBB` backend.  Regular testing happens with
Thrust 1.8.

## Obtaining portage

The latest release of [portage](https://github.com/laristra/portage)
lives on GitHub.  portage makes use of git submodules, so it must be
cloned recursively:

```sh
git clone --recursive https://github.com/laristra/portage
```

## Building

portage uses the CMake build system.  In the simplest case where you
want to build a serial version of the code, and CMake knows where to
find your Boost and LAPACKE installations, one can do

```sh
portage/ $ mkdir build
portage/ $ cd build
portage/build/ $ cmake ..
portage/build/ $ make
```

This will build a serial version of the code into a library (without
any tests).  A more complete build with MPI, Thrust (for on-node
parallelism), unit and application test support, documentation
support, and support for both [Jali](https://github.com/lanl/jali) and
the Burton [specialization](https://github.com/laristra/flecsi-sp) of
the [FleCSI](https://github.com/laristra/flecsi) library would look
like:

~~~sh
portage/ $ mkdir build
portage/ $ cd build
portage/build/ $ cmake -DCMAKE_C_COMPILER=`which mpicc` \
                       -DCMAKE_CXX_COMPILER=`which mpiCC` \
					   -DENABLE_APP_TESTS=True -DENABLE_UNIT_TESTS=True \
					   -DENABLE_MPI=True \
					   -DENABLE_THRUST=True -DTHRUST_DIR=/path/to/thrust/include/directory \
					   -DJali_DIR=path/ to/Jali/lib \
					   -DENABLE_FleCSI=True \
					   -DCMAKE_PREFIX_PATH="/path/to/FleCSI/install;/path/to/FleCSI-sp/install" \
					   -DENABLE_DOXYGEN=True \
					   -DPC_LAPACKE_NCLUDE_DIRS=/path/to/LAPACKE/include \
					   -DPC_LAPACKE_LIBRARY_DIRS=/path/to/LAPACKE/install \
					   ..
portage/build/ $ make           # builds the library and tests
portage/build/ $ make test      # runs the tests
portage/build/ $ make doxygen   # builds this HTML and a PDF form of the documentation
portage/build/ $ make install   # installs the portage library and headers into CMAKE_INSTALL_PREFIX
~~~

## Useful CMake Flags
Below is a non-exhaustive list of useful CMake flags for building
portage.

| CMake flag:type | Description | Default |
|:----------|:------------|:--------|
| `CMAKE_BUILD_TYPE:STRING`| `Debug` or optimized `Release` build | `Debug` |
| `CMAKE_INSTALL_PREFIX:PATH` | Location for the portage library and headers to be installed | `/usr/local` |
| `CMAKE_PREFIX_PATH:PATH` | Locations where CMake can look for packages; needs to be set to the FleCSI and FleCSI-SP locations if using FleCSI | "" |
| `ENABLE_APP_TESTS:BOOL` | Turn on compilation and test harness of application tests | `False` |
| `ENABLE_DOXYGEN:BOOL` | Create a target to build this documentation | `False` |
| `ENABLE_FleCSI:BOOL` | Turn on support for the FleCSI Burton specialization; must set `CMAKE_PREFIX_PATH` to a location where _both_ FleCSI and FleCSI-SP can be found | `False` |
| `ENABLE_MPI:BOOL` | Build with support for MPI; | `False` |
| `ENABLE_TCMALLOC:BOOL` | Build with support for TCMalloc | `False` |
| `ENABLE_THRUST:BOOL` | Turn on Thrust support for on-node parallelism | `False` |
| `ENABLE_UNIT_TESTS:BOOL` | Turn on compilation and test harness of unit tests | `False` |
| `ENABLE_FleCSI:BOOL` | Turn on support for FleCSI; _requires C++14-compatible compiler_ | `False` |
| `Jali_DIR:PATH` | Hint location for CMake to find Jali | "" |
| `PC_LAPACKE_NCLUDE_DIRS:PATH` | Hint location for CMake to find LAPACKE include files if `pkg_config` fails | "" |
| `PC_LAPACKE_LIBRARY_DIRS:PATH` | Hint location for CMake to find LAPACKE library files if `pkg_config` fails | "" |
| `TCMALLOC_LIB:PATH` | The TCMalloc library to use | `${HOME}` |
| `THRUST_DIR:PATH` | Directory of the Thrust install | "" |
| `THRUST_BACKEND:STRING` | Backend to use for Thrust | `"THRUST_DEVICE_SYSTEM_OMP"` |
