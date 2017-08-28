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

---

# Details and Requirements

At a minimum, portage requires:
- A C++-11 compatible compiler; regular testing is performed with gcc
  5.3+ and Intel 17+.
- MPI; regular testing is performed with OpenMPI 1.10.3+
- CMake 3.0+
- LAPACKE (3.7.1+)
- Boost (1.53.0+) **__or__** Thrust (1.6.0+)

Distributed parallelism of portage is currently supported through MPI.  Setting
the `ENABLE_MPI` and `ENABLE_MPI_CXX_BINDINGS` CMake flags to `True` (see the
@ref quickstart for more details).

On-node parallelism is exposed through the [Thrust](https://thrust.github.io)
library.  Enabling Thrust within portage requires setting at least two CMake
variables: `ENABLE_THRUST=True` and `THRUST_DIR=<path_to_thrust_directory>`.
Additionally, one can specify the Thrust backend to utilize, with the default
being the OpenMP backend `THRUST_BACKEND=THRUST_DEVICE_SYSTEM_OMP`.  portage
also supports the `THRUST_DEVICE_SYSTEM_TBB` backend, and currently has limited
support for the CUDA backend.  Regular testing happens with Thrust 1.8.

## Obtaining portage

The latest release of [portage](https://github.com/laristra/portage)
lives on GitHub.  portage makes use of git submodules, so it must be
cloned recursively:

```sh
git clone --recursive https://github.com/laristra/portage
```

## Building

portage uses the CMake build system.  Below is an (incomplete) list of
useful CMake flags for building portage.

| CMake flag | Description | Default |
| ---------- |:------------|--------:|
| `ENABLE_MPI` | Build with support for MPI; note most of the code requires this | False |
| `ENABLE_MPI_CXX_BINDINGS` | Tell CMake that the MPI language is C++, not C; this appears to be an OpenMPI thing only | False |
| `ENABLE_THRUST` | Turn on Thrust support for on-node parallelism | False |
| `ENABLE_UNIT_TESTS` | Turn on compilation and test harness of unit tests | False |
| `ENABLE_APP_TESTS` | Turn on compilation and test harness of application tests | False |
| `ENABLE_FleCSI` | Turn on support for FleCSI; _requires C++14-compatible compiler_ | False |
| `Jali_Dir` | Hint location for CMake to find Jali | "" |
| `PC_LAPACKE_NCLUDE_DIRS` | Hint location for CMake to find LAPACKE include files if `pkg_config` fails | "" |
| `PC_LAPACKE_LIBRARY_DIRS` | Hint location for CMake to find LAPACKE library files if `pkg_config` fails | "" |
| `THRUST_DIR` | Directory of the Thrust install | "" |
| `THRUST_BACKEND` | Backend to use for Thrust | `"THRUST_DEVICE_SYSTEM_OMP"` |
