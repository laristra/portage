# Welcome to portage!   {#mainpage}

portage is a framework that computational physics applications can use
to build a highly customize, hybrid parallel (MPI+X) conservative
remapping library for transfer of field data between meshes.

We aim to provide:
- A modern, modular design - pick and choose your preferred search, intersect, and interpolate method.
- A conservative, intersection-based remap capability on general polytopal
meshes.
- A pick-and-choose framework for search, intersect, and interpolate algorithms.
- High flexibility for application customization.
- Algorithms that take advantage of both distrubuted and on-node parallelism.
- An _Open Source Community_ for these tools!
- Use of client application's native mesh and state (data/field) data structures.

---

# Details and Requirements

At a minimum, portage requires:
- a C++-11 compatible compiler; regular testing is performed with gcc 5+
and Intel 15+.
- CMake 3.0+

Distributed parallelism of portage is currently supported through MPI.  Setting
the `ENABLE_MPI` and `ENABLE_MPI_CXX_BINDINGS` CMake flags to `True` (see the
@ref quickstart for more details).  Regular testing happens with OpenMPI 1.6.5+.

On-node parallelism is exposed through the [Thrust](https://thrust.github.io)
library.  Enabling Thrust within portage requires setting at least two CMake
variables: `ENABLE_THRUST=True` and `THRUST_DIR=<path_to_thrust_directory>`.
Additionally, one can specify the Thrust backend to utilize, with the default
being the OpenMP backend `THRUST_BACKEND=THRUST_DEVICE_SYSTEM_OMP`.  portage
also supports the `THRUST_DEVICE_SYSTEM_TBB` backend, and currently has limited
support for the CUDA backend.  Regular testing happens with Thrust 1.8.
