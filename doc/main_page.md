# Introduction {#mainpage} 
<a name="Introduction"></a>

Portage is a framework for building a **highly customized, hybrid
parallel (MPI+X) remapping library** for the transfer of
field data between *meshes* or [*particle
swarms*](concepts.html#meshfree remap) in computational physics applications.

Application developers can:
- **use one of the included drivers** with a mix of available and custom
components to readily deploy a powerful remapping capability into
their application, or
- **write a custom remapping driver** and use it with a mix of available
and custom components to create a remapping capability uniquely
tailored to their application needs.

<img src="portage-tangram-diagram.png" alt="Portage design" width="100%">
<br>

The drivers furnished with Portage provide:
- [conservative remapping of multi-material fields between general
  polygonal/polyhedral meshes](concepts.html#mesh-mesh remap).
- [higher-order, non-conservative interpolation between particle
swarms](concepts.html#meshfree remap) as well as between meshes and
particle swarms.
- a modern design templated on the major components - mix and match from
  the furnished suite or use a custom component.
- direct (no-copy) use of client application's native mesh/particle and 
  field data structures whenever possible (see [distributed remap](@ref distributed_concepts)).
- built-in distributed and on-node parallelism even with custom
components (see [scaling](performance.html#scaling) results).

They also have the following known limitations:
- Only remapping of scalar fields is supported; vectors and tensors
must be remapped component by component.
- Multi-material remap is not yet supported for distributed meshes.
- Conservative remapping of compound quantities such as
momentum is not supported when they are not the primary fields.
- Field fixup in the presence of mismatched mesh boundaries is not
implemented for multi-material fields.
- Remapping involving particles is currently non-conservative.

<br>

See the [Concepts](@ref concepts) page for a high-level discussion of
the methods used within Portage.

See the [Example Use](@ref example) page for a simple example of
hooking Portage up to a mesh and state manager.

---

<a name="Requirements"></a>
# Requirements

At a minimum, Portage requires:
- A C++-11 compatible compiler; regular testing is performed with GCC
  6.3+ and Intel 17+.
- [CMake](https://cmake.org) 3.8+
- [LAPACKE](https://https://github.com/Reference-LAPACK/lapack/tree/master/LAPACKE) 3.8.0+
- [Boost](https://www.boost.org) 1.68.0+ **or** [Thrust](https://thrust.github.io/) 1.8.1+

Portage depends on the [Wonton](https://github.com/laristra/wonton)
library to provide mesh/state wrapper interfaces and some common
defintions and functionality. Portage is known to work with Wonton
mesh/state wrappers for version 1.0.0 of the
[Jali](https://github.com/lanl/jali) and \#e78c594 of
[flecsi-sp](https://github.com/laristra/flecsi-sp) (along with
\#374b56b of [FleCSI](https://github.com/laristra/flecsi)).

The supplied mesh-mesh remap driver of Portage also requires the use
of the [R3D](https://www.github.com/devonmpowell/r3d) library for
intersection of polyhedra.

**Distributed parallelism** of Portage is currently supported through MPI;
regular testing is performed with OpenMPI 1.10.3+ .  Most application
tests and all of the units tests are currently only built if MPI is
used.  MPI is enabled in Portage by setting the CMake variable
`ENABLE_MPI=True`.

**On-node parallelism** is enabled through the
[Thrust](https://thrust.github.io) library. The default on-node
parallelism mechanism is OpenMP (YMMV for TBB, no support
yet for CUDA). Regular testing happens with Thrust 1.8. **If you turn
on Thrust for multi-threading-enabled executables, the team strongly
recommends linking to the TCMalloc library available in [Google
Performance Tools](https://github.com/gperftools/gperftools) to see the expected scaling.**

**Multi-material remapping** requires use of the
[Tangram](https://github.com/laristra/tangram) library for interface
reconstruction and optionally, the
[XMOF2D](https://github.com/laristra/xmof2d) library. Portage has been
tested with version 0.9.1 of both Tangram and the X-MOF library.

<a name="Downloading"></a>
## Downloading

The latest release of [Portage](https://github.com/laristra/Portage)
lives on GitHub.  Portage makes use of git submodules, so it must be
cloned recursively:

```sh
git clone --recursive https://github.com/laristra/Portage
```

<a name="Building"></a>
## Building

Portage uses the [CMake](https://cmake.org) build system. In the
simplest case where you want to build a serial version of the code,
and CMake knows where to find your Boost and LAPACKE installations,
one can do

```sh
Portage/ $ mkdir build
Portage/ $ cd build
Portage/build/ $ cmake ..
Portage/build/ $ make
```

This will build a serial version of the code into a library (without
any tests).  A more complete build with MPI, Thrust/OpenMP, unit and
application test support, documentation support, and support for both
[Jali](https://github.com/lanl/jali) and the Burton
[specialization](https://github.com/laristra/flecsi-sp) of the
[FleCSI](https://github.com/laristra/flecsi) library would look like:

~~~sh
Portage/ $ mkdir build
Portage/ $ cd build
Portage/build/ $ cmake -DENABLE_APP_TESTS=True \
                       -DENABLE_UNIT_TESTS=True \
                       -DENABLE_MPI=True \
                       -DENABLE_THRUST=True -DTHRUST_DIR=/path/to/thrust/include/directory \
                       -DENABLE_TCMALLOC=True -DTCMALLOC_LIB=path/to/Optional/TCMalloc/lib \
                       -DJali_DIR=path/to/Optional/Jali/lib \
                       -DENABLE_FleCSI=True \
                       -DCMAKE_PREFIX_PATH="/path/to/Optional/FleCSI/install;/path/to/FleCSI-sp/install" \
                       -DENABLE_DOXYGEN=True \
                       -DLAPACKE_DIR=/path/to/LAPACKE
					   ..
Portage/build/ $ make           # builds the library and tests
Portage/build/ $ make test      # runs the tests
Portage/build/ $ make doxygen   # builds this HTML and a PDF form of the documentation
Portage/build/ $ make install   # installs the Portage library and headers into CMAKE_INSTALL_PREFIX
~~~

<a name="Useful CMake Flags"></a>
## Useful CMake Flags
Below is a non-exhaustive list of useful CMake flags for building
Portage.

| CMake flag:type | Description | Default |
|:----------|:------------|:--------|
| `CMAKE_BUILD_TYPE:STRING`| `Debug` or optimized `Release` build | `Debug` |
| `CMAKE_INSTALL_PREFIX:PATH` | Location for the Portage library and headers to be installed | `/usr/local` |
| `CMAKE_PREFIX_PATH:PATH` | Locations where CMake can look for packages; needs to be set to the FleCSI and FleCSI-SP locations if using FleCSI | NO_DEFAULT |
| `ENABLE_APP_TESTS:BOOL` | Turn on compilation and test harness of application tests | `False` |
| `ENABLE_DOXYGEN:BOOL` | Create a target to build this documentation | `False` |
| `ENABLE_FleCSI:BOOL` | Turn on support for the FleCSI Burton specialization; must set `CMAKE_PREFIX_PATH` to a location where _both_ FleCSI and FleCSI-SP can be found. Both FleCSI packages are under constant development. | `False` |
| `ENABLE_MPI:BOOL` | Build with support for MPI | `False` |
| `ENABLE_TCMALLOC:BOOL` | Build with support for TCMalloc | `False` |
| `ENABLE_THRUST:BOOL` | Turn on Thrust support for on-node parallelism | `False` |
| `ENABLE_UNIT_TESTS:BOOL` | Turn on compilation and test harness of unit tests | `False` |
| `ENABLE_FleCSI:BOOL` | Turn on support for FleCSI; _requires C++14-compatible compiler_ | `False` |
| `Jali_DIR:PATH` | Hint location for CMake to find Jali | NO_DEFAULT |
| `LAPACKE_DIR:PATH` | Hint location for CMake to find LAPACKE include and library files | NO_DEFAULT |
| `TCMALLOC_LIB:PATH` | The TCMalloc library to use (recommended if enabling Thrust) | NO_DEFAULT |
| `THRUST_DIR:PATH` | Directory of the Thrust install | NO_DEFAULT |
| `THRUST_BACKEND:STRING` | Backend to use for Thrust | `"THRUST_DEVICE_SYSTEM_OMP"` |
