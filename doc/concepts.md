# Quickstart Guide      {#concepts}

## Obtaining portage
portage is maintained within a git [repository](https://github.com/laristra/portage)
on GitHub that contains a couple of git submodules.  To fully checkout the code
and its submodules into a new directory called `portage`, execute the following
~~~
$ git clone --recursive https://github.com/laristra/portage.git
~~~

## Building the Code
portage uses [CMake](https://cmake.org) to configure its Makefiles, and
currently requires CMake version 3.0 or greater.  Additionally, portage does not
allow for in-source builds, which implies the first step is to create a build
directory
~~~
$ cd portage
$ mkdir build
$ cd build
~~~
The next step is to invoke CMake to build our Makefiles.  On the command line,
one can enter (at a minimum)
~~~
$ cmake \
      -D CMAKE_C_COMPILER=<path_to_your_C_compiler> \
	  -D CMAKE_CXX_COMPILER=<path_to_your_C++_compiler> \
	  -D ENABLE_UNIT_TESTS=True \
	  -D ENABLE_APP_TESTS=True \
	  ..
~~~
The flags `ENABLE_UNIT_TESTS` and `ENABLE_APP_TESTS` are not strictly
_required_ but will provide for running sample tests below.  As a
concrete example, our local systems use modulefiles to load
environments, and we would typically do something like the following to
build with MPI
~~~
$ module load intel/15.0.3 openmpi/1.6.5 cmake
$ cmake \
      -D CMAKE_C_COMPILER=`which mpicc` \
	  -D CMAKE_CXX_COMPILER=`which mpiCC` \
	  -D ENABLE_MPI=True \
	  -D ENABLE_MPI_CXX_BINDINGS=True \
	  -D ENABLE_UNIT_TESTS=True \
	  -D ENABLE_APP_TESTS=True \
	  ..
~~~
Alternatively, you could run the interactive `ccmake` command, which exposes
many more variables that can be set within the CMake build step.

After one of the above commands, CMake will generate several Makefiles.  Simply
running
~~~
$ make
~~~
will compile the source code.

## Running Sample Tests
After the code is compiled, running
~~~
$ make test
~~~
will run several tests.  The number and type of tests actually run depends on
the options turned on during the CMake configuration step.  For example,
`ENABLE_APP_TESTS` will run the sample applications of remapping that are
declared within the `portage/test` directory; the actual source for the app
tests live in the `portage/app` directory.

Under the hood, we use Google Test to configure our testing framework.
Additionally, we use CMake to conditionally turn on specific tests based on
settings passed at configure time.  For example, if `ENABLE_MPI=True` at
configure time, then we will include both unit and app tests that utilize
MPI; otherwise they are ignored.  *NOTE* If you request MPI and run the tests,
make sure that you actually have access to MPI --- on some clusters, the
frontend or login nodes don't actually have the ability to call `MPI_INIT`,
and you may get obscure errors if you try.

## What's next?
The apps directory is a good place to look for _how_ one uses the methods
within the portage library to actually perform a remap.  Out of the box,
portage provides wrappers to three mesh and state manager frameworks:

1. `Simple_Mesh` and `Simple_State` are examples of a very simple mesh and
state manager framework that is included with portage.  It is not intended to
be production quality, but to demonstrate how one might hook their own mesh and
state manager framework into portage.  A detailed description of `Simple_Mesh`
and `Simple_State` and how they are used is shown on the @ref example page.

2. [`FleCSI`](https://github.com/losalamos/flecsi) is a modern, flexible,
computational framework designed for next-generation architectures.
In particular, portage currently uses the _Burton Mesh_ specialization.  FleCSI
is undergoing change at fast pace, and portage's support of FleCSI features is
a little behind. *portage is currently known to be consistent with hash*
`170b54a` of FleCSI.  If you wish to run portage with FleCSI, please checkout
that version of FleCSI and build it.  Then to enable within portage, set the
`FLECSI_INSTALL_DIR:FILEPATH=<path_to_your_flecsi_installation>` CMake variable.
Note that FleCSI requires a C++-14 compatible compiler.  We are in the process
of updating portage's FleCSI support to be more up-to-date.

3. `Jali` is a stripped-down version of the mesh infrastructure of the
[Amanzi](https://software.lanl.gov/ascem/amanzi/) code.  It is currently not
open-source and is mainly used internally for testing.
