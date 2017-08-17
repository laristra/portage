

# Installation instructions

Below we list copy & paste instructions for several machines. You can easily
adapt them for other machines.

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
    -D ENABLE_MPI_CXX_BINDINGS=True \
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
    -D ENABLE_MPI_CXX_BINDINGS=True \
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
LAPACKE_INCLUDE_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-intel-17.0.1/LAPACKE/include
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
    -D ENABLE_MPI_CXX_BINDINGS=True \
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
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.8-gcc-5.3.0-openmpi-1.10.3
LAPACKE_INCLUDE_DIR=/usr/local/codes/ngc/private/lapack/lapack-3.7.1-gcc-5.3.0/LAPACKE/include
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
    -D ENABLE_MPI_CXX_BINDINGS=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D PC_LAPACKE_NCLUDE_DIRS=$LAPACKE_INCLUDE_DIR \
    -D PC_LAPACKE_LIBRARY_DIRS=$LAPACKE_LIBRARY_DIR \
    ..
make -j2
ctest -j2 --output-on-failure
```
