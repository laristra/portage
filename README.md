# Installation instructions

Below we list copy & paste instructions for several machines. You can easily
adapt them for other machines.

## Darwin

Execute the following from the portage root directory:

```c++
# darwin-fe
module load openmpi/1.10.0-intel_15.0.3 cmake
JALI_INSTALL_PREFIX=/path/to/inst-jali
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

Execute the following from the Jali root directory:

```c++
# barugon
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/15.0.3 openmpi/1.6.5 cmake
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.2-intel-15.0.3-openmpi-1.6.5-rh6
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
make -j2
ctest -j2 --output-on-failure
```

---

If you want to build an app that uses
[FleCSI](https://github.com/losalamos/flecsi), you can link against a built
verison of FleCSI on Varan.  An example is below:

```c++
# barugon::flecsi
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load gcc/5.3.0 openmpi/1.6.5 cmake
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.2-gcc-5.3.0-openmpi-1.6.5-rh6
FLECSI_INSTALL_DIR=/usr/local/codes/ngc/private/flecsi-gcc
mkdir build-flecsi
cd build-flecsi
cmake \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D CMAKE_CXX_COMPILER=`which mpiCC` \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D ENABLE_MPI_CXX_BINDINGS=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D FLECSI_INSTALL_DIR:FILEPATH=$FLECSI_INSTALL_DIR \
    ..
make -j2
ctest -j2 --output-on-failure
```
