#!/usr/bin/env bash

#Exit on error
set -e
#Echo each command
set -x


#General NGC directory
NGC=/usr/local/codes/ngc

# General NGC include directory
NGC_INCLUDE_DIR=$NGC/private/include

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load gcc/5.3.0
module load openmpi/1.10.3

JALI_INSTALL_PREFIX=$NGC/private/jali/0.9.8-gcc-5.3.0-openmpi-1.10.3
TANGRAM_INSTALL_PREFIX=$NGC/private/tangram/133c1db580f-gcc-5.3.0-openmpi-1.10.3
XMOF2D_INSTALL_PREFIX=$NGC/private/xmof2d/0.9-gcc-5.3.0-openmpi-1.10.3

echo $WORKSPACE
cd $WORKSPACE

# Build Portage

mkdir build
cd build

cmake \
  -D CMAKE_C_FLAGS="-coverage" \
  -D CMAKE_CXX_FLAGS="-coverage" \
  -D CMAKE_BUILD_TYPE=Debug \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX//lib \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  -D TANGRAM_DIR:FILEPATH=$TANGRAM_INSTALL_PREFIX \
  -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
  -D ENABLE_THRUST=True \
  ..

make -j2
make test

pwd 
#gcovr -r .. -x -e cinch/* -e build/*  -e '.*test.*' -e '.*Test.*'  -e '.*clipper.*' > coverage.xml
gcovr -r .. -x  > coverage.xml

