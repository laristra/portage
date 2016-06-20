#!/usr/bin/env bash

#Exit on error
set -e
#Echo each command
set -x


# General NGC include directory
NGC_INCLUDE_DIR=/usr/local/codes/ngc/private/include

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load gcc/5.3.0
module load openmpi/1.6.5

JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-0.9.0-gcc-5.3.0

echo $WORKSPACE
cd $WORKSPACE

# Build Portage

mkdir build
cd build

cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_C_FLAGS="-coverage" \
  -D CMAKE_CXX_FLAGS="-coverage" \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Debug \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX//lib \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  -D ENABLE_THRUST=True \
  ..

make -j2
make test

pwd 
#gcovr -r .. -x -e cinch/* -e build/*  -e '.*test.*' -e '.*Test.*'  -e '.*clipper.*' > coverage.xml
gcovr -r .. -x  > coverage.xml

