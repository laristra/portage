#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry.sh <compiler> <build_type>
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

compiler=$1
build_type=$2

# set modules and install paths

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include

# compiler-specific settings
if [[ $compiler == "intel" ]]; then
  cxxmodule=intel/15.0.3
  jali_install_dir=$NGC/private/jali-0.7.1-intel
elif [[ $compiler == "gcc" ]]; then
  cxxmodule=gcc/4.9.2
  jali_install_dir=$NGC/private/jali-0.7.1-gcc
elif [[ $compiler == "gcc53" ]]; then
  cxxmodule=gcc/5.3.0
  jali_install_dir=$NGC/private/jali-0.7.1-gcc53
  flecsi_install_dir=$NGC/private/flecsi-gcc
fi
  
extra_flags=
if [[ $build_type == "thrust" ]]; then
  extra_flags="-D ENABLE_THRUST=True"
elif [[ $build_type == "flecsi" ]]; then
  extra_flags="-D FLECSI_INSTALL_DIR:FILEPATH=$flecsi_install_dir"
elif [[ $build_type == "coverage" ]]; then
  extra_flags="-D CMAKE_C_FLAGS='-coverage' \
               -D CMAKE_CXX_FLAGS='-coverage'"
fi

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
module load openmpi/1.6.5

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Debug \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  -D Jali_DIR:FILEPATH=$jali_install_dir/lib \
  -D NGC_INCLUDE_DIR:FILEPATH=$ngc_include_dir \
  $extra_flags \
  ..
make -j2
ctest --output-on-failure

if [[ $build_type == "coverage" ]]; then
  gcovr -r .. -x  >coverage.xml
fi

