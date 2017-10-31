#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry.sh <compiler> <build_type>
#
# The exit code determines if the test succeeded or failed.
# Note that the environment variable WORKSPACE must be set (Jenkins
# will do this automatically).

# Exit on error
set -e
# Echo each command
set -x

compiler=$1
build_type=$2

# set modules and install paths

jali_version=0.9.8
openmpi_version=1.10.5

export NGC=/usr/projects/ngc
ngc_include_dir=$NGC/private/include

# compiler-specific settings
if [[ $compiler == "intel" ]]; then
  cxxmodule=intel/17.0.1
  jali_install_dir=$NGC/private/jali/${jali_version}-intel-17.0.1-openmpi-${openmpi_version}
elif [[ $compiler == "gcc" ]]; then
  cxxmodule=gcc/5.3.0
  jali_install_dir=$NGC/private/jali/${jali_version}-gcc-5.3.0-openmpi-${openmpi_version}
fi
  
cmake_build_type=Release
extra_flags=
jali_flags="-D Jali_DIR:FILEPATH=$jali_install_dir/lib"
mpi_flags="-D ENABLE_MPI=True"

if [[ $build_type == "debug" ]]; then
  cmake_build_type=Debug
elif [[ $build_type == "serial" ]]; then
  mpi_flags=
  # jali is not available in serial
  jali_flags=
elif [[ $build_type == "thrust" ]]; then
  extra_flags="-D ENABLE_THRUST=True"
fi

export SHELL=/bin/sh

. /usr/share/lmod/lmod/init/sh
module load $cxxmodule
module load cmake # 3.0 or higher is required

if [[ -n "$mpi_flags" ]] ; then
  module load openmpi/${openmpi_version}
  mpi_flags+=" -D CMAKE_C_COMPILER=`which mpicc` \
               -D CMAKE_CXX_COMPILER=`which mpiCC`"
fi

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=$cmake_build_type \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_APP_TESTS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D NGC_INCLUDE_DIR:FILEPATH=$ngc_include_dir \
  $mpi_flags \
  $jali_flags \
  $extra_flags \
  ..
make -j2
ctest --output-on-failure

