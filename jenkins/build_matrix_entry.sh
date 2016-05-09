#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry.sh <compiler> <flecsi_flag>
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

compiler=$1
flecsi_flag=$2

# set modules and install paths

export NGC=/usr/local/codes/ngc

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
  
if [[ $flecsi_flag == "yes" ]]; then
  flecsi_opts="-D FLECSI_INSTALL_DIR:FILEPATH=$flecsi_install_dir"
else
  flecsi_opts=
fi

# General NGC include directory
ngc_include_dir=/usr/local/codes/ngc/private/include

#git config user.email ""
#git config user.name "Jenkins"

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
  $flecsi_opts \
  ..
make -j2
ctest --output-on-failure
