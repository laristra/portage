#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     bash $WORKSPACE/jenkins/build_pr.sh
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

git config user.email ""
git config user.name "Jenkins"
git merge origin/master

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/15.0.3
module load openmpi/1.6.5

echo $WORKSPACE
cd $WORKSPACE
rm -rf build
mkdir build
cd build

cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Debug \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=/usr/local/codes/ngc/private/jali-0.6.0-intel/lib \
  ..
make
ctest --output-on-failure
