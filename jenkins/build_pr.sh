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

# Tag or git commit hash of Jali version to build and use for this PR:
JALI_VERSION=v0.6.0

# Where to find Jali's TPLs:
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-0.6.0-tpl-intel


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
git clean -dfx

# Build Jali

git clone ssh://git@xcp-stash.lanl.gov:7999/ngc/jali.git jali-repo
cd jali-repo
git checkout $JALI_VERSION
JALI_INSTALL_PREFIX=`pwd`/jali-inst
mkdir build
cd build

cmake \
  -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
  -D CMAKE_BUILD_TYPE=Debug \
  -D CMAKE_CXX_FLAGS='-std=c++11' \
  -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
  -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
  -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
  -D ENABLE_MSTK_Mesh:BOOL=TRUE \
  -D ENABLE_STK_Mesh:BOOL=FALSE \
  -D ENABLE_MOAB_Mesh:BOOL=FALSE \
  ..
make -j16
ctest -j16 --output-on-failure
make install


# Build Portage

cd $WORKSPACE
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
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  ..
make -j16
ctest --output-on-failure
