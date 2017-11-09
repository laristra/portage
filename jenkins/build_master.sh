#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     bash $WORKSPACE/jenkins/build_master.sh
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

# Tag or git commit hash of Jali version to build and use for this PR:
JALI_VERSION=v0.9.8
openmpi_version=1.10.5

# Where to find Jali's TPLs:
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-tpl/1.0.9-intel-17.0.1-openmpi-${openmpi_version}


git config user.email ""
git config user.name "Jenkins"
git merge origin/master

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/17.0.1
module load openmpi/${openmpi_version}
module load cmake

echo $WORKSPACE
cd $WORKSPACE
git clean -dfx

# Build Jali

git clone ssh://git@xcp-stash.lanl.gov:7999/laristra/jali.git jali-repo
cd jali-repo
git checkout $JALI_VERSION
JALI_INSTALL_PREFIX=`pwd`/jali-inst
mkdir build
cd build

cmake \
  -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS='-std=c++11' \
  -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
  -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
  -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
  -D ENABLE_MSTK_Mesh:BOOL=TRUE \
  -D ENABLE_STK_Mesh:BOOL=FALSE \
  -D ENABLE_MOAB_Mesh:BOOL=FALSE \
  ..
make -j2
ctest -j2 --output-on-failure
make install


# Build Portage

cd $WORKSPACE
mkdir build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  ..
make -j2
ctest --output-on-failure
