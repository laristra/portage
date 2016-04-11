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

python $WORKSPACE/jenkins/parseREADME.py $WORKSPACE/README.md $WORKSPACE

# TEMPORARY FIX:
# Exit at this point.  We've already tested the PR on varan against 
# a fixed Jali release, so we know the PR is good.
exit

# The remaining part of this script currently runs on vulpix, and
# our Jali TPL 1.0.5 install does not work there since vulpix and
# varan have different environments.
# TODO:  Either modify the below to run on varan, or delete it if we
#        decide we don't need it.  Need a team discussion to decide.

#
# Tag or git commit hash of Jali version to build and use for this PR:
JALI_VERSION=239d3f314e1ebc73fa0a16dbea7a0156a5e06544

# Where to find Jali's TPLs:
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-1.0.2-tpl-intel

# General NGC include directory
NGC_INCLUDE_DIR=/usr/local/codes/ngc/private/include

git config user.email ""
git config user.name "Jenkins"

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
make -j2
ctest -j2 --output-on-failure
make install


# Build Portage with Thrust

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
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  -D ENABLE_THRUST=True \
  ..
make -j2
ctest --output-on-failure
for test_name in intersectClipper test_matfuncs; do
  valgrind --gen-suppressions=all --suppressions=../jenkins/portage_valgrind.supp --error-exitcode=1 --leak-check=full test/$test_name
done


# Build Portage without Thrust

cd $WORKSPACE
mkdir build-nothrust
cd build-nothrust

cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Debug \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  -D ENABLE_THRUST=False \
  ..
make -j2
ctest --output-on-failure
