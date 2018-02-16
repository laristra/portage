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

# Put a couple of settings in place to generate test output even if
# the README doesn't ask for it.
export CTEST_OUTPUT_ON_FAILURE=1
CACHE_OPTIONS="-D ENABLE_JENKINS_OUTPUT=True"
sed "s/^ *cmake/& $CACHE_OPTIONS/g" $WORKSPACE/README.md >$WORKSPACE/README.md.1

# Run build/test commands from README
python $WORKSPACE/jenkins/parseREADME.py $WORKSPACE/README.md.1 $WORKSPACE

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
JALI_VERSION=0.9.8

# Where to find Jali's TPLs:
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-tpl/1.0.9-intel-17.0.1-openmpi-1.10.5

#General NGC directory
NGC=/usr/local/codes/ngc

# General NGC include directory
NGC_INCLUDE_DIR=$NGC/private/include

git config user.email ""
git config user.name "Jenkins"

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/17.0.1
module load openmpi/1.10.5

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


TANGRAM_INSTALL_PREFIX=$NGC/private/tangram/475b813919f-intel-17.0.1-openmpi-1.10.5
XMOF_INSTALL_PREFIX=$NGC/private/xmof2d/0.9-intel-17.0.1-openmpi-1.10.5/share/cmake

# Build Portage with Thrust

cd $WORKSPACE
mkdir build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  -D TANGRAM_DIR:FILEPATH=$TANGRAM_INSTALL_PREFIX \
  -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
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
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  -D TANGRAM_DIR:FILEPATH=$TANGRAM_INSTALL_PREFIX \
  -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  -D ENABLE_THRUST=False \
  ..
make -j2
ctest --output-on-failure
