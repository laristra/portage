#!/usr/bin/env bash
#
# Script to build the doxygen docs.  Here, we don't really care about which
# compiler or build options are set.

# Exit on error
set -e

# Echo each command
set -x

JALI_VERSION=jali-0.9.0

# location on XLAN
NGC_DIR=/usr/local/codes/ngc

JALI_INST=${NGC_DIR}/private/${JALI_VERSION}-gcc-5.3.0

git config user.email ""
git config user.name "Jenkins"
git merge origin/master

export SHELL=/bin/sh
export MODULEPATH=""

# Setup modules
. /opt/local/packages/Modules/default/init/sh
module load gcc/5.3.0 openmpi cmake

# the system doxygen and LaTeX are too old; use these instead
export PATH=/usr/local/codes/ngc/home/cmalone/texlive/2016/bin:$PATH
DOXY_EXE=/home/cmalone/code/doxygen/build/bin/doxygen

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
    -D CMAKE_CXX_COMPILER=`which mpiCC` \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D Jali_DIR:FILEPATH=${JALI_INST}/lib \
    -D ENABLE_DOXYGEN=True \
    -D DOXYGEN_EXECUTABLE=$DOXY_EXE \
    ..
make doxygen

# build a PDF copy of the docs
pushd doc/doxygen/latex
make pdf
popd


