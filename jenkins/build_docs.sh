#!/usr/bin/env bash
#
# Script to build the doxygen docs.  Here, we don't really care about which
# compiler or build options are set.

# Exit on error
set -e

# Echo each command
set -x

JALI_VERSION=1.0.4
openmpi_version=2.1.2

# location on XLAN
NGC_DIR=/usr/local/codes/ngc

JALI_INST=${NGC_DIR}/private/jali/${JALI_VERSION}-gcc-6.4.0-openmpi-${openmpi_version}
TANGRAM_INSTALL_PREFIX=${NGC_DIR}/private/tangram/0.9.3-gcc-6.4.0-openmpi-${openmpi_version}
XMOF_INSTALL_PREFIX=${NGC_DIR}/private/xmof2d/0.9.4-gcc-6.4.0-openmpi-${openmpi_version}/share/cmake


git config user.email ""
git config user.name "Jenkins"
git merge origin/master

export SHELL=/bin/sh
export MODULEPATH=""

# Setup modules
. /opt/local/packages/Modules/default/init/sh
module load gcc/6.4.0 openmpi/${openmpi_version} cmake

# the system doxygen and LaTeX are too old; use these instead
export PATH=/usr/local/codes/ngc/home/cmalone/texlive/2016/bin/x86_64-linux/:$PATH
DOXY_EXE=/home/cmalone/code/doxygen/build/bin/doxygen

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
    -D Jali_DIR:FILEPATH=${JALI_INST}/lib \
    -D TANGRAM_DIR:FILEPATH=${TANGRAM_INSTALL_PREFIX} \
    -D XMOF2D_DIR:FILEPATH=${XMOF2D_INSTALL_PREFIX} \
    -D ENABLE_DOXYGEN=True \
    -D DOXYGEN_EXECUTABLE=$DOXY_EXE \
    ..
make doxygen

# build a PDF copy of the docs
pushd doc/doxygen/latex
make pdf
popd


