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

# special case for README builds
if [[ $build_type == "readme" ]]; then
  # Put a couple of settings in place to generate test output even if
  # the README doesn't ask for it.
  export CTEST_OUTPUT_ON_FAILURE=1
  CACHE_OPTIONS="-D ENABLE_JENKINS_OUTPUT=True"
  sed "s/^ *cmake/& $CACHE_OPTIONS/g" $WORKSPACE/README.md >$WORKSPACE/README.md.1
  python2 $WORKSPACE/jenkins/parseREADME.py \
      $WORKSPACE/README.md.1 \
      $WORKSPACE \
      varan
  exit
fi

# set modules and install paths

jali_version=1.0.0
tangram_version=34fc65e27e0
xmof2d_version=529f2dcdbe4
lapack_version=3.8.0
 

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include

# compiler-specific settings
if [[ $compiler == "intel" ]]; then
  intel_version=18.0.1
  cxxmodule=intel/${intel_version}
  # openmpi version that libs were built against
  openmpi_version=2.1.2
  # openmpi module for compiling and linking
  mpi_module=openmpi/2.1.2
  jali_install_dir=$NGC/private/jali/${jali_version}-intel-${intel_version}-openmpi-${openmpi_version}
  tangram_install_dir=$NGC/private/tangram/${tangram_version}-intel-${intel_version}-openmpi-${openmpi_version}
  xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}-intel-${intel_version}
  lapacke_dir=$NGC/private/lapack/${lapack_version}-patched-intel-${intel_version}
elif [[ $compiler == "gcc" ]]; then
  gcc_version=6.4.0
  cxxmodule=gcc/${gcc_version}
  # openmpi version that libs were built against
  openmpi_version=2.1.2
  # openmpi module for compiling and linking
  mpi_module=openmpi/2.1.2
  jali_install_dir=$NGC/private/jali/${jali_version}-gcc-${gcc_version}-openmpi-${openmpi_version}
  flecsi_install_prefix=$NGC/private/flecsi/374b56b-gcc-6.4.0
  flecsisp_install_prefix=$NGC/private/flecsi-sp/e78c594-gcc-6.4.0
  tangram_install_dir=$NGC/private/tangram/${tangram_version}-gcc-${gcc_version}-openmpi-${openmpi_version}
  xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}-gcc-${gcc_version}
  lapacke_dir=$NGC/private/lapack/${lapack_version}-patched-gcc-${gcc_version}
fi
  
cmake_build_type=Release
extra_flags=
jali_flags="-D Jali_DIR:FILEPATH=$jali_install_dir/lib"
tangram_flags="-D TANGRAM_DIR:FILEPATH=$tangram_install_dir"
xmof2d_flags="-D XMOF2D_DIR:FILEPATH=$xmof2d_install_dir/share/cmake"
mpi_flags="-D ENABLE_MPI=True"
lapacke_flags="-D LAPACKE_DIR:FILEPATH=$lapacke_dir"

if [[ $build_type == "debug" ]]; then
  cmake_build_type=Debug
elif [[ $build_type == "serial" ]]; then
  mpi_flags=
  # jali is not available in serial
  jali_flags=
elif [[ $build_type == "thrust" ]]; then
  extra_flags="-D ENABLE_THRUST=True"
elif [[ $build_type == "flecsi" ]]; then
  extra_flags="-D CMAKE_PREFIX_PATH='$flecsi_install_prefix;$flecsisp_install_prefix' \
               -D ENABLE_FleCSI=True"
elif [[ $build_type == "coverage" ]]; then
  extra_flags="-D CMAKE_C_FLAGS='-coverage' \
               -D CMAKE_CXX_FLAGS='-coverage'"
  cmake_build_type=Debug
  export PATH=$NGC/private/bin:${PATH}
fi

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
module load cmake # 3.0 or higher is required

if [[ -n "$mpi_flags" ]] ; then
  module load ${mpi_module}
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
  $tangram_flags \
  $xmof2d_flags \
  $lapacke_flags \
  $extra_flags \
  ..
make -j2
ctest --output-on-failure

if [[ $build_type == "coverage" ]]; then
  gcovr -r .. -x  >coverage.xml
fi

