#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry_hpc.sh <compiler> <build_type>
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
	    sn-fey
    exit

fi

# set modules and install paths

openmpi_version=2.1.2
wonton_version=new-cmake
tangram_version=new-cmake

export NGC=/usr/projects/ngc
ngc_include_dir=$NGC/private/include
ngc_tpl_dir=$NGC/private

# compiler-specific settings
if [[ $compiler == "intel18" ]]; then

    compiler_type=intel
    compiler_version=18.0.5
    cxxmodule=intel/${compiler_version}

elif [[ $compiler == "gcc6" ]]; then

    compiler_type=gcc
    compiler_version=6.4.0
    cxxmodule=gcc/${compiler_version}

elif [[ $compiler == "gcc7" ]]; then

    compiler_type=gcc
    compiler_version=7.4.0
    cxxmodule=gcc/${compiler_version}

fi

wonton_install_dir_base=${ngc_tpl_dir}/wonton/${wonton_version}-${compiler_type}-${compiler_version}
wonton_install_suffix="-openmpi-${openmpi_version}"
tangram_install_dir_base=${ngc_tpl_dir}/tangram/${tangram_version}-${compiler_type}-${compiler_version}
tangram_install_suffix="-openmpi-${openmpi_version}"

cmake_build_type=Release
extra_flags=
mpi_flags="-D PORTAGE_ENABLE_MPI=True"

if [[ $build_type == "debug" ]]; then

    cmake_build_type=Debug

elif [[ $build_type == "serial" ]]; then

    mpi_flags=
    tangram_install_suffix="-nompi"

elif [[ $build_type == "thrust" ]]; then

    extra_flags="-D PORTAGE_ENABLE_THRUST=True"
    wonton_install_suffix="${wonton_install_suffix}-thrust"
    tangram_install_suffix="${tangram_install_suffix}-thrust"

fi

tangram_flags="-D PORTAGE_ENABLE_TANGRAM=True \
-D TANGRAM_DIR:FILEPATH=${tangram_install_dir_base}${tangram_install_suffix}"

export SHELL=/bin/sh

. /usr/share/lmod/lmod/init/sh
module load ${cxxmodule}
module load cmake/3.14.6 # 3.13 or higher is required

if [[ -n "$mpi_flags" ]] ; then
  module load openmpi/${openmpi_version}
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
  $mpi_flags \
  $wonton_flags \
  $tangram_flags \
  $extra_flags \
  ..
make -j2
ctest --output-on-failure

