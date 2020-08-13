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

echo "--------------------------------------------------------------"
echo "Running configuration $COMPILER $BUILD_TYPE on `hostname`"
echo "--------------------------------------------------------------"

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

wonton_version=1.2.2
tangram_version=1.0.1

export NGC=/usr/projects/ngc
ngc_include_dir=$NGC/private/include
ngc_tpl_dir=$NGC/private

# compiler-specific settings
if [[ $compiler == "intel18" ]]; then

    compiler_version=18.0.5
    cxxmodule=intel/${compiler_version}
    compiler_suffix="-intel-${compiler_version}"

    openmpi_version=2.1.2
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"
    
elif [[ $compiler =~ "gcc" ]]; then

    openmpi_version=2.1.2
    if [[ $compiler == "gcc6" ]]; then
	compiler_version=6.4.0
    elif [[ $compiler == "gcc7" ]]; then
	compiler_version=7.4.0
    fi
    
    cxxmodule=gcc/${compiler_version}
    compiler_suffix="-gcc-${compiler_version}"

    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"

fi

mpi_flags="-D PORTAGE_ENABLE_MPI=True"
if [[ $build_type == "serial" ]]; then
    mpi_flags=
    mpi_suffix=
fi

cmake_build_type=Release
if [[ $build_type == "debug" ]]; then
    cmake_build_type=Debug
fi

thrust_flags=
thrust_suffix=
if [[ $build_type == "thrust" ]]; then
    thrust_flags="-D PORTAGE_ENABLE_THRUST=True"
    thrust_suffix="-thrust"
fi

cov_flags=
if [[ $build_type == "coverage" ]]; then
    cov_flags="-D CMAKE_C_FLAGS='-coverage' -D CMAKE_CXX_FLAGS='-coverage' -D CMAKE_EXE_LINKER_FLAGS=-coverage"
    cmake_build_type=Debug
    export PATH=$NGC/private/bin:${PATH}
fi


wonton_install_dir=${ngc_tpl_dir}/wonton/${wonton_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}
wonton_flags="-D WONTON_ROOT:PATH=$wonton_install_dir"

tangram_install_dir=${ngc_tpl_dir}/tangram/${tangram_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}
tangram_flags="-D PORTAGE_ENABLE_TANGRAM=True -D TANGRAM_ROOT:PATH=$tangram_install_dir"

flecsi_flags="-D PORTAGE_ENABLE_FleCSI:BOOL=False"  # Not building with FleCSI for HPC builds

if [[ $build_type != "serial" ]]; then
    jali_flags="-D PORTAGE_ENABLE_Jali:BOOL=True"  # Jali found through Wonton
fi


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
  $thrust_flags \
  $jali_flags \
  $flecsi_flags \
  $cov_flags \
  ..
make -j8
ctest -j36 --output-on-failure  && true #keep going if tests fail so that we get coverage report 
status=$?
if [[ $build_type == "coverage" ]]; then                     
    echo 'building coverage reports'
    export PYTHONPATH=/usr/projects/ngc/private/gcovr/var/lib/perceus/vnfs/asc-fe/rootfs/usr/lib/python2.7/site-packages
    /usr/projects/ngc/private/gcovr/var/lib/perceus/vnfs/asc-fe/rootfs/usr/bin/gcovr -f "$(readlink -f ..)"  -e '.*googletest' -e '.*exprtk.hpp' -e '.*json.h' -e '.*CMakeFiles' -x >coverage.xml
fi
exit $status #return the status of the ctest build so that jenkins knows whether tests past or fail
