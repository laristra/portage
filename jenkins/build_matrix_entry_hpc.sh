#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry_hpc.sh BUILD_TYPE <VER> <WONTON_VER>
#
# BUILD_TYPE - pr, nightly, install
#
# if VER is abset, the HEAD of the master branch will be taken. If
# WONTON_VER is absent, the HEAD of the master branch of wonton will
# be taken. If BUILD_TYPE is 'install' and VER is specified, it will
# install it to /install_prefix/tangram/$VER-blah-blah; if VER is not
# specified, it will install to /install_prefix/wonton/dev-blah-blah
#
# WORKSPACE   -  where the code is checked out
# CONFIG_TYPE -  base, debug, serial, readme, thrust, kokkos
# COMPILER    -  intel, gcc6, gcc7
# BRANCH_NAME -  master
#
# The exit code determines if the test succeeded or failed.


# Exit on error
set -e
# Echo each command
set -x

# set umask so installations will have group rwx permission
umask 007

BUILD_TYPE=$1
version=$2
if [[ $version == "" ]]; then
    version=dev
fi
tangram_version=$3
if [[ $tangram_version == "" ]]; then
    tangram_version=dev
fi
wonton_version=$3
if [[ $wonton_version == "" ]]; then
    wonton_version=dev
fi

echo "inside build_matrix on PLATFORM=$PLATFORM with BUILD_TYPE=$BUILD_TYPE $CONFIG_TYPE=$CONFIG_TYPE COMPILER=$COMPILER"

# special case for README builds
if [[ $BUILD_TYPE != "install" && $CONFIG_TYPE == "readme" ]]; then

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

export NGC=/usr/projects/ngc
ngc_include_dir=$NGC/private/include


# compiler-specific settings
if [[ $COMPILER =~ "intel" ]]; then

    compiler_version=18.0.5
    cxxmodule=intel/${compiler_version}
    compiler_suffix="-intel-${compiler_version}"

    openmpi_version=2.1.2
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"
    
elif [[ $COMPILER =~ "gcc" ]]; then

    openmpi_version=2.1.2
    if [[ $COMPILER == "gcc6" ]]; then
	compiler_version=6.4.0
    elif [[ $COMPILER == "gcc7" ]]; then
	compiler_version=7.4.0
    fi
    
    cxxmodule=gcc/${compiler_version}
    compiler_suffix="-gcc-${compiler_version}"

    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"

fi

# Jali
jali_flags="-D PORTAGE_ENABLE_Jali::BOOL=True"

# FleCSI
flecsi_flags="-D PORTAGE_ENABLE_FleCSI:BOOL=False"  # Not building with FleCSI for HPC builds

# THRUST
thrust_flags=
thrust_suffix=
if [[ $CONFIG_TYPE == "thrust" ]]; then
    thrust_flags="-D PORTAGE_ENABLE_THRUST=True"
    thrust_suffix="-thrust"
fi

# MPI or not
mpi_flags="-D PORTAGE_ENABLE_MPI=True"
if [[ $CONFIG_TYPE == "serial" ]]; then
    mpi_flags="-D PORTAGE_ENABLE_MPI=False"
    mpi_suffix=
    jali_flags=
    flecsi_flags=
fi

# Debug or Optimized build
cmake_build_type=Release
debug_suffix=
if [[ $CONFIG_TYPE == "debug" ]]; then
    cmake_build_type=Debug
    debug_suffix="-debug"
fi

# WONTON
wonton_install_dir=$NGC/private/wonton/${wonton_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}${kokkos_suffix}${debug_suffix}
wonton_flags="-D WONTON_ROOT:FILEPATH=$wonton_install_dir"

# TANGRAM
tangram_install_dir=$NGC/private/tangram/${tangram_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}${kokkos_suffix}${debug_suffix}
tangram_flags="-D PORTAGE_ENABLE_TANGRAM=True -D TANGRAM_ROOT:FILEPATH=$tangram_install_dir"

# Build up an install dir name
portage_install_dir=$NGC/private/portage/${version}${compiler_suffix}${mpi_suffix}${thrust_suffix}${kokkos_suffix}${debug_suffix}


# Coverage
cov_flags=""
if [[ $CONFIG_TYPE == "coverage" ]]; then
    cov_flags="-D CMAKE_C_FLAGS='-coverage' -D CMAKE_CXX_FLAGS='-coverage' -D CMAKE_EXE_LINKER_FLAGS=-coverage"
    cmake_build_type=Debug
    export PATH=$NGC/private/bin:${PATH}
fi


#Rely on default user environment to load modules; these scripts can be found in /etc/profile.d/ as of 8/25/20 the scripts that set up modules on snow are  /etc/profile.d/z00_lmod.sh; /etc/profile.d/00-modulepath.sh; /etc/profile.d/z01-modules.lanl.sh;
module load $cxxmodule
if [[ -n "$mpi_flags" ]]; then
    module load ${mpi_module}
fi
module load cmake/3.14.0 # 3.13 or higher is required

echo "JENKINS WORKSPACE = $WORKSPACE"
cd $WORKSPACE

rm -fr build
mkdir build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=$cmake_build_type \
  -D CMAKE_CXX_FLAGS="-Wall -Werror" \
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

make -j36
ctest -j36 --output-on-failure  && true #keep going if tests fail so that we get coverage report 
status=$?

if [[ $CONFIG_TYPE == "coverage" ]]; then                     
    echo 'building coverage reports'
    export PYTHONPATH=/usr/projects/ngc/private/gcovr/var/lib/perceus/vnfs/asc-fe/rootfs/usr/lib/python2.7/site-packages
    /usr/projects/ngc/private/gcovr/var/lib/perceus/vnfs/asc-fe/rootfs/usr/bin/gcovr -f "$(readlink -f ..)"  -e '.*googletest' -e '.*exprtk.hpp' -e '.*json.h' -e '.*CMakeFiles' -x >coverage.xml
fi

if [[ $BUILD_TYPE == "config" ]]; then
    make install
fi

exit $status #return the status of the ctest build so that jenkins knows whether tests past or fail
