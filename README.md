# Installation instructions

Below we list copy & paste instructions for several machines. You can easily
adapt them for other machines.

## Darwin

Execute the following from the portage root directory:

    module load openmpi/1.10.0-intel_15.0.3 cmake
    JALI_INSTALL_PREFIX=/path/to/inst-jali
    mkdir build
    cd build
    cmake \
        -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpiCC` \
        -D CMAKE_BUILD_TYPE=Debug \
        -D ENABLE_UNIT_TESTS=True \
        -D ENABLE_MPI=True \
        -D ENABLE_MPI_CXX_BINDINGS=True \
        -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
        ..
    make -j16
    ctest -j16
