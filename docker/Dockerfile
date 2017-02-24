FROM laristra/portage-buildenv:latest

COPY portage /home/portage/portage
USER root
RUN chown -R portage:portage /home/portage/portage 
USER portage
WORKDIR /home/portage/portage
RUN mkdir build
WORKDIR build
RUN cmake -D CMAKE_C_COMPILER=`which mpicc`  -D CMAKE_CXX_COMPILER=`which mpiCC` -D CMAKE_BUILD_TYPE=Debug -D ENABLE_UNIT_TESTS=True -D ENABLE_APP_TESTS=True -D ENABLE_MPI=True -D ENABLE_MPI_CXX_BINDINGS=True ..
RUN make
RUN make test
