FROM laristra/portage-buildenv:latest

ARG MPI

COPY portage /home/portage/portage
USER root
RUN chown -R portage:portage /home/portage/portage 
USER portage
WORKDIR /home/portage/portage
RUN mkdir build
WORKDIR build
RUN cmake -D CMAKE_BUILD_TYPE=Debug -D ENABLE_UNIT_TESTS=True -D ENABLE_APP_TESTS=True -D ENABLE_MPI=${MPI} -D ENABLE_MPI_CXX_BINDINGS=True ..
RUN make
RUN make test
