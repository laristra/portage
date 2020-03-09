FROM laristra/portage-buildenv:ubuntu

ARG MPI
ARG COVERAGE

#for coverage
ARG CI
ARG TRAVIS
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT
ARG TRAVIS_OS_NAME

# for docs
ARG DOCS

COPY portage /home/portage/portage

USER root
RUN chown -R portage:portage /home/portage/portage
USER portage

# Build tangram
WORKDIR /home/portage
RUN git clone https://github.com/laristra/tangram.git 
WORKDIR tangram
RUN git submodule update --init --recursive
RUN git checkout v0.9.7 
RUN mkdir build
WORKDIR build
ENV TANGRAM_DIR /home/portage/tangram/install 
RUN cmake -D CMAKE_BUILD_TYPE=Debug -D CMAKE_INSTALL_PREFIX=$TANGRAM_DIR \ 
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True -D ENABLE_MPI=${MPI} \
    -D ENABLE_MPI_CXX_BINDINGS=True .. 
RUN make && make test 
RUN make install
RUN echo $TANGRAM_DIR

# Build MSTK and Jali 
WORKDIR ../..
RUN wget https://github.com/MeshToolkit/mstk/archive/3.3.2.tar.gz &&
     tar xvzf 3.3.2.tar.gz 
WORKDIR MSTK-3.3.2 
RUN mkdir build
WORKDIR build
ENV MSTK_DIR /home/portage/MSTK-3.3.2/install
RUN  cmake -D CMAKE_BUILD_TYPE=Release \
     -D CMAKE_INSTALL_PREFIX=$MSTK_DIR \
     -D ENABLE_PARALLE=$PARALLEL \
     -D ENABLE_METIS=yes \
     -D ENABLE_ExodusII=yes \
     -D ENABLE_Zoltan=yes \
     -D ENABLE_Tests=yes .. 
RUN make && make test 
RUN make install 
RUN echo $MSTK_DIR

WORKDIR ../..
RUN git clone https://github.com/lanl/jali.git
WORKDIR jali
RUN git checkout 1.0.5 
RUN mkdir build
WORKDIR build
ENV JALI_DIR /home/portage/jali/install
RUN cmake -D CMAKE_PREFIX_PATH=$MSTK_DIR \
     -D CMAKE_INSTALL_PREFIX=$JALI_DIR \
     -D CMAKE_BUILD_TYPE=Release \
     -D ENABLE_MSTK=yes \
     -D ENABLE_Tests=yes .. \
RUN make && make test 
RUN make install 
RUN echo $JALI_DIR 

# Build Portage
WORKDIR ../../portage
RUN mkdir build
WORKDIR build
RUN cmake -D CMAKE_BUILD_TYPE=Debug -D CMAKE_PREFIX_PATH=$TANGRAM_DIR \ 
    -D ENABLE_UNIT_TESTS=True -D ENABLE_APP_TESTS=True \ 
    -D ENABLE_MPI=${MPI} -D ENABLE_MPI_CXX_BINDINGS=True \
    ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} -D ENABLE_DOXYGEN=True \
    -D TANGRAM_DIR=$TANGRAM_DIR \
    -D Jali_DIR:FILEPATH=$JALI_DIR/lib ..
RUN make VERBOSE=1 -j2
RUN make test
RUN cd .. && if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov; \
  fi; \
fi && cd -
RUN if [ ${DOCS} ] || [ -n "${TRAVIS_TAG}" ]; then \
  make doxygen; \
fi
