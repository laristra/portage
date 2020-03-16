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
ENV TANGRAM_INSTALL_DIR /home/portage/tangram/install 
RUN cmake -D CMAKE_BUILD_TYPE=Debug \ 
    -D CMAKE_INSTALL_PREFIX=$TANGRAM_INSTALL_DIR \ 
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True -D ENABLE_MPI=${MPI} \
    -D ENABLE_MPI_CXX_BINDINGS=True .. 
RUN make && make test 
RUN make install
RUN ls -l $TANGRAM_INSTALL_DIR

#Build XMOF2D
WORKDIR ../..
RUN git clone https://github.com/laristra/XMOF2D.git
WORKDIR XMOF2D
RUN mkdir build
WORKDIR build
ENV XMOF2D_INSTALL_DIR /home/portage/XMOF2D/install
RUN cmake -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_CXX_FLAGS:STRING="-std=c++11" \ 
    -D INSTALL_DIR:FILEPATH=$XMOF2D_INSTALL_DIR ..  
RUN make && make install 
RUN ls -l $XMOF2D_INSTALL_DIR 

# Build Portage
RUN echo $JALI_DIR
RUN ls -l $JALI_DIR/lib
WORKDIR ../../portage
RUN mkdir build
WORKDIR build
RUN cmake -D CMAKE_BUILD_TYPE=Debug -D CMAKE_PREFIX_PATH=$TANGRAM_DIR \ 
    -D ENABLE_UNIT_TESTS=True -D ENABLE_APP_TESTS=True \ 
    -D ENABLE_MPI=${MPI} -D ENABLE_MPI_CXX_BINDINGS=True \
    ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} -D ENABLE_DOXYGEN=True \
    -D TANGRAM_DIR=$TANGRAM_INSTALL_DIR \
    -D XMOF2D_DIR=$XMOF2D_INSTALL_DIR/share/cmake \
    -D Jali_DIR:FILEPATH=$JALI_DIR/lib/cmake ..
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
