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
WORKDIR /home/portage/portage
RUN mkdir build
WORKDIR build
RUN cmake -D CMAKE_BUILD_TYPE=Debug -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True -D ENABLE_MPI=${MPI} \
    -D ENABLE_MPI_CXX_BINDINGS=True ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} \
    -D ENABLE_DOXYGEN=True ..
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
