FROM laristra/portage-buildenv:latest

ARG MPI
ARG COVERAGE
ARG SONARQUBE
ARG SONARQUBE_TOKEN
ARG SONARQUBE_GITHUB_TOKEN

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
RUN cmake -D CMAKE_BUILD_TYPE=Debug -D ENABLE_UNIT_TESTS=True -D ENABLE_APP_TESTS=True -D ENABLE_MPI=${MPI} -D ENABLE_MPI_CXX_BINDINGS=True ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} ..
RUN ${SONARQUBE:+build-wrapper-linux-x86-64 --out-dir bw-output} make VERBOSE=1 -j2
RUN make test
RUN cd .. && if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov; \
  fi; \
fi && cd -
RUN if [ ${DOCS} ]; then \
  make doxygen; \
fi
RUN cd .. && if [ ${SONARQUBE} ]; then \
  sonar-scanner \
    -Dsonar.projectKey=${TRAVIS_REPO_SLUG##*/} \
    -Dsonar.projectName=${TRAVIS_REPO_SLUG#*/} \
    -Dsonar.projectVersion=${TRAVIS_COMMIT} \
    -Dsonar.branch=/${TRAVIS_BRANCH} \
    -Dsonar.links.homepage=http://${TRAVIS_REPO_SLUG%%/*}.github.io/${TRAVIS_REPO_SLUG#*/} \
    -Dsonar.links.ci=https://travis-ci.org/${TRAVIS_REPO_SLUG} \
    -Dsonar.links.scm=https://github.com/${TRAVIS_REPO_SLUG} \
    -Dsonar.links.issue=https://github.com/${TRAVIS_REPO_SLUG}/issues \
    -Dsonar.sources=. \
    -Dsonar.exclusions=build/CMakeFiles/**,cinch/**,build/doc/doxygen/**,build/include/** \
    -Dsonar.cfamily.build-wrapper-output=build/bw-output \
    -Dsonar.cfamily.gcov.reportsPath=. \
    -Dsonar.host.url=https://sonarcloud.io \
    -Dsonar.organization=${TRAVIS_REPO_SLUG%%/*} \
    $([ ${TRAVIS_PULL_REQUEST} != false ] && echo \
       -Dsonar.analysis.mode=preview \
       -Dsonar.github.pullRequest=${TRAVIS_PULL_REQUEST} \
       -Dsonar.github.repository=${TRAVIS_REPO_SLUG} \
       -Dsonar.github.oauth=${SONARQUBE_GITHUB_TOKEN}) \
    -Dsonar.login=${SONARQUBE_TOKEN}; \
fi && cd -
