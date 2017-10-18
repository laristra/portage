FROM ubuntu

#to run this script: docker build -t docker-portage .
#to start the container: docker run -h simpleportage -it -v /Users:/Users --rm docker-portage
#to run as root in it: docker ls; docker exec -itu root <id>  bash

# use and arg to build the directory
ARG	user=portage

# set the proxies so we can get outside of LANL
ENV http_proxy="http://proxyout.lanl.gov:8080" \
	https_proxy="https://proxyout.lanl.gov:8080"

# build steps: get packages, then create a new group/user
RUN 	apt-get update && apt-get install -y \
		ccache texlive-latex-base texlive-fonts-recommended texlive-latex-recommended texlive-font-utils \
		cmake \
		g++ \
		libboost-dev \
		liblapack-dev liblapacke-dev \
		libopenmpi-dev openmpi-bin \
		python-pip git \
		wget curl lcov doxygen && \
 	groupadd -r ${user} && \
 	useradd -r -m -g ${user} ${user}

#RUN apt-get install -y openjdk-8-jdk unzip
#RUN wget https://sonarsource.bintray.com/Distribution/sonar-scanner-cli/sonar-scanner-2.8.zip https://sonarqube.com/static/cpp/build-wrapper-linux-x86.zip
#RUN unzip sonar-scanner-2.8.zip -d /sonarqube/
#RUN unzip build-wrapper-linux-x86.zip -d /sonarqube/
#ENV PATH=${PATH}${PATH:+:}/sonarqube/build-wrapper-linux-x86:/sonarqube/sonar-scanner-2.8/bin
#RUN sonar-scanner -h

USER 	${user}
WORKDIR /home/${user}
#COPY 	.ssh .ssh

# install coverage tools as user, and get portage
ENV 	PATH=/usr/lib/ccache:${PATH}
RUN 	git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/~shevitz/portage.git
	
	#pip install --user codecov coverxygen && \
 	 
	#&& \
 	##git clone --recursive https://github.com/laristra/portage && \
	#mkdir -p portage/build && \
	#cd portage/build && \
	#cmake -DENABLE_APP_TESTS=True .. && \
	#make

# leave the user in ~/portage/build
#WORKDIR ~/portage/build


#USER root
#RUN cd /home/portage/portage/build
#RUN make install
#USER portage


