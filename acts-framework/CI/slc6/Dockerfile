FROM cern/slc6-base:latest

LABEL maintainer="Moritz Kiehn <msmk@cern.ch>"
LABEL version="1.0"
LABEL description="Scientific Linux CERN 6 with minimal ACTS dependencies"

# second set of packages is needed by Geant4 and ROOT from cvmfs
RUN yum -y update \
    && yum -y install \
       autoconf \
       automake \
       binutils \
       cmake \
       doxygen \
       gcc \
       gcc-c++ \
       gcc-fortran \
       git \
       gettext \
       libtool \
       make \
       patch \
       pkgconfig \
       which \
    && yum -y install \
       freeglut-devel \
       openmotif-devel \
       libXpm-devel \
    && yum -y clean all
