FROM cern/cc7-base:latest

LABEL maintainer="Moritz Kiehn <msmk@cern.ch>"
LABEL version="1.0"
LABEL description="CERN CentOS 7 with minimal ACTS dependencies"

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
       motif-devel \
       ninja-build \
       patch \
       pkgconfig \
       which \
    && yum -y clean all
