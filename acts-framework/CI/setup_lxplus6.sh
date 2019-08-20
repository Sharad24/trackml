#!/bin/sh -ex
#
# setup requirements on lxplus6

platform=x86_64-slc6-gcc62-opt
view=/cvmfs/sft.cern.ch/lcg/views/LCG_90a/${platform}

source ${view}/setup.sh
# additional variables that are not set automatically
export BOOST_ROOT="${view}"
export DD4hep_DIR="${view}"
export EIGEN_INCLUDE_DIR="${view}/include/eigen3"
export PYTHIA8_INCLUDE_DIR="${view}/include"
export PYTHIA8_LIBRARY_DIR="${view}/lib"
