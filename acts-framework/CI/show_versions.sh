#!/bin/sh
#
# show available versions of external libraries

BOOST_INCLUDE_DIR=`cmake --find-package -DNAME=Boost -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=COMPILE | sed -e 's/-I//' | awk '{print $1}'`

BOOST_VERSION=`grep "#define BOOST_LIB_VERSION" ${BOOST_INCLUDE_DIR}/boost/version.hpp 2> /dev/null | awk '{print $3}' | sed -e 's/\"//g' -e 's/_/./g'`
CLANG_VERSION=`clang++ --version | head -n 1 | cut -d ' ' -f 3`
CLHEP_VERSION=`clhep-config --version 2> /dev/null | awk '{print $2}'`
CMAKE_VERSION=`cmake --version  2> /dev/null | head -n 1 | awk '{print $3}'`
DD4HEP_VERSION=`grep 'set ( DD4hep_VERSION' ${DD4hep_DIR}/DD4hepConfig.cmake 2> /dev/null | awk '{print $4}' | sed -e 's/"//g'`
DOXYGEN_VERSION=`doxygen --version  2> /dev/null`
EIGEN_WORLD_VERSION=`grep "#define EIGEN_WORLD_VERSION" ${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h 2> /dev/null | awk '{print $3}'`
EIGEN_MAJOR_VERSION=`grep "#define EIGEN_MAJOR_VERSION" ${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h 2> /dev/null | awk '{print $3}'`
EIGEN_MINOR_VERSION=`grep "#define EIGEN_MINOR_VERSION" ${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h 2> /dev/null | awk '{print $3}'`
EIGEN_VERSION="${EIGEN_WORLD_VERSION}.${EIGEN_MAJOR_VERSION}.${EIGEN_MINOR_VERSION}"
GCC_VERSION=`g++ -dumpversion 2> /dev/null`
GEANT4_VERSION=`geant4-config --version 2> /dev/null`
GIT_VERSION=`git --version 2> /dev/null | awk '{print $3}'`
MAKE_VERSION=`make --version  2> /dev/null | head -n 1 | awk '{print $3}'`
PYTHON_VERSION=`python --version  2>&1 | awk '{print $2}'`
ROOT_VERSION=`root-config --version 2> /dev/null`

BOOST_VERSION=${BOOST_VERSION:-"not found"}
CLANG_VERSION=${CLANG_VERSION:-"not found"}
CLHEP_VERSION=${CLHEP_VERSION:-"not found"}
CMAKE_VERSION=${CMAKE_VERSION:-"not found"}
DD4HEP_VERSION=${DD4HEP_VERSION:-"not found"}
DOXYGEN_VERSION=${DOXYGEN_VERSION:-"not found"}
EIGEN_VERSION=${EIGEN_VERSION:-"not found"}
GCC_VERSION=${GCC_VERSION:-"not found"}
GEANT4_VERSION=${GEANT4_VERSION:-"not found"}
GIT_VERSION=${GIT_VERSION:-"not found"}
MAKE_VERSION=${MAKE_VERSION:-"not found"}
PYTHON_VERSION=${PYTHON_VERSION:-"not found"}
ROOT_VERSION=${ROOT_VERSION:-"not found"}

echo ' *   __ *   __  ____* __ '
echo '    / \\   / _||_  _|/ _|'
echo '*  / / \\  ||    ||  \_\*'
echo '  / /===\\ ||_ * ||   _\\'
echo ' / /  *  \\\__|  || *|_ /'
echo
echo "****************************"
echo "* program      version     *"
echo "****************************"
echo "boost          $BOOST_VERSION"
echo "clang++        $CLANG_VERSION"
echo "clhep          $CLHEP_VERSION"
echo "cmake          $CMAKE_VERSION"
echo "dd4hep         $DD4HEP_VERSION"
echo "doxygen        $DOXYGEN_VERSION"
echo "eigen          $EIGEN_VERSION"
echo "g++            $GCC_VERSION"
echo "geant4         $GEANT4_VERSION"
echo "git            $GIT_VERSION"
echo "make           $MAKE_VERSION"
echo "python         $PYTHON_VERSION"
echo "root           $ROOT_VERSION"
echo "****************************"

