#!/bin/bash

#-----------------------------------------------------------------------
# Default options
#-----------------------------------------------------------------------

unset G4P_USE_CLHEP ; G4P_USE_CLHEP=0
unset G4P_CXX       ; G4P_CXX=/usr/bin/c++

#-----------------------------------------------------------------------
# Create the directory structure
#-----------------------------------------------------------------------
umask 0002

WORK_DIR="/scratch/lima/geant"
BUILD_DIR="${WORK_DIR}/geant4.9.6.p01/build"
INSTALL_DIR="${WORK_DIR}/geant4.9.6.p01"

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

#-----------------------------------------------------------------------
# configure with cmake
#-----------------------------------------------------------------------
XERCESC_DIR=/products/xerces_c/v3_1_1/Linux64bit+2.6-2.12-gcc46-prof
export XERCESC_DIR

cmake -DCMAKE_CXX_COMPILER=${G4P_CXX} \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG" \
      -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O2 -g -fno-omit-frame-pointer" \
      -DGEANT4_USE_SYSTEM_CLHEP=0 \
      -DGEANT4_INSTALL_DATA=0 \
      -DGEANT4_USE_GDML=ON \
      -DGEANT4_USE_OPENGL_X11=ON \
      -DXERCESC_ROOT_DIR=${XERCESC_DIR} \
       .. 
#      ${INSTALL_DIR}/source ${INSTALL_DIR}
#      -DCMAKE_BUILD_TYPE=RelWithDebInfo \

#      -DNO_BACKPORT_CUDA_REMOVALS=1 \

#-----------------------------------------------------------------------
# build and install
#-----------------------------------------------------------------------
make -j20
make install

