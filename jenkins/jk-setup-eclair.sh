#!/bin/bash -x

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


THIS=$(dirname ${BASH_SOURCE[0]})

# first arguments is the source directory
if [ $# -ge 6 ]; then
  LABEL=$1 ; shift
  COMPILER=$1 ; shift
  BUILDTYPE=$1 ; shift
  EXTERNALS=$1 ; shift
  WORKSPACE=$1 ; shift
  BACKEND=$1 ; shift
else
  echo "$0: expecting 5 arguments [LABEL]  [COMPILER] [BUILDTYPE] [EXTERNALS] [WORKSPACE] [BACKEND]"
  return
fi

if [ $LABEL == slc6 ] || [ $LABEL == cc7 ]
then
  export PATH=/afs/cern.ch/sw/lcg/contrib/CMake/3.0.0/Linux-i386/bin:${PATH}
else
  export EXTERNALDIR=$HOME/ROOT-externals/
fi

if [[ $COMPILER == *gcc* ]]
then
  gcc47version=4.7
  gcc48version=4.8
  gcc49version=4.9
  COMPILERversion=${COMPILER}version

  ARCH=$(uname -m)
  . /afs/cern.ch/sw/lcg/contrib/gcc/${!COMPILERversion}/${ARCH}-${LABEL}/setup.sh
  #export FC=gfortran
  #export CXX=`which g++`
  #export CC=`which gcc`

  export CMAKE_SOURCE_DIR=$WORKSPACE/geant
  export CMAKE_BINARY_DIR=$WORKSPACE/geant/builds
  export CMAKE_BUILD_TYPE=$BUILDTYPE

#  export CTEST_BUILD_OPTIONS="-DUSE_VECGEOM_NAVIGATOR=OFF '-DCMAKE_CXX_FLAGS=-O2 -std=c++11' -DUSE_ROOT=ON -DCTEST=ON "
  export CTEST_BUILD_OPTIONS=" '-DCMAKE_CXX_FLAGS=-O2 -std=c++11' -DUSE_ROOT=ON -DCTEST=ON "
  export CMAKE_INSTALL_PREFIX=$WORKSPACE/geant/installation
  export BACKEND=$BACKEND
#  export BACKEND=VecGeomNavigator
  export LD_LIBRARY_PATH=$WORKSPACE/geant/lib:$LD_LIBRARY_PATH

fi

echo ${THIS}/setup.py -o ${LABEL} -c ${COMPILER} -b ${BUILDTYPE} -v ${EXTERNALS} -w ${WORKSPACE}
eval `${THIS}/setup.py -o ${LABEL} -c ${COMPILER} -b ${BUILDTYPE} -v ${EXTERNALS} -w ${WORKSPACE}`

############## Eclair setings ###################

#set -e

export PROJECT_ROOT=${WORKSPACE}
export ANALYSIS_DIR=${WORKSPACE}/geant/jenkins
export OUTPUT_DIR=${WORKSPACE}/GeantV
mkdir -p ${OUTPUT_DIR}
export PB_OUTPUT=${OUTPUT_DIR}/REPORT.@FRAME@.pb
export ECLAIR_DIAGNOSTICS_OUTPUT=${OUTPUT_DIR}/DIAGNOSTICS.txt
export ECL_CONFIG_FILE=${ANALYSIS_DIR}/GeantV.ecl
rm -f ${ECLAIR_DIAGNOSTICS_OUTPUT}
rm -f ${OUTPUT_DIR}/REPORT.*.pb

export CC=/usr/bin/cc
export CXX=/usr/bin/c++
export AS=/usr/bin/as
export LD=/usr/bin/ld
export AR=/usr/bin/ar

eclair_env +begin +end -eval-file=${ECL_CONFIG_FILE} -- 'cmake ../ $CTEST_BUILD_OPTION && make -j 24' || true

cd ${OUTPUT_DIR}
rm -f REPORTS.db
eclair_report -create-db=REPORTS.db *.pb -load
rm -rf eclair_output
eclair_report -db=REPORTS.db -output=eclair_output/@TAG@.etr -reports
