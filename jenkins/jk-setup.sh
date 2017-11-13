#!/bin/bash -x

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

THIS=$(dirname ${BASH_SOURCE[0]})

# first arguments is the source directory
if [ $# -ge 4 ]; then
  LABEL=$1 ; shift
  COMPILER=$1 ; shift
  BUILDTYPE=$1 ; shift
  EXTERNALS=$1 ; shift
#  WORKSPACE=$1 ; shift
#  TYPE=$1 ; shift
#  BACKEND=$1 ; shift
else
  echo "$0: expecting 4 arguments [LABEL] [COMPILER] [BUILDTYPE] [EXTERNALS]"
  return
fi

PLATFORM=`$THIS/getPlatform.py`
COMPATIBLE=`$THIS/getCompatiblePlatform.py $PLATFORM`
ARCH=$(uname -m)

export BUILDTYPE
export COMPILER

# Set up the externals against devgeantv in CVMFS
if [ -a /cvmfs/sft.cern.ch/lcg/views/devgeantv/latest/$PLATFORM ]; then
  source /cvmfs/sft.cern.ch/lcg/views/devgeantv/latest/$PLATFORM/setup.sh
elif [ -a /cvmfs/sft.cern.ch/lcg/views/devgeantv/latest/$COMPATIBLE ]; then
  source /cvmfs/sft.cern.ch/lcg/views/devgeantv/latest/$COMPATIBLE/setup.sh
else
  echo "No externals for $PLATFORM in $EXTERNALDIR/$EXTERNALS"
fi

if [ $LABEL == slc6 ] || [ $LABEL == gvslc6 ] || [ $LABEL == cc7 ] || [ $LABEL == cuda7 ] || [ $LABEL == slc6-physical ] || [ $LABEL == lcgapp-SLC6_64b ] || [  $LABEL == continuous-sl6 ] || [  $LABEL == continuous-cuda7 ] || [ $LABEL == continuous-xeonphi ]
then
  export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.7.0/Linux-$ARCH/bin/:${PATH}
  kinit sftnight@CERN.CH -5 -V -k -t /ec/conf/sftnight.keytab
elif [ $LABEL == xeonphi ]
then
  export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.7.0/Linux-$ARCH/bin:${PATH}
  kinit sftnight@CERN.CH -5 -V -k -t /data/sftnight/ec/conf/sftnight.keytab
fi

if [[ $COMPILER == *gcc* ]]; then
  gcc47version=4.7
  gcc48version=4.8
  gcc49version=4.9
  COMPILERversion=${COMPILER}version
  ARCH=$(uname -m)
  if [ $LABEL == cuda7 ] || [ $LABEL == gvslc6 ] || [ $LABEL == slc6-physical ] ||  [ $LABEL == lcgapp-SLC6_64b ] || [  $LABEL == continuous-sl6 ] || [  $LABEL == continuous-cuda7 ]; then
    . /cvmfs/sft.cern.ch/lcg/contrib/gcc/${!COMPILERversion}/${ARCH}-slc6/setup.sh
  else
    . /cvmfs/sft.cern.ch/lcg/contrib/gcc/${!COMPILERversion}/${ARCH}-${LABEL}/setup.sh
  fi
  export FC=gfortran
  export CXX=`which g++`
  export CC=`which gcc`
elif [[ $COMPILER == *native* && $PLATFORM == *mac* ]]; then
  export LD_LIBRARY_PATH=/usr/local/gfortran/lib
  export PATH=/usr/bin:/usr/local/bin:/opt/X11/bin
  export CC=`which clang`
  export CXX=`which clang++`
  export FC=`which gfortran`
elif [[ $PLATFORM == *native* ]]; then
  export CC=`which gcc`
  export CXX=`which g++`
  export FC=`which gfortran`
elif [[ $COMPILER == *icc* ]]; then
  iccyear=2013
  icc14year=2013
  icc15year=2015
  icc16year=2016
  COMPILERyear=${COMPILER}year
  iccgcc=4.9
  icc14gcc=4.9
  icc15gcc=4.9
  icc16gcc=4.9
  GCCversion=${COMPILER}gcc
  ARCH=$(uname -m)
  . /afs/cern.ch/sw/lcg/contrib/gcc/${!GCCversion}/${ARCH}-slc6/setup.sh
  . /afs/cern.ch/sw/IntelSoftware/linux/setup.sh
  . /afs/cern.ch/sw/IntelSoftware/linux/${ARCH}/xe${!COMPILERyear}/bin/ifortvars.sh intel64
  . /afs/cern.ch/sw/IntelSoftware/linux/${ARCH}/xe${!COMPILERyear}/bin/iccvars.sh intel64
  export CC=icc
  export CXX=icc
  export FC=ifort
elif [[ $COMPILER == *icc17* ]]; then
    icc17gcc=6.2
    . /cvmfs/sft.cern.ch/lcg/contrib/gcc/${!GCCversion}/${ARCH}-${LABEL_COMPILER}/setup.sh
    . /cvmfs/projects.cern.ch/intelsw/psxe/linux/all-setup.sh
    export CXX=`which icpc`
    export FC=`which gfortran`
    export LDFLAGS="-lirc -limf"
elif [[ $COMPILER == *clang* ]]; then
  clang34version=3.4
  clang35version=3.5
  clang36version=3.6
  clang37version=3.7
  clang38version=3.8
  COMPILERversion=${COMPILER}version
  clang34gcc=48
  clang35gcc=49
  clang37gcc=49
  clang38gcc=49
  GCCversion=${COMPILER}gcc
 . /cvmfs/sft.cern.ch/lcg/contrib/llvm/${!COMPILERversion}/${ARCH}-${LABEL_COMPILER}/setup.sh
  export CC=`which clang`
  export CXX=`which clang++`
  export FC=`which gfortran`
fi

export CMAKE_SOURCE_DIR=$WORKSPACE/geant
export CMAKE_BINARY_DIR=$WORKSPACE/geant/builds
export CMAKE_BUILD_TYPE=$BUILDTYPE

export CTEST_BUILD_OPTIONS=" -DCMAKE_CXX_STANDARD=14 -DUSE_ROOT=ON -DCTEST=ON ${ExtraCMakeOptions}"
export CMAKE_INSTALL_PREFIX=$WORKSPACE/geant/installation
export LD_LIBRARY_PATH=$WORKSPACE/lib:$LD_LIBRARY_PATH