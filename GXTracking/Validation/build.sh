#!/bin/sh
#
# File:    GXTracking/Validation/build.sh
# Purpose:
#
#  A script to build the validation libraries and binaries, including
#  the choice of verbosity and randomness.
#
#  Running this script is also the easiest way to change modes, since
#  it starts by completely removing previous setups.  The user should
#  select the running mode by commenting out only the appropriate
#  cmake command below.
#

if [ ${GXVALID}"x" == "x" ]; then
  source ./setmeup.sh
fi

rm -rf ${BUILDDIR} 
mkdir -p ${BUILDDIR}
cd ${BUILDDIR}

# if cmake cleanup is needed
rm -rf CMakeFiles CMakeCache.txt Makefile cmake_install.cmake

#.. build-debug + nonRandom (GPUNONRANDOM) + printouts (GPUDEBUG) + plots (GPUPLOTS)
cmake -DGeant4_DIR=${G4DIR}/build-debug \
  -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DCMAKE_BUILD_TYPE=Debug \
  -DGPUDEBUG=ON -DGPUNONRANDOM=ON -DGPUPLOTS=ON -DCMAKE_BUILD_TYPE=Debug ${SRCDIR}

#.. build-debug + printouts (GPUDEBUG) + plots (GPUPLOTS)
#cmake -DGeant4_DIR=${G4DIR}/build-debug \
#  -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DCMAKE_BUILD_TYPE=Debug \
#  -DGPUDEBUG=ON -DGPUNONRANDOM=OFF -DGPUPLOTS=ON -DCMAKE_BUILD_TYPE=Debug ${SRCDIR}

#.. build-normal + plots (GPUPLOTS) - no printouts
#cmake -DGeant4_DIR=${G4DIR}/build-normal \
#  -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DCMAKE_BUILD_TYPE=Debug \
#  -DGPUDEBUG=OFF -DGPUNONRANDOM=OFF -DGPUPLOTS=ON  ${SRCDIR}

make -j4 install
cd ${RUNDIR}
