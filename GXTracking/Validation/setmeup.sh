
#.. TOPDIR is the top git area (where .git directory is)
export TOPDIR=/home/lima/work/g4hpcbenchmarks
export SRCDIR=${TOPDIR}/GXTracking/Validation

#.. BUILDDIR is where cmake and make commands are run
export BUILDDIR=/scratch/lima/valid/build

#.. INSTALLDIR is where "make install" will install libraries, binaries, cmake scripts, etc.
export INSTALLDIR=/scratch/lima/valid/install

#.. Provide locations for Geant4 and Root
export G4DIR=/scratch/lima/geant/geant4.9.6.p02
export ROOTSYS=/home/lima/work/root/v5.34.10

export PATH=${ROOTSYS}/bin:${INSTALLDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib

#.. setup cmake
source /products/setup
setup cmake

#.. where to launch validaiton jobs from
export RUNDIR=/home/lima/work/cmake/run

#.. useful aliases
alias gomake='cd ${BUILDDIR} && /usr/bin/make -j4 install; cd $RUNDIR'
alias goclean='cd ${BUILDDIR} && /usr/bin/make clean; cd $RUNDIR'

#.. Geant4 rebuilds
alias g4makeNormal='cd ${G4DIR}/build-normal && make -j20; cd ${RUNDIR}'
alias g4makeDebug='cd ${G4DIR}/build-debug && make -j20; cd ${RUNDIR}'

#.. for Geant4 jobs
export G4DATA=/home/g4p/pbs/download/g4data
export G4LEDATA=${G4DATA}/G4EMLOW6.32
export G4LEVELGAMMADATA=${G4DATA}/PhotonEvaporation2.3
export G4NEUTRONHPDATA=${G4DATA}/G4NDL4.3
export G4RADIOACTIVEDATA=${G4DATA}/RadioactiveDecay3.6
export G4ABLADATA=${G4DATA}/G4ABLA3.0
export G4REALSURFACEDATA=${G4DATA}/RealSurface1.0
export G4NEUTRONXSDATA=${G4DATA}/G4NEUTRONXS1.2
export G4PIIDATA=${G4DATA}/G4PII1.3
export G4SAIDXSDATA=${G4DATA}/G4SAIDDATA1.1

alias g4find='find /scratch/lima/geant/geant4.9.6.p01-orig/source'

#.. setup flag
export GXVALID=OK

cd ${RUNDIR}
