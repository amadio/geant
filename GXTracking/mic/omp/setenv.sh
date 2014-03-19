#!/bin/bash
#
# compilers and compilation options
# -DCMAKE_TOOLCHAIN_FILE=<path_to_this_file> or set env inside CMAkeLists.txt
#
export CC=icc
export CXX=icpc
export FC=ifort
#export CFLAGS="-openmp -O3 -opt-threads-per-core=2"
export CFLAGS="-openmp -O3 -DUSE_MIC"
export CXXFLAGS=$CFLAGS
export FFLAGS=$CFLAGS
export MPI_C=mpiicc
export MPI_CXX=mpiicpc
#
# setup env for the intel compiler (PATH, LD_LIBRARY_PATH, etc.) 
#
if [ -e /opt/intel/bin/iccvars.sh ]; then
   .    /opt/intel/bin/iccvars.sh intel64 
fi
