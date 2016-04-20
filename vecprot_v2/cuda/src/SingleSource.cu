// This file is use to gather all the source file
// that ought to be compiled when building in the
// single compilation unit nvcc mode.

#include "../src/GeantTrack.cxx"
#include "../src/GeantTaskData.cxx"
#include "../src/ScalarNavInterfaceVGM.cxx"
#include "GeantTrack.cu"
#include "PropagateTracks.cu"
#include "GeantCudaUtils.cu"

// #include "test_compile.cu"
