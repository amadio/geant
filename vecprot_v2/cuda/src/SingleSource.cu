// This file is use to gather all the source file
// that ought to be compiled when building in the
// single compilation unit nvcc mode.

#include "../src/GeantTaskData.cxx"
#include "../src/GeantTrack.cxx"
#include "../src/GeantTrackVec.cxx"
#include "../src/GeantTrackGeo.cxx"
#include "../src/GeantEvent.cxx"
#include "../src/ScalarNavInterfaceVGM.cxx"
#include "../src/GeantPropagator.cxx"

//#include "../src/WorkloadManager.cxx"
//#include "../src/GeantRunManager.cxx"


#include "GeantCudaUtils.cu"
#include "GeantTrack.cu"
#include "PropagateTracks.cu"

// #include "test_compile.cu"
