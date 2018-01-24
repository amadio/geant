// This file is use to gather all the source file
// that ought to be compiled when building in the
// single compilation unit nvcc mode.

#include "../src/GeantTaskData.cxx"
#include "../src/GeantTrack.cxx"
#include "../src/GeantTrackVec.cxx"
#include "../src/GeantTrackGeo.cxx"
#include "../src/GeantEvent.cxx"
#include "../src/ScalarNavInterfaceVGM.cxx"
#include "../src/VectorNavInterface.cxx"
#include "../src/GeantPropagator.cxx"

#include "../src/PreStepStage.cxx"
#include "../src/SimulationStage.cxx"
#include "../src/SteppingActionsStage.cxx"
#include "../src/GeomQueryStage.cxx"
#include "../src/ContinuousProcStage.cxx"
//#include "../src/LinearPropagationStage.cxx"
//#include "../src/FieldPropagationStage.cxx"
#include "../src/XSecSamplingStage.cxx"
#include "../src/DiscreteProcStage.cxx"
#include "../src/PropagationStage.cxx"

#include "../src/PreStepHandler.cxx"
#include "../src/SteppingActionsHandler.cxx"
#include "../src/GeomQueryHandler.cxx"
#include "../src/ContinuousProcHandler.cxx"
#include "../src/LinearPropagationHandler.cxx"
#include "../src/FieldPropagationHandler.cxx"
#include "../src/XSecSamplingHandler.cxx"
#include "../src/DiscreteProcHandler.cxx"

#include "../src/NumaUtils.cxx"
#include "../src/Handler.cxx"
#include "../src/Basket.cxx"

//#include "../src/WorkloadManager.cxx"
//#include "../src/GeantRunManager.cxx"


#include "GeantCudaUtils.cu"
#include "GeantTrack.cu"
#include "PropagateTracks.cu"

// #include "test_compile.cu"
