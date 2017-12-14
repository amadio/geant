#include "TestProcess.h"

#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "globals.h"
#include "GeantTrackVec.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#else
#include "TGeoBranchArray.h"
#endif


TestProcess::TestProcess() : Geant::PhysicsProcessOld() {}


void TestProcess::Initialize() {
  std::cout << "TestProcess::Initialize : Start" << std::endl;
  std::cout << "TestProcess::Initialize : End" << std::endl;
}


void TestProcess::Eloss( Material_t */* mat */, int ntracks, GeantTrack_v &tracks, int &/* nout */, GeantTaskData */* td */ ) {

  //  std::cout << "TestProcess::Eloss : Start" << std::endl;

  for ( int i = 0; i < ntracks; i++ ) {
    tracks.fEdepV[i] = 10;
  }
  //  std::cout << "TestProcess::Eloss : --- End ---" << std::endl;
}


void TestProcess::ComputeIntLen( Material_t */* mat */, int ntracks, GeantTrack_v &tracks,
				    GeantTaskData */* td */ ) {
  //  std::cout << "TestProcess::ComputeIntLen : Start : ntracks=" << ntracks << std::endl;
  for ( int i = 0; i < ntracks; i++ ) {
    //    std::cout << " TestProcess::ComputeIntLen : proposedStepLengths[" << i << "]=" << 10
    //             << " ; position=("  << tracks.fXposV[i] << "," << tracks.fYposV[i] << "," << tracks.fZposV[i]
    //             << ")" << std::endl;  // Debug
    tracks.fPstepV[i] = 10;
    tracks.fEindexV[i] = 1000;  // Forced to be always treated as continuous & discrete.
  }
  //  std::cout << "TestProcess::ComputeIntLen : --- End ---" << std::endl;
}


void TestProcess::PostStepTypeOfIntrActSampling( Material_t */* mat */, int /* ntracks */, GeantTrack_v &/* tracks */,
                                                    GeantTaskData */* td */ ) {
  //  std::cout << "TestProcess::PostStepTypeOfIntrActSampling : Start & End" << std::endl;
}


void TestProcess::PostStepFinalStateSampling( Material_t */* mat */, int /* ntracks */, GeantTrack_v &/* tracks */, int &/* nout */,
                                                 GeantTaskData */* td */ ) {
  //  std::cout << "TestProcess::PostStepTypeOfIntrActSampling : Start" << std::endl;
  //  std::cout << "TestProcess::PostStepTypeOfIntrActSampling : --- End ---" << std::endl;
}
