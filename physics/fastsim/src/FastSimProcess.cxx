#include "FastSimProcess.h"

#include "Smearer.h"
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

FastSimProcess::FastSimProcess() : PhysicsProcessOld() , fSmearer( nullptr ) {}


FastSimProcess::~FastSimProcess() {
  delete fSmearer;
}


void FastSimProcess::Initialize(GeantTaskData *td) {
  //std::cout << "FastSimProcess::Initialize : Start" << std::endl;  // Debug
  if ( fSmearer ) delete fSmearer;
  fSmearer = new Smearer(td);
  //std::cout << "FastSimProcess::Initialize : --- End ---" << std::endl;  // Debug
}


void FastSimProcess::ComputeIntLen( Material_t * /* mat */, int ntracks, GeantTrack_v &tracks,
                                    double * /* lengths */, GeantTaskData *td ) {
  //std::cout << "FastSimProcess::ComputeIntLen : Start : ntracks=" << ntracks << std::endl;  // Debug
  std::vector< double > proposedStepLengths =
    fSmearer->StepLengthProposedByParameterisation( ntracks, tracks, *td );
  for ( int i = 0; i < ntracks; i++ ) {
    //std::cout << " FastSimProcess::ComputeIntLen : proposedStepLengths[" << i
    //          << "]=" << proposedStepLengths[i] << " ; position=("  << tracks.fXposV[i]
    //          << "," << tracks.fYposV[i] << "," << tracks.fZposV[i] << ")"
    //          << std::endl;  // Debug
    tracks.fPstepV[i] = proposedStepLengths[i];
    tracks.fEindexV[i] = 1000;  // Forced to be always treated as continuous & discrete.
  }
  //std::cout << "FastSimProcess::ComputeIntLen : --- End ---" << std::endl;  // Debug
}


void FastSimProcess::Eloss( Material_t * /* mat */, int ntracks, GeantTrack_v &tracks,
                            int & /* nout */, GeantTaskData *td ) {
  //std::cout << "FastSimProcess::Eloss : Start" << std::endl;  // Debug
  fSmearer->ApplyParameterisation( ntracks, tracks, *td, true );
  //std::cout << "FastSimProcess::Eloss : --- End ---" << std::endl;  // Debug
}


void FastSimProcess::PostStepTypeOfIntrActSampling( Material_t * /* mat */, int /* ntracks */,
                                                    GeantTrack_v & /* tracks */,
                                                    GeantTaskData * /* td */ ) {
  //std::cout << "FastSimProcess::PostStepTypeOfIntrActSampling : Start & End" << std::endl;  // Debug
}


void FastSimProcess::PostStepFinalStateSampling( Material_t * /* mat */, int ntracks,
                                                 GeantTrack_v &tracks, int & /* nout */,
                                                 GeantTaskData *td ) {
  //std::cout << "FastSimProcess::PostStepTypeOfIntrActSampling : Start" << std::endl;  // Debug
  fSmearer->ApplyParameterisation( ntracks, tracks, *td );
  //std::cout << "FastSimProcess::PostStepTypeOfIntrActSampling : --- End ---" << std::endl;  // Debug
}
