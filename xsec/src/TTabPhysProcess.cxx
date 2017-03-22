// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "TTabPhysProcess.h"

#include "TTabPhysMgr.h"
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
using std::string;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
TTabPhysProcess::TTabPhysProcess()
  : PhysicsProcessOld(), fMgr(0), fXsecFileName(), fFinalSFileName() {
  // I/O ctor
  SetType(kDiscrete);
}

//______________________________________________________________________________
TTabPhysProcess::TTabPhysProcess(const char *name, const char *fxsec, const char *ffstate)
    : PhysicsProcessOld(name), fMgr(0), fXsecFileName(fxsec), fFinalSFileName(ffstate) {
  // Normal ctor
  SetType(kDiscrete);
}

//______________________________________________________________________________
void TTabPhysProcess::Initialize() {
  // Initialize physics.
   fMgr = TTabPhysMgr::Instance(fXsecFileName.c_str(), fFinalSFileName.c_str());
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
void TTabPhysProcess::ApplyMsc(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/,
                               GeantTaskData * /*td*/) {
  // Apply multiple scattering
  //   fMgr->ApplyMsc(mat, ntracks, tracks, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::ApplyMsc(Material_t * /*mat*/, TrackVec_t & /*tracks*/,
                               GeantTaskData * /*td*/) {
  // Apply multiple scattering
  //   fMgr->ApplyMsc(mat, tracks, td);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
void TTabPhysProcess::Eloss(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout, GeantTaskData *td) {
  // Fill energy loss for the tracks according their fStepV

  nout = fMgr->Eloss(mat, ntracks, tracks, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::Eloss(GeantTrack *track, int &nout, TrackVec_t &output, GeantTaskData *td) {
  // Fill energy loss for the tracks according their fStepV

  nout = fMgr->Eloss(track, output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::Eloss(TrackVec_t &tracks, int &nout, TrackVec_t &output, GeantTaskData *td) {
  // Fill energy loss for the tracks according their fStepV

  nout = fMgr->Eloss(tracks, output, td);
}

//______________________________________________________________________________
void TTabPhysProcess::ComputeIntLen(Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                    GeantTaskData *td) {
  // Tabulated cross section generic process computation of interaction length.

  fMgr->ProposeStep(mat, ntracks, tracks, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::ComputeIntLen(GeantTrack *track, GeantTaskData *td) {
  // Tabulated cross section generic process computation of interaction length.

  fMgr->ProposeStep(track, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::ComputeIntLen(TrackVec_t &tracks, GeantTaskData *td) {
  // Tabulated cross section generic process computation of interaction length.

  fMgr->ProposeStep(tracks, td);
}

//______________________________________________________________________________
void TTabPhysProcess::PostStepTypeOfIntrActSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                                    GeantTaskData *td) {
  // # smapling: target atom and type of the interaction for each primary tracks
  //             all inf. regarding output of sampling is stored in the tracks
  int imat = -1;
  if (mat)
    imat = mat->GetIndex();
  fMgr->SampleTypeOfInteractions(imat, ntracks, tracks, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::PostStepTypeOfIntrActSampling(GeantTrack *track, GeantTaskData *td) {
  // # smapling: target atom and type of the interaction for each primary tracks
  //             all inf. regarding output of sampling is stored in the tracks
  fMgr->SampleTypeOfInteractions(track, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::PostStepTypeOfIntrActSampling(TrackVec_t &tracks, GeantTaskData *td) {
  // # smapling: target atom and type of the interaction for each primary tracks
  //             all inf. regarding output of sampling is stored in the tracks
  fMgr->SampleTypeOfInteractions(tracks, td);
}

//______________________________________________________________________________
void TTabPhysProcess::PostStepFinalStateSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                                 GeantTaskData *td) {
  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled in SampleTypeOfInteractionsInt;
  // # upadting primary track properties and inserting secondary tracks;
  // # return: number of inserted secondary tracks
  int imat = -1;
  if (mat)
    imat = mat->GetIndex();
  nout = fMgr->SampleFinalStates(imat, ntracks, tracks, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::PostStepFinalStateSampling(GeantTrack *track, int &nout,
                                                 TrackVec_t &output, GeantTaskData *td) {
  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled in SampleTypeOfInteractionsInt;
  // # upadting primary track properties and inserting secondary tracks;
  // # return: number of inserted secondary tracks
  nout = fMgr->SampleFinalStates(track, output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TTabPhysProcess::PostStepFinalStateSampling(TrackVec_t &tracks, int &nout,
                                                 TrackVec_t &output, GeantTaskData *td) {
  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled in SampleTypeOfInteractionsInt;
  // # upadting primary track properties and inserting secondary tracks;
  // # return: number of inserted secondary tracks
  nout = fMgr->SampleFinalStates(tracks, output, td);
}

//______________________________________________________________________________
void TTabPhysProcess::AtRest(int /*ntracks*/, GeantTrack_v & /*tracks*/, int & /*nout*/, GeantTaskData * /*td*/) {
  // Do at rest actions on particle after generic tabxsec process.
  // Daughter tracks copied in trackout.
}

//______________________________________________________________________________
void TTabPhysProcess::AtRest(TrackVec_t & /*tracks*/, int & /*nout*/, GeantTaskData * /*td*/) {
  // Do at rest actions on particle after generic tabxsec process.
  // Daughter tracks copied in trackout.
}

} // GEANT_IMPL_NAMESPACE
} // Geant
