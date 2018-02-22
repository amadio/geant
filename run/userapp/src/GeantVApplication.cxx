#include "GeantVApplication.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantVApplication::GeantVApplication(RunManager *runmgr):fRunMgr(runmgr)
{
  // Ctor..

}

//______________________________________________________________________________
void GeantVApplication::SetRunManager(RunManager *runmgr) {
  fRunMgr = runmgr;
}

//______________________________________________________________________________
void GeantVApplication::BeginTrack(TrackVec_t &tracks, GeantTaskData *td)
{
  // Invoke by default the scalar version
  for (auto track : tracks)
    BeginTrack(*track, td);
}

//______________________________________________________________________________
void GeantVApplication::FinishTrack(TrackVec_t &tracks, GeantTaskData *td)
{
  // Invoke by default the scalar version
  for (auto track : tracks)
    FinishTrack(*track, td);
}

//______________________________________________________________________________
void GeantVApplication::SteppingActions(TrackVec_t &tracks, GeantTaskData *td)
{
  // Invoke by default the scalar version
  for (auto track : tracks)
    SteppingActions(*track, td);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
