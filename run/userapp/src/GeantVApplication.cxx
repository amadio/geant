#include "GeantVApplication.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantVApplication::GeantVApplication(GeantRunManager *runmgr):fRunMgr(runmgr)
{
  // Ctor..

}

//______________________________________________________________________________
void GeantVApplication::SetRunManager(GeantRunManager *runmgr) {
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
