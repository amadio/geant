#include "UserApplication.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
UserApplication::UserApplication(RunManager *runmgr):fRunMgr(runmgr)
{
  // Ctor..

}

//______________________________________________________________________________
void UserApplication::SetRunManager(RunManager *runmgr) {
  fRunMgr = runmgr;
}

//______________________________________________________________________________
void UserApplication::BeginTrack(TrackVec_t &tracks, TaskData *td)
{
  // Invoke by default the scalar version
  for (auto track : tracks)
    BeginTrack(*track, td);
}

//______________________________________________________________________________
void UserApplication::FinishTrack(TrackVec_t &tracks, TaskData *td)
{
  // Invoke by default the scalar version
  for (auto track : tracks)
    FinishTrack(*track, td);
}

//______________________________________________________________________________
void UserApplication::SteppingActions(TrackVec_t &tracks, TaskData *td)
{
  // Invoke by default the scalar version
  for (auto track : tracks)
    SteppingActions(*track, td);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
