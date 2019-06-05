#include "Geant/Event.h"
#include "Geant/RunManager.h"
#include <iostream>

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
int Event::AddTrack()
{
  // Thread safe track addition
  int ntracks   = ++fNtracks;
  int ninflight = ntracks - fNdone.load();
  //  int nmax = fNmax.load();
  // The fNmax does not need to be accurate, avoid expensive compare_exchange
  if (fNmax < ninflight) fNmax.store(ninflight);
  // Thread safe update of the max number of tracks in flight
  //  while ( fNmax < ninflight && !fNmax.compare_exchange_weak(nmax, ninflight) )
  //    ;
  return (ntracks - 1);
}

//______________________________________________________________________________
void Event::Clear(TaskData *td)
{
  // Clear the event.
  fPrioritize  = false;
  fTransported = false;
  fEvent       = 0;
  fSlot        = 0;
  fNprimaries  = 0;
  fNtracks.store(0);
  fNdone.store(0);
  fNmax.store(0);
  fNmax.store(0);
  fLock.clear();
  fNfilled.store(0);
  // Release primary tracks
  if (td)
    for (auto track : fPrimaries)
      td->ReleaseTrack(*track);
  fPrimaries.clear();
  fNdispatched.store(0);
  fTransported = false;
}

//______________________________________________________________________________
bool Event::StopTrack(RunManager *runmgr, TaskData *td)
{
// Mark one track as stopped. Check if event has to be prioritized and return
// true in this case.
#ifdef VECCORE_CUDA
  assert(0 && "StopTrack not implemented yet for CUDA host/device code.");
#else
  // td->fNinflight--;
  int ndone     = fNdone.fetch_add(1) + 1;
  int ntracks   = fNtracks.load();
  int ninflight = ntracks - ndone;
  if ((ndone >= fNprimaries) && (ninflight == 0)) {
    fTransported = true;
    // Notify run manager that event is transported
    runmgr->EventTransported(this, td);
    FinalActions();
    return false;
  }
  if (!fPrioritize) {
    // Check if the event needs to be prioritized
    int npriority = runmgr->GetNpriority();
    if (npriority < runmgr->GetNthreads()) {
      if (ninflight < fPriorityThr * GetNmax()) {
        fPrioritize = true;
        //      std::cout << "### Event " << fEvent << " prioritized at " <<
        //        100.*fPriorityThr << " % threshold (npri=" << npriority << ")" << std::endl;
        return true;
      }
    }
  }
#endif
  return false;
}

//______________________________________________________________________________
void Event::Print(const char *) const
{
  // Print events content
  std::cout << "Event " << fEvent << ": " << GetNtracks() << " tracks transported, max in flight " << GetNmax()
            << std::endl;
}

//______________________________________________________________________________
bool Event::Prioritize()
{
  // Prioritize the event
  if (fLock.test_and_set(std::memory_order_acquire) || fPrioritize) return false;
  if (GetNinflight()) {
    //    std::cout << "### Event " << fEvent << " forced prioritized" << std::endl;
    fPrioritize = true;
  }
  fLock.clear(std::memory_order_release);
  return true;
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
