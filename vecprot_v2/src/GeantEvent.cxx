#include "GeantEvent.h"
#include "GeantRunManager.h"
#include <iostream>

//______________________________________________________________________________
int GeantEvent::AddTrack() {
  // Thread safe track addition
  int ntracks = ++fNtracks;
  int ninflight = ntracks - fNdone.load();
  int nmax = fNmax.load();
  // Thread safe update of the max number of tracks in flight
  while ( fNmax < ninflight && !fNmax.compare_exchange_weak(nmax, ninflight) )
    ;
  return ntracks;  
}

//______________________________________________________________________________
bool GeantEvent::StopTrack(GeantRunManager *runmgr) {
  // Mark one track as stopped. Check if event has to be prioritized and return
  // true in this case.
  fNdone++;
  int npriority = runmgr->GetNpriority();
  if (!fPrioritize && (npriority < runmgr->GetNthreads())) {
    if (GetNinflight() < fPriorityThr*GetNmax()) {
      fPrioritize = true;
//      std::cout << "### Event " << GetEvent() << " prioritized at " <<
//        100.*fPriorityThr << " % threshold (npri=" << npriority << ")" << std::endl;
      return true;  
    }  
  }
  return false;
}

//______________________________________________________________________________
void GeantEvent::Print(const char *) const {
  // Print events content
  std::cout << "Event " << GetEvent() << ": " << GetNtracks() << 
    " tracks transported, max in flight " <<  GetNmax() << std::endl;
}

//______________________________________________________________________________
bool GeantEvent::Prioritize() {
  // Prioritize the event
  if (fLock.test_and_set(std::memory_order_acquire) || fPrioritize) return false;
  if (GetNinflight()) {
//    std::cout << "### Event " << GetEvent() << " forced prioritized" << std::endl;
    fPrioritize = true;
  }
  fLock.clear(std::memory_order_release);
  return true;
}
  
