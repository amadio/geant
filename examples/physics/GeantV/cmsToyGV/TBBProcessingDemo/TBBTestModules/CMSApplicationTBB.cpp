//===--- CMSApplicationTBB.cpp - Geant-V ------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file CMSApplicationTBB.cpp
 * @brief Implementation of simple scoring for CMS geometry 
 * @author Guilherme Lima, Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#include "CMSApplicationTBB.h"
#include "Geant/Event.h"
#include "TaskData.h"
#include "globals.h"
#include "Geant/Error.h"

using namespace geant;

namespace cmsapp {
//// User actions in terms of TBB tasks
void CMSApplicationTBB::SetEventContinuationTask(int ievt, tbb::task *pTask) {
  std::lock_guard<std::mutex> lock(fMapLock);
  fPostSimTaskMap.insert( std::pair<int,tbb::task*>(ievt, pTask) );
  printf("CMSApplicTBB::SetEvtContTask: ievt=%i, pTask=<%p>, map.size=%lu\n", ievt, pTask, fPostSimTaskMap.size());
}


// //______________________________________________________________________________
// // any user tasks after simulation is complete (and before digitization)
// tbb::task *CMSApplicationTBB::SpawnUserEventFeeder(EventServer *evserv) {

//   printf("-- CMSAppTBB::SpawnUserEventFeeder() called...");
//   //return static_cast<CMSApplicationTBB*>(fRunMgr->GetUserApplication())->SpawnEventLoopTask(/*nevents*/);
//   return 0;
// }

// //______________________________________________________________________________
// // any user tasks after simulation is complete (and before digitization)
//  tbb::task *CMSApplicationTBB::SpawnUserEndRunTask() {
//    // any user tasks after simulation is complete (and before digitization)
//    printf("-- CMSAppTBB::SpawnUserEndRunTask() called...");
//    //return static_cast<CMSApplicationTBB*>(fRunMgr->GetUserApplication())->SpawnEndRunTask();
//    return 0;
// }

//______________________________________________________________________________
void CMSApplicationTBB::FinishEvent(geant::Event *event) {
  // find next tbb::task and decrement its ref count
  CMSFullApp::FinishEvent(event);
  std::lock_guard<std::mutex> lock(fMapLock);
  auto iter = fPostSimTaskMap.find(event->GetEvent());
  if (iter != fPostSimTaskMap.end()) {
    tbb::task* pTask = iter->second;
    printf("CMSAppTBB::FinishEvent(%i,%i), iter=<%p>, map.size=%lu\n", event->GetEvent(), event->GetSlot(), pTask, fPostSimTaskMap.size());
    pTask->decrement_ref_count();
    printf("CMSAppTBB::FinishEvent: pTask ref count=%i\n", pTask->ref_count());
  }
}

}
