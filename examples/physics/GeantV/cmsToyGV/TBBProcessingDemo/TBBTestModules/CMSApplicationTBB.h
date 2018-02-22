//===--- CMSApplicationTBB.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file CMSApplicationTBB.h
 * @brief Implementation of simple scoring for CMS geometry 
 * @author Guilherme Lima, Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CMSApplication
#define GEANT_CMSApplication

#include "tbb/task.h"
#include <mutex>
#include <map>
#include "CMSFullApp.h"
#include "Geant/Config.h"
#include "Geant/Typedefs.h"

#include "base/Vector.h"

using GeantTaskData = geant::GeantTaskData;
using GeantTrack    = geant::GeantTrack;
using TrackVec_t    = geant::GeantTrack_v;

namespace cmsapp {
/** @brief CMSApplication class */
class CMSApplicationTBB : public CMSFullApp {

private:
  std::map<int,tbb::task*> fPostSimTaskMap; /** Map of in-flight event numbers to post-simulation task */
  std::mutex fMapLock;

  CMSApplicationTBB(const CMSApplicationTBB &) = delete;
  CMSApplicationTBB &operator=(const CMSApplicationTBB &) = delete;

public:

  /** @brief Constructor CMSApplicationTBB */
  CMSApplicationTBB(geant::GeantRunManager *runmgr, CMSParticleGun* gun)
    : CMSFullApp(runmgr, gun) {}

  /** @brief Destructor CMSApplicationTBB */
  virtual ~CMSApplicationTBB() {}

  /** @brief Interace method that is called when the transportation of an event (including all primary and their
    *        secondary particles) is completed .*/
  virtual void FinishEvent(geant::GeantEvent *event);

  /**
   * @brief  Receive a pointer to the tbb::task to be run once a given event has been fully simulated.
   * @details Once the simulation of that event is complete, the task's counter gets decremented, which
   *   triggers the task to be spawned.
   */
  void SetEventContinuationTask(int ievt, tbb::task *pTask);
};
}
#endif
