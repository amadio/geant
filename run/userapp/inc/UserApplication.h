//===--- UserApplication.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file UserApplication.h
 * @brief Implementation of user application in Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VAPPLICATION
#define GEANT_VAPPLICATION

#include "GeantFwd.h"
#include "Track.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class RunManager;
class Event;

/** @brief UserApplication class */
class UserApplication {
public:
  RunManager *fRunMgr; /*taskData*/

  /** @brief UserApplication constructor */
  UserApplication(RunManager *runmgr);

  void SetRunManager(RunManager *runmgr);

  /** @brief UserApplication destructor */
  virtual ~UserApplication() {}

  /** @brief Function of initialization */
  virtual bool Initialize() = 0;

  //=== N E W   I N T E R F A C E S ===//

  /**
   * @brief Method called at initialization allowing to attach user data to the
   * task data whiteboard. Use handles provided with TDManager::RegisterUserData
   */
  virtual void AttachUserData(TaskData *) {}

  /** @brief Use TDManager::DeleteUserData providing user data handles */
  virtual void DeleteUserData(TaskData *) {}

  /**
   * @brief Begin a new event.
   * @details The slot number is evt%ninflight, for easier user data management.
   */
  virtual void BeginEvent(int /*evt*/, int /*islot*/) {}
  /**
   * @brief  User FinishEvent function.
   * @details The slot released is evt%ninflight, for easier user data management.
   */
  virtual void FinishEvent(Event */*event*/) {}

  /** @brief User FinishRun function */
  virtual void FinishRun() {}

  /** @brief Begin new track(s). */
  virtual void BeginTrack(Track &/*track*/, TaskData */*td*/) {}
  virtual void BeginTrack(TrackVec_t &/*tracks*/, TaskData */*td*/);

  /** @brief Finish track(s). */
  virtual void FinishTrack(Track &/*track*/, TaskData */*td*/) {}
  virtual void FinishTrack(TrackVec_t &/*tracks*/, TaskData */*td*/);

  /** @brief User stepping actions */
  virtual void SteppingActions(Track &/*track*/, TaskData */*td*/) {}
  virtual void SteppingActions(TrackVec_t &/*tracks*/, TaskData */*td*/);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
