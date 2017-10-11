//===--- GeantVApplication.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantVApplication.h
 * @brief Implementation of user application in Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VAPPLICATION
#define GEANT_VAPPLICATION

#include "GeantFwd.h"
#include "GeantTrack.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantRunManager;
class GeantEvent;

/** @brief GeantVApplication class */
class GeantVApplication {
public:
  GeantRunManager *fRunMgr; /*taskData*/

  /** @brief GeantVApplication constructor */
  GeantVApplication(GeantRunManager *runmgr);

  void SetRunManager(GeantRunManager *runmgr);

  /** @brief GeantVApplication destructor */
  virtual ~GeantVApplication() {}

  /** @brief Function of initialization */
  virtual bool Initialize() = 0;

  /**
   * @brief User callback function for scoring
   *
   * @param tid  Thread ID
   * @param npart Number of tracks
   * @param tracks Set of tracks
   */
  virtual void StepManager(int /*npart*/, const GeantTrack_v & /*tracks*/, GeantTaskData * /*td*/) {}

  //=== N E W   I N T E R F A C E S ===//

  /**
   * @brief Method called at initialization allowing to attach user data to the
   * task data whiteboard. Use handles provided with TDManager::RegisterUserData
   */
  virtual void AttachUserData(GeantTaskData *) {}

  /** @brief Use TDManager::DeleteUserData providing user data handles */
  virtual void DeleteUserData() {}

  /**
   * @brief Begin a new event.
   * @details The slot number is evt%ninflight, for easier user data management.
   */
  virtual void BeginEvent(int /*evt*/, int /*islot*/) {}
  /**
   * @brief  User FinishEvent function.
   * @details The slot released is evt%ninflight, for easier user data management.
   */
  virtual void FinishEvent(GeantEvent */*event*/) {}

  /** @brief User FinishRun function */
  virtual void FinishRun() {}

  /** @brief Begin new track(s). */
  virtual void BeginTrack(GeantTrack &/*track*/, GeantTaskData */*td*/) {}
  virtual void BeginTrack(TrackVec_t &/*tracks*/, GeantTaskData */*td*/);

  /** @brief Finish track(s). */
  virtual void FinishTrack(GeantTrack &/*track*/, GeantTaskData */*td*/) {}
  virtual void FinishTrack(TrackVec_t &/*tracks*/, GeantTaskData */*td*/);

  /** @brief User stepping actions */
  virtual void SteppingActions(GeantTrack &/*track*/, GeantTaskData */*td*/) {}
  virtual void SteppingActions(TrackVec_t &/*tracks*/, GeantTaskData */*td*/);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
