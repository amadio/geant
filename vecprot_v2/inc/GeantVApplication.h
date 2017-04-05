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

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantRunManager;
class GeantEvent;

/** @brief GeantVApplication class */
class GeantVApplication {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

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
  virtual void StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) = 0;

  /**
   * @brief Function of digitization
   * 
   * @param event Event for digitization
   */
  virtual void Digitize(GeantEvent *event) = 0;

  /** @brief User FinishRun function */
  virtual void FinishRun() = 0;

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
