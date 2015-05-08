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

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GeantHitBlock;
class GeantTrack_v;
class GeantTaskData;

/** @brief GeantVApplication class */
class GeantVApplication : public TObject {
public:
  
  /** @brief GeantVApplication constructor */	
  GeantVApplication();

  /** @brief GeantVApplication destructor */
  virtual ~GeantVApplication() {}

  /** @brief Function of initialization */
  virtual Bool_t Initialize() = 0;

  /**
   * @brief User callback function for scoring
   * 
   * @param tid  Thread ID
   * @param npart Number of tracks
   * @param tracks Set of tracks
   */
  virtual void StepManager(Int_t npart, const GeantTrack_v &tracks, GeantTaskData *td) = 0;

  /**
   * @brief Function of digitization
   * 
   * @param event Event for digitization
   */
  virtual void Digitize(Int_t event) = 0;

  /** @brief User FinishRun function */
  virtual void FinishRun() = 0;

  ClassDef(GeantVApplication, 1) // User application
};
#endif
