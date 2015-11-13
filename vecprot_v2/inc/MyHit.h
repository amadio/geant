//===--- MyHit.h - Geant-V --------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file MyHit.h
 * @brief Implementation of hits in user application for Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_MYHIT
#define GEANT_MYHIT

#ifndef GEANTV_MIC
#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif
#endif

/**
 * @brief Class MyHit
 * 
 */
class MyHit {
public:
  double fX; /** X position */
  double fY; /** Y position */
  double fZ; /** Z position */
  double fEdep; /** Energy loss */
  double fTime; /** Time */
  int fEvent; /* Event number */
  int fTrack; /* Track number */
  int fVolId;   /** Volume Id */
  int fDetId;   /** Replica (segmentation) */
  int fStatus;  /** track status */

  /** @brief MyHit constructor */
  MyHit() : fX(0), fY(0), fZ(0), fEdep(0),
            fTime(0), fEvent(0), fTrack(0), fVolId(0), fDetId(0), fStatus(0) {}
  
  /**
   * @brief MyHit constructor
   * 
   * @param x X position 
   * @param y Y position
   * @param z Z position
   * @param edep Energy loss
   * @param volid Volume Id
   * @param detid Replica (segmentation)
   */
  MyHit(double x, double y, double z, double edep, double time, int event, int track, int volid, int detid);
  
  /** @brief MyHit destructor*/
  ~MyHit() {}
  
  /**
   * @brief Reset function 
   * @details Set default values for MyHit object
   */
  void Reset() {
    fX = 0.;
    fY = 0.;
    fZ = 0.;
    fEdep = 0.;
    fTime = 0.;
    fEvent = 0;
    fTrack = 0;
    fVolId = 0;
    fDetId = 0;
    fStatus = 0;
  }

  /** @brief Function that add user hit */
  void AddHit();
#ifndef GEANTV_MIC
  ClassDefNV(MyHit, 1) // User hit
#endif
};
#endif
