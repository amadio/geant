//===--- MCTruthMgr.h - Geant-V ---------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file MCTruthMgr.h
 * @brief Implementation of MCTruthMgr for Geant-V prototype
 * @author Witek Pokorski
 */
//===----------------------------------------------------------------------===//

#ifndef GEANTV_MCTruthMgr_h
#define GEANTV_MCTruthMgr_h

#include "base/Global.h"

#include "Geant/Geant/Fwd.h"
#include "Track.h"

#include "cuckoohash_map.hh"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

struct MCParticle
{
  int        pid;        /** PDG ID */
  
  int        motherid;   /** ID of the mother particle */

  double fXpos;          /** X position */
  double fYpos;          /** Y position */
  double fZpos;          /** Z position */
  double fTime;          /** Time */
  
  double fPx;            /** X direction */
  double fPy;            /** Y direction */
  double fPz;            /** Z direction */
  double fE;             /** Energy */

  bool has_end;          /** flag saying if end vertex exists */
  double fXend;          /** end X position */
  double fYend;          /** end Y position */
  double fZend;          /** end Z position */
  double fTend;          /** end Time */
};

struct MCEvent
{
  int event_id;
  cuckoohash_map<int, MCParticle*> particles;
};

/**
 * @brief Class of MC truth manager
 */
class MCTruthMgr {

protected:
  cuckoohash_map<int, MCEvent*> events_map;

public:
  MCTruthMgr() {}
  virtual ~MCTruthMgr() {}

  /**
   * @brief Pure virtual function of initialization of truth manager
   * @details Set truth handling properties
   */
  virtual void InitMCTruthMgr() = 0;

  /**
   * @brief Pure virtual function checking the conditions of the track to be stored
   * @details Looks at different tracks properties and returns yes/no to be stored
   */
  virtual bool CheckTrack(geant::Track &gtrack, MCEvent* evt) = 0;

  /**
   * @brief Function that adds a track
   *
   * @param gtrack track
   */
  void AddTrack(geant::Track &gtrack);

  /**
   * @brief Function that indicates a track stopping
   *
   * @param track Track to be stopped
   */
  void EndTrack(Track *track);

    /**
   * @brief Function that opens an event
   *
   * @param evID ID of event
   */
  void OpenEvent(int evID);

  /**
   * @brief Pure virtual function that closes an event
   *
   * @param evID ID of event
   */
  virtual void CloseEvent(int evID) = 0;
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
