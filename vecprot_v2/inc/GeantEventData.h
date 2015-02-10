//===--- GeantEventData.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEventData.h
 * @brief Implementation of event data organized per thread for Geant-V prototype
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_THREADDATA
#define GEANT_THREADDATA

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include <vector>

class GeantEvent;
class GeantHitBlock;

/** @brief Class GeantEventData responsible for event data organized per thread */
class GeantEventData : public TObject {
public:
  GeantEvent *fEvent;                 /** Event pointer */
  vector<GeantHitBlock *> fHitBlocks; /** Hit blocks */
  // GeantDigitBlock   **fDigitBlocks;           //! Digit blocks
public:
  
  /** @brief GeantEventData constructor */
  GeantEventData();

  /** @brief GeantEventData destructor */
  virtual ~GeantEventData();
  
  /**
   * @brief Function of addition block of hits 
   * 
   * @param block Block of hits
   */
  void AddHitBlock(GeantHitBlock *block);

  /** @brief Function of clearing hits */
  void ClearHits();

  ClassDef(GeantEventData, 1) // Data per event
};
#endif
