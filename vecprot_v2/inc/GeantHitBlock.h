//===--- GeantHitBlock.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTaskManager.h
 * @brief Implementation of a fixed-size block of user hits stored contiguously in memory
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_HITBLOCK
#define GEANT_HITBLOCK

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class TClass;

/**
 * @brief GeantHitBlock class
 * @details Fixed-size block of user hits stored contiguously in memory
 * The class is used to efficiently manage hit user data in a multithreaded environment.
 */
class GeantHitBlock : public TObject {
public:
  Int_t fSize;    /** Fixed size */
  Int_t fNext;    /** number of hits */
  TClass *fClass; /** class type for the hits (MyHit::Class()) */
  Int_t fClSize;  /** Size of hit class stored */
  char *fSupport; /** Support array */
  void *fArray;   /** Hits array */
  
  /**
   * @brief GeantHitBlock constructor
   * 
   * @param cl Hit class
   * @param size Size of hit class
   * @param storage Storage function object
   */
  GeantHitBlock(TClass *cl, Int_t size, void *storage = 0);
  
  /** @brief GeantHitBlock destructor */
  virtual ~GeantHitBlock();
  
  /**
   * @brief Function of checking if block is full
   * @return Boolean value if fNext == fSize
   */
  Bool_t Full() const { return (fNext == fSize); }
  
  /** @brief Function of next hit object */
  TObject *NextHit();
  
  /** @brief Function of cleaning hits */
  void ClearHits();

  /**
   * @brief Function of hit copy
   * 
   * @param ihit Hit for copy 
   * @param other GeantHitBlock object to copy 
   */
  void CopyHit(Int_t ihit, GeantHitBlock &other);

  ClassDef(GeantHitBlock, 0) // A block of hits
};
#endif
