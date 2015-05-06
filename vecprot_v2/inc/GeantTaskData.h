//===--- GeantTaskData.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTaskData.h
 * @brief Implementation of data organized per thread Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKDATA
#define GEANT_TASKDATA

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_TRACK
#include "GeantTrack.h"
#endif

#include <deque>

class TGeoVolume;
class TRandom;
class GeantBasketMgr;
class GeantBasket;

/**
 * @brief Class GeantTaskData
 * @details Class descripting data organized per thread
 * 
 */
class GeantTaskData : public TObject {
public:
  Int_t fTid;          /** Thread unique id */
  Int_t fNthreads;     /** Number of transport threads */
  Int_t fMaxDepth;     /** Maximum geometry depth */
  Int_t fSizeBool;     /** Size of bool array */
  Int_t fSizeDbl;      /** Size of dbl array */
  Bool_t fToClean;     /** Flag set when the basket queue is to be cleaned */
  TGeoVolume *fVolume; /** Current volume per thread */
  TRandom *fRndm;      /** Random generator for thread */
  Bool_t *fBoolArray;  /** [fSizeBool] Thread array of bools */
  Double_t *fDblArray; /** [fSizeDbl] Thread array of doubles */
  GeantTrack fTrack;   /** Track support for this thread */
  VolumePath_t *fPath; /** Volume path for the thread */
  GeantBasketMgr *fBmgr; /** Basket manager collecting mixed tracks */
  std::deque<GeantBasket*> fPool; /** Pool of empty baskets */  

public:

  /** @brief GeantTaskData constructor */
  GeantTaskData();

  /** @brief GeantTaskData destructor */
  virtual ~GeantTaskData();

  /**
   * @brief Function that return double array
   * 
   * @param size Size of double array
   */
  Double_t *GetDblArray(Int_t size);

  /**
   * @brief Function that return boolean array
   * 
   * @param size Size of boolean array
   */
  Bool_t *GetBoolArray(Int_t size);
  
  /**
   * @brief Get the cleared storedtrack
   */
  GeantTrack &GetTrack() { fTrack.Clear(); return fTrack; }
   

  /**
   * @brief Get next free basket or null if not available
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket();

  /** @brief Setter for the toclean flag */
  void SetToClean(Bool_t flag) { fToClean = flag; }

  /** @brief Getter for the toclean flag */
  bool NeedsToClean() const { return fToClean; }

  /**
   * @brief Recycles a given basket
   * 
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b);

  /**
   * @brief Function cleaning a number of free baskets
   * 
   * @param ntoclean Number of baskets to be cleaned
   * @return Number of baskets actually cleaned
   */
  Int_t CleanBaskets(size_t ntoclean);

private:

  /**
   * @brief Constructor GeantTaskData
   * @todo Still not implemented
   */
  GeantTaskData(const GeantTaskData &);

  /**
   * @brief Operator &operator=
   * @todo Still not implemented
   */
  GeantTaskData &operator=(const GeantTaskData &);

  ClassDef(GeantTaskData, 1) // Stateful data organized per thread
};
#endif
