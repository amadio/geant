//===--- GeantThreadData.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantThreadData.h
 * @brief Implementation of data organized per thread Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_THREADDATA
#define GEANT_THREADDATA

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
 * @brief Class GeantThreadData
 * @details Class descripting data organized per thread
 * 
 */
class GeantThreadData : public TObject {
public:
  Int_t fTid;          /** Thread unique id */
  Int_t fNthreads;     /** Number of transport threads */
  Int_t fMaxDepth;     /** Maximum geometry depth */
  Int_t fSizeBool;     /** Size of bool array */
  Int_t fSizeDbl;      /** Size of dbl array */
  TGeoVolume *fVolume; /** Current volume per thread */
  TRandom *fRndm;      /** Random generator for thread */
  Bool_t *fBoolArray;  /** [fSizeBool] Thread array of bools */
  Double_t *fDblArray; /** [fSizeDbl] Thread array of doubles */
  GeantTrack fTrack;   /** Track support for this thread */
  VolumePath_t *fPath; /** Volume path for the thread */
  GeantBasketMgr *fBmgr; /** Basket manager collecting mixed tracks */
  std::deque<GeantBasket*> fPool; /** Pool of empty baskets */  

public:

  /** @brief GeantThreadData constructor */
  GeantThreadData();

  /** @brief GeantThreadData destructor */
  virtual ~GeantThreadData();

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
   * @brief Get next free basket or create one with requested capacity
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket(Int_t capacity);

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
   */
  void CleanBaskets(size_t ntoclean);

private:

  /**
   * @brief Constructor GeantThreadData
   * @todo Still not implemented
   */
  GeantThreadData(const GeantThreadData &);

  /**
   * @brief Operator &operator=
   * @todo Still not implemented
   */
  GeantThreadData &operator=(const GeantThreadData &);

  ClassDef(GeantThreadData, 1) // Stateful data organized per thread
};
#endif
