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

class TGeoVolume;
class TRandom;
class GeantBasketMgr;

/**
 * @brief Class GeantThreadData
 * @details Class descripting data organized per thread
 * 
 */
class GeantThreadData : public TObject {
public:
  Int_t fMaxPerBasket; /** Max number of tracks per basket */
  Int_t fNprocesses;   /** Number of physics processes */
  Int_t fSizeDbl;      /** Size of array of doubles */
  Int_t fSizeBool;     /** Size of bool array */
  TGeoVolume *fVolume; /** Current volume per thread */
  TRandom *fRndm;      /** Random generator for thread */
  Bool_t *fBoolArray;  /** [5*fMaxPerBasket] Thread array of bools */
  Double_t *fDblArray; /** [5*fMaxPerBasket] Thread array of doubles */
  Double_t *fProcStep; /** [fNprocesses*fMaxPerBasket] */
  GeantTrack fTrack;   /** Track support for this thread */
  VolumePath_t *fPath; /** Volume path for the thread */
  GeantBasketMgr *fBmgr; /** Basket manager collecting mixed tracks */  

public:

  /** @brief GeantThreadData constructor */
  GeantThreadData();

  /**
   * @brief GeantThreadData parameterized constructor
   * 
   * @param maxperbasket Maximum tracks per basket 
   * @param maxprocesses Maximum physics processes 
   */
  GeantThreadData(Int_t maxperbasket, Int_t maxprocesses);

  /** @brief GeantThreadData destructor */
  virtual ~GeantThreadData();

  /**
   * @brief Function that return process step
   * 
   * @param iproc ID process
   * @todo  Doxygen details on return function
   * @return fNprocesses*fMaxPerBasket + iproc*fMaxPerBasket ?????
   */
  Double_t *GetProcStep(Int_t iproc) { return fProcStep + iproc * fMaxPerBasket; }

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
