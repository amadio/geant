//===--- MyApplication.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file MyApplication.h
 * @brief Implementation of typical application for Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_MYAPPLICATION
#define GEANT_MYAPPLICATION

#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif

class GeantTrack_v;

/** @brief MyApplication class */
class MyApplication : public GeantVApplication {
  static const Int_t kNlayers = 10;
  static const Int_t kMaxThreads = 36;

private:
  Bool_t fInitialized;                       /** Initialized flag */
  Int_t fIdGap;                              /** ID for the gap volume */
  Int_t fIdAbs;                              /** ID for the absorber volume */
  Float_t fEdepGap[kNlayers][kMaxThreads];   /** Energy deposition per layer */
  Float_t fLengthGap[kNlayers][kMaxThreads]; /** step length in every layer */
  Float_t fEdepAbs[kNlayers][kMaxThreads];   /** Energy deposition per layer */
  Float_t fLengthAbs[kNlayers][kMaxThreads]; /** Step length in every layer */
  GeantFactory<MyHit> *fFactory;             /** Hits factory */
  
  /**
   * @brief Copy constructor MyApplication
   * * @todo Still not implemented
   */
  MyApplication(const MyApplication &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  MyApplication &operator=(const MyApplication &);
public:

  /** @brief Constructor MyApplication */
  MyApplication();

  /** @brief Destructor MyApplication */
  virtual ~MyApplication() {}

  /**
   * @brief Function of initialization
   */
  virtual Bool_t Initialize();

  /**
   * @brief Function that provides step manager 
   * 
   * @param tid ?????
   * @param npart ?????
   * @param tracks GeantV tracks
   */
  virtual void StepManager(Int_t tid, Int_t npart, const GeantTrack_v &tracks);

  /**
   * @brief Function of digitization
   * 
   * @param event Event that should be digitized
   */
  virtual void Digitize(Int_t event);

  ClassDef(MyApplication, 1) // User application
};
#endif
