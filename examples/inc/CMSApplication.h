//===--- CMSApplication.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file CMSApplication.h
 * @brief Implementation of simple scoring for CMS geometry 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CMSApplication
#define GEANT_CMSApplication

#include <map>
#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

class GeantTrack_v;

/** @brief CMSApplication class */
class CMSApplication : public GeantVApplication {
  static const Int_t kMaxThreads = 36;
  static const Int_t kNvolumes     = 4156;
  static const Int_t kNECALModules = 36;
  static const Int_t kNHCALModules = 100;

private:
  Bool_t fInitialized;                            /** Initialized flag */
  Bool_t  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  Float_t fEdepECAL[kNECALModules][kMaxThreads];  /** Energy deposition in ECAL */
  Float_t fEdepHCAL[kNHCALModules][kMaxThreads];  /** Energy deposition in HCAL */
  Int_t fECALid[kNECALModules];                   /** ECAL volume id's */
  Int_t fHCALid[kNHCALModules];                   /** HCAL volume id's */
  std::map<int,int> fECALMap;                     /** Map of ECAL modules */
  std::map<int,int> fHCALMap;                     /** Map of ECAL modules */
  
  /**
   * @brief Copy constructor CMSApplication
   * * @todo Still not implemented
   */
  CMSApplication(const CMSApplication &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  CMSApplication &operator=(const CMSApplication &);
public:

  /** @brief Constructor CMSApplication */
  CMSApplication();

  /** @brief Destructor CMSApplication */
  virtual ~CMSApplication() {}

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

  ClassDef(CMSApplication, 1) // User application
};
#endif
