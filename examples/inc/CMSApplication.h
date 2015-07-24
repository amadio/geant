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
#include <mutex>
#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

class TH1F;
class TH1D;
class TProfile;

#include "Geant/Config.h"
#include "Geant/Typedefs.h"

/** @brief CMSApplication class */
class CMSApplication : public GeantVApplication {
  static const Int_t kMaxThreads = 36;
  static const Int_t kNvolumes     = 4156;
  static const Int_t kNECALModules = 36;
  static const Int_t kNHCALModules = 100;

public:
enum EScoreType {
  kNoScore = 0,
  kScore
};  

private:
  Bool_t fInitialized;                            /** Initialized flag */
  Bool_t  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  Float_t fEdepECAL[kNECALModules][kMaxThreads];  /** Energy deposition in ECAL */
  Float_t fEdepHCAL[kNHCALModules][kMaxThreads];  /** Energy deposition in HCAL */
  Int_t fECALid[kNECALModules];                   /** ECAL volume id's */
  Int_t fHCALid[kNHCALModules];                   /** HCAL volume id's */
  std::map<int,int> fECALMap;                     /** Map of ECAL modules */
  std::map<int,int> fHCALMap;                     /** Map of ECAL modules */
  std::mutex fMHist;                              /** Mutex for concurrent histogram filling */
  EScoreType fScore;                              /** Entity for scoring */
  TH1F   *fFluxElec;                              /** Flux histogram for electrons */
  TH1F   *fFluxGamma;                             /** Flux histogram for gammas */
  TH1F   *fFluxP;                                 /** Flux histogram for protons */
  TH1F   *fFluxPi;                                /** Flux histogram for pions */
  TH1F   *fFluxK;                                 /** Flux histogram for kaons */
  TH1F   *fEdepElec;                              /** Edep histogram for electrons */
  TH1F   *fEdepGamma;                             /** Edep histogram for gammas */
  TH1F   *fEdepP;                                 /** Edep histogram for protons */
  TH1F   *fEdepPi;                                /** Edep histogram for pions */
  TH1F   *fEdepK;                                 /** Edep histogram for kaons */
  
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

  /** @brief Initialization function */
  virtual Bool_t Initialize();

  /** @brief Set scoring type */
  void SetScoreType(EScoreType type) { fScore = type; }
  
  /**
   * @brief Callback function for user scoring 
   * 
   * @param tid Thread id.
   * @param npart Number of tracks
   * @param tracks GeantV track container
   */
  virtual void StepManager(Int_t npart, const GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Function of digitization
   * 
   * @param event Event that should be digitized
   */
  virtual void Digitize(Int_t event);

  /** @brief User FinishRun function */
  virtual void FinishRun();

  ClassDef(CMSApplication, 1) // User application
};
#endif
