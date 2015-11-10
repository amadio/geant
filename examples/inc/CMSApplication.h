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
#ifdef USE_VECGEOM_NAVIGATOR
#include "base/Map.h"
#else
#include <map>
#endif

#include <mutex>
#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

class TH1F;
class TH1D;
class TProfile;

#include "Geant/Config.h"
#include "Geant/Typedefs.h"


#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif


/** @brief CMSApplication class */
class CMSApplication : public GeantVApplication {
  static const int kMaxThreads = 36;
  static const int kNvolumes     = 4500;
  static const int kNECALModules = 36;
  static const int kNHCALModules = 112;

public:
enum EScoreType {
  kNoScore = 0,
  kScore
};  

private:
  bool fInitialized;                            /** Initialized flag */
  bool  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  float fEdepECAL[kNECALModules][kMaxThreads];  /** Energy deposition in ECAL */
  float fEdepHCAL[kNHCALModules][kMaxThreads];  /** Energy deposition in HCAL */
  int fECALid[kNECALModules];                   /** ECAL volume id's */
  int fHCALid[kNHCALModules];                   /** HCAL volume id's */
#ifdef USE_VECGEOM_NAVIGATOR 
  vecgeom::map<int,int> fECALMap;                     /** Map of ECAL modules */
  vecgeom::map<int,int> fHCALMap;                     /** Map of ECAL modules */
#else 
  std::map<int,int> fECALMap;                     /** Map of ECAL modules */
  std::map<int,int> fHCALMap;                     /** Map of ECAL modules */
#endif
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
  GeantFactory<MyHit> *fFactory;                  /** Hits factory */

  
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
  virtual bool Initialize();

  /** @brief Set scoring type */
  void SetScoreType(EScoreType type) { fScore = type; }
  
  /**
   * @brief Callback function for user scoring 
   * 
   * @param tid Thread id.
   * @param npart Number of tracks
   * @param tracks GeantV track container
   */
  virtual void StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Function of digitization
   * 
   * @param event Event that should be digitized
   */
  virtual void Digitize(int event);

  /** @brief User FinishRun function */
  virtual void FinishRun();

  ClassDef(CMSApplication, 1) // User application
};
#endif
