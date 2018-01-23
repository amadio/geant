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
#ifdef VECCORE_CUDA
#include "base/Map.h"
#else
#include <map>
#endif

#include <mutex>
#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif
#ifdef USE_ROOT
class TH1F;
class TH1D;
class TProfile;
#endif
#include "Geant/Config.h"
#include "Geant/Typedefs.h"


#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif

#include "base/Vector.h"

/** @brief CMSApplication class */
class CMSApplication : public Geant::GeantVApplication {
  static const int kMaxThreads = 36;
  static const int kNvolumes     = 4500;
  static const int kNECALModules = 36;
  static const int kNHCALModules = 112;
  using GeantRunManager = Geant::GeantRunManager;
  using GeantEvent = Geant::GeantEvent;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

public:

enum EScoreType {
  kNoScore = 0,
  kScore
};  

private:
#if defined(USE_VECGEOM_NAVIGATOR) && defined(VECCORE_CUDA)
  using Map_t = vecgeom::map<int,int>;
#else
  using Map_t = std::map<int,int>;
#endif
  bool fInitialized;                            /** Initialized flag */
  bool  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  float fEdepECAL[kNECALModules][kMaxThreads];  /** Energy deposition in ECAL */
  float fEdepHCAL[kNHCALModules][kMaxThreads];  /** Energy deposition in HCAL */
  int fECALid[kNECALModules];                   /** ECAL volume id's */
  int fHCALid[kNHCALModules];                   /** HCAL volume id's */
  Map_t fECALMap;                               /** Map of ECAL modules */
  Map_t fHCALMap;                                /** Map of ECAL modules */
  std::mutex fMHist;                              /** Mutex for concurrent histogram filling */
  EScoreType fScore;                              /** Entity for scoring */
#ifdef USE_ROOT
  TH1F   *fFluxElec;                              /** Flux histogram for electrons */
  TH1F   *fFluxGamma;                             /** Flux histogram for gammas */
  TH1F   *fFluxP;                                 /** Flux histogram for protons */
  TH1F   *fFluxPi;                                /** Flux histogram for pions */
  TH1F   *fFluxK;                                 /** Flux histogram for kaons */
  TH1F   *fEdepElec;                              /** Edep histogram for electrons */
  TH1F   *fEdepGamma;                             /** Edep histogram for gammas */
  TH1F   *fEdepP;                                 /** Edep histogram for protons */
  TH1F   *fEdepPi;                                /** Edep histogram for pions */
  TH1F   *fEdepK; 
#endif                                /** Edep histogram for kaons */
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
  CMSApplication(GeantRunManager *runmgr);

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

  /** @brief  User FinishEvent function.*/
  virtual void FinishEvent(GeantEvent *event);

  /** @brief User FinishRun function */
  virtual void FinishRun();
};
#endif
