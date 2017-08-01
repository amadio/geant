//===--- CaloApplicationRP.h - Geant-V Using real physics----------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file CaloApplicationRP.h
 * @brief Implementation of Calo geometry application for Geant-V prototype
 * @author R. Schmitz -- based on example by M. Novak
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CaloApplicationRP
#define GEANT_CaloApplicationRP

#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif

#include "Geant/Typedefs.h"
#include "CaloDetectorConstruction.h"
#include "GeantFwd.h"

#include <vector>
/** @brief CaloApplicationRP class */
class CaloApplicationRP : public Geant::GeantVApplication {
  static const int maxLayers = 100;
  static const int maxAbsorbers = 5;
  int kNlayers;
  int kNumAbsorbers; // lead and scintillator
  using GeantRunManager = Geant::GeantRunManager;
  using CaloDetectorConstruction = Geant::CaloDetectorConstruction;
  using GeantEvent      = Geant::GeantEvent;
  using GeantTrack_v    = Geant::GeantTrack_v;
  using GeantTaskData   = Geant::GeantTaskData;
  using GeantTrack      = Geant::GeantTrack;

public:
  // data structures to store per absorber, per working thread data
  struct DataPerAbsorber {
    std::vector<double> fEdep;   /** Energy deposit in one absorber per layer  */
    std::vector<double> fLength; /** Step length in one absorber per layer */
  };
  struct DataPerThread {
    std::vector<DataPerAbsorber> fListDataPerAbsorber; /** as many as absorbers in this application */
  };
  std::vector<DataPerThread> fListDataPerThread; /** as many as working thread */

private:
  bool fInitialized;                       /** Initialized flag */
//  int  fIdGap[kNlayers];                             /** ID for the gap volume */
  int  fIdAbs[maxLayers][maxAbsorbers];                             /** ID for the absorber volume */

  CaloDetectorConstruction *calo; 
  int previousLayer;
  int  fNumWThreads;                       /** number of working threads */
  // map working thread ID-s to index on [0,#working-thread)
  std::vector<int>  fWThreadIdToIndexMap;  /** as many as thread */

  GeantFactory<MyHit> *fFactory;             /** Hits factory */

  /**
   * @brief Copy constructor CaloApplicationRP
   * * @todo Still not implemented
   */
  CaloApplicationRP(const CaloApplicationRP &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  CaloApplicationRP &operator=(const CaloApplicationRP &);
public:

  /** @brief Constructor CaloApplicationRP */
  CaloApplicationRP(GeantRunManager *runmgr, CaloDetectorConstruction *userCalo);

  /** @brief Destructor ExN03ApplicationRP */
  virtual ~CaloApplicationRP() {}

  /**
   * @brief Function of initialization
   */
  virtual bool Initialize();

  /**
   * @brief Function that provides step manager
   *
   * @param tid ?????
   * @param npart ?????
   * @param tracks GeantV tracks
   */
  virtual void StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Function of digitization
   *
   * @param event Event that should be digitized
   */
  virtual void Digitize(GeantEvent * /*event*/) {}

  /** @brief User FinishRun function */
  virtual void FinishRun();

  //
  // NEW INTERFACE METHODS WITH V3

  // (V3 scalar version)
  virtual void SteppingActions(GeantTrack &track, GeantTaskData *td);
  //

};
#endif
