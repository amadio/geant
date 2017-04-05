//===--- ExN03ApplicationRP.h - Geant-V Using real physics----------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ExN03ApplicationRP.h
 * @brief Implementation of Geant4 example N03 application for Geant-V prototype
 * @author M. Novak
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_ExN03ApplicationRP
#define GEANT_ExN03ApplicationRP

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

#include "GeantFwd.h"

#include <vector>

/** @brief ExN03ApplicationRP class */
class ExN03ApplicationRP : public Geant::GeantVApplication {
  static const int kNlayers = 15;
  static const int kNumAbsorbers = 2; // lead and scintillator
  using GeantRunManager = Geant::GeantRunManager;
  using GeantEvent      = Geant::GeantEvent;
  using GeantTrack_v    = Geant::GeantTrack_v;
  using GeantTaskData   = Geant::GeantTaskData;

public:
  // data structures to store per absorber, per working thread data
  struct DataPerAbsorber {
    std::vector<double> fEdep;   /** Energy deposit in one absorber per layer  */
    std::vector<double> fLength; /** Step length in one absorber per layer */
  };
  struct DataPerThread {
    std::vector<DataPerAbsorber> fListDataPerAbsorber; /** as many as absorber (2) in this application */
  };
  std::vector<DataPerThread> fListDataPerThread; /** as many as working thread */

private:
  bool fInitialized;                       /** Initialized flag */
  int  fIdGap;                             /** ID for the gap volume */
  int  fIdAbs;                             /** ID for the absorber volume */
  int  fNumWThreads;                       /** number of working threads */
  // map working thread ID-s to index on [0,#working-thread)
  std::vector<int>  fWThreadIdToIndexMap;  /** as many as thread */

  GeantFactory<MyHit> *fFactory;             /** Hits factory */

  /**
   * @brief Copy constructor ExN03ApplicationRP
   * * @todo Still not implemented
   */
  ExN03ApplicationRP(const ExN03ApplicationRP &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  ExN03ApplicationRP &operator=(const ExN03ApplicationRP &);
public:

  /** @brief Constructor ExN03ApplicationRP */
  ExN03ApplicationRP(GeantRunManager *runmgr);

  /** @brief Destructor ExN03ApplicationRP */
  virtual ~ExN03ApplicationRP() {}

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

};
#endif
