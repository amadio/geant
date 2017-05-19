//===--- TestEM3ApplicationRP.h - Geant-V Using real physics----------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TestEM3ApplicationRP.h
 * @brief Implementation of Geant4 TestEM3 like application for Geant-V prototype
 * @author M. Novak
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TestEM3ApplicationRP
#define GEANT_TestEM3ApplicationRP

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

/** @brief TestEM3ApplicationRP class */
class TestEM3ApplicationRP : public Geant::GeantVApplication {
  static const int kNlayers      = 15;
  static const int kNumAbsorbers = 2; // lead and scintillator
  using GeantRunManager = Geant::GeantRunManager;
  using GeantEvent      = Geant::GeantEvent;
  using GeantTrack_v    = Geant::GeantTrack_v;
  using GeantTaskData   = Geant::GeantTaskData;
  using GeantTrack      = Geant::GeantTrack;

public:
  // data structures to store per absorber, per working thread data
  struct DataPerAbsorber {
    double fEdep;   /** Energy deposit in one absorber per layer  */
    double fLength; /** Step length in one absorber per layer */
  };
  struct DataGloabl {
    double fNumGamma;
    double fNumElectron;
    double fNumPositron;
    double fNumChargedSteps;
    double fNumNeutralSteps;
  };
  struct DataPerThread {
    std::vector<DataPerAbsorber> fListDataPerAbsorber; /** as many as absorber (2) in this application */
    DataGloabl                   fDataGlobal;
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
   * @brief Copy constructor TestEM3ApplicationRP
   * * @todo Still not implemented
   */
  TestEM3ApplicationRP(const TestEM3ApplicationRP &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  TestEM3ApplicationRP &operator=(const TestEM3ApplicationRP &);
public:

  /** @brief Constructor TestEM3ApplicationRP */
  TestEM3ApplicationRP(GeantRunManager *runmgr);

  /** @brief Destructor TestEM3ApplicationRP */
  virtual ~TestEM3ApplicationRP() {}

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
