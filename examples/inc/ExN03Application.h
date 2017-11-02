//===--- ExN03Application.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file ExN03Application.h
 * @brief Implementation of Geant4 example N03 application for Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_ExN03Application
#define GEANT_ExN03Application

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
#include "GeantTaskData.h"
#include "EventSet.h"
#include "ExN03ScoringData.h"

class PrimaryGenerator;

/** @brief ExN03Application class */
class ExN03Application : public Geant::GeantVApplication {
  static const int kNlayers = 15;
  static const int kMaxThreads = 36;
  using GeantRunManager = Geant::GeantRunManager;
  using PrimaryGenerator = Geant::PrimaryGenerator;
  using GeantEvent = Geant::GeantEvent;
  using EventSet = Geant::EventSet;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTrack = Geant::GeantTrack;
  using GeantTaskData = Geant::GeantTaskData;
  using DigitsHandle = Geant::TaskDataHandle<ExN03ScoringData>;

private:
  bool fInitialized;                       /** Initialized flag */
  PrimaryGenerator *fGenerator;            /** Generator used in external loop mode */
  int fIdGap;                              /** ID for the gap volume */
  int fIdAbs;                              /** ID for the absorber volume */
  // These are needed only for v2
  float fEdepGap[kNlayers][kMaxThreads];   /** Energy deposition per layer */
  float fLengthGap[kNlayers][kMaxThreads]; /** step length in every layer */
  float fEdepAbs[kNlayers][kMaxThreads];   /** Energy deposition per layer */
  float fLengthAbs[kNlayers][kMaxThreads]; /** Step length in every layer */

  DigitsHandle *fDigitsHandle = nullptr;   /** Handler for digits */
  
  GeantFactory<MyHit> *fFactory;           /** Hits factory */
  
  /**
   * @brief Copy constructor ExN03Application
   * * @todo Still not implemented
   */
  ExN03Application(const ExN03Application &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  ExN03Application &operator=(const ExN03Application &);
public:

  /** @brief Constructor ExN03Application */
  ExN03Application(GeantRunManager *runmgr);

  /** @brief Destructor ExN03Application */
  virtual ~ExN03Application() {}

  /** @brief In external event loop mode, this aplication can use an internal generator */
  void SetGenerator(PrimaryGenerator *gen);

  /** @brief Generate an event set to be processed by a single task.
       Not required as application functionality, the event reading or generation
       can in the external event loop.
  */
  EventSet* GenerateEventSet(size_t nevents, GeantTaskData *td);

  /**
   * @brief Method called at initialization allowing to attach user data to the
   * task data whiteboard. Called by every worker. in the initialization phase.
   */
  virtual void AttachUserData(GeantTaskData *td);

  /**
   * @brief Function of initialization
   */
  virtual bool Initialize();

  /** @brief User scoring per step */
  virtual void SteppingActions(GeantTrack &/*track*/, GeantTaskData */*td*/);

  /** @brief  User FinishEvent function.*/
  virtual void FinishEvent(GeantEvent *event);

  /** @brief User FinishRun function */
  virtual void FinishRun() {}

};
#endif
