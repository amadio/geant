
#ifndef LHCbFULLAPP_H
#define LHCbFULLAPP_H

#ifndef GEANT_VAPPLICATION
  #include "GeantVApplication.h"
#endif

#include "Geant/Typedefs.h"
#include "GeantFwd.h"
#include "GeantTaskData.h"

namespace GEANT_IMPL_NAMESPACE {
  namespace Geant {
    class GeantRunManager;
    class GeantTaskDataHandle;
    class GeantEvent;
    class GeantTrack;
  }
}


#include "LHCbData.h"

#include "GeantFactory.h"
#include "MyHit.h"

#include <mutex>
#include <vector>

#include "TTree.h"
#include "ROOT/TBufferMerger.hxx"
#include "ROOT/TTaskGroup.hxx"

#include "TFile.h"
#include "TROOT.h"

namespace lhcbapp {

class LHCbParticleGun;

class LHCbFullApp : public Geant::GeantVApplication {
public:


  LHCbFullApp(Geant::GeantRunManager *runmgr, LHCbParticleGun* gun);
  virtual ~LHCbFullApp();

  /** @brief Interface method to allow registration of user defined thread local data. */
  virtual void AttachUserData(Geant::GeantTaskData *td);

  /** @brief Applications creating data per thread have to clean it up */
  virtual void DeleteUserData(Geant::GeantTaskData *td);

  /** @brief Interface method to initialize the application. */
  virtual bool Initialize();

  /** @brief Interace method that is called at the end of each simulation step. */
  virtual void SteppingActions(Geant::GeantTrack &track, Geant::GeantTaskData *td);

  /** @brief Interace method that is called when the transportation of an event (including all primary and their
    *        secondary particles) is completed .*/
  virtual void FinishEvent(Geant::GeantEvent *event);

  /** @brief Interface method that is called at the end of the simulation (when the transportation of all events are
    *        are completed). */
  virtual void FinishRun();


  void         SetPerformanceMode(bool val) { fIsPerformance = val; }

private:
  /** @brief Copy constructor TestEm5 (deleted) */
  LHCbFullApp(const LHCbFullApp&) = delete;
  /** @brief Operator= for TestEm5 (deleted) */
  LHCbFullApp &operator=(const LHCbFullApp&) = delete;


private:
  bool        fIsPerformance;
  bool        fInitialized;
  int         fNumPrimaryPerEvent;
  int         fNumBufferedEvents;
  // user defined thread local data structure handlers to obtain the thread local data structures (defined and
  // registered by the user) during the simulation (in the SteppingActions(i.e. at the end of each simulation step),
  // Digitization(i.e. at the end of an event) and FinishRun (i.e. at the end of the simulation):
  //
  // 1. merged from all working threads when transportation of an event (with all corresponding primary and secondary
  //    particles) are completed
  Geant::TaskDataHandle<LHCbThreadDataEvents>  *fDataHandlerEvents;
  // a unique, run-global user defined data structure to store cumulated quantities per primary particle type
  // during the simulation
  LHCbData   *fData = nullptr;

  // mutex to prevent multiple threads writing into the unique, run-global, unique LHCbData object (in the FinishEvent
  // after the merge of the user defined per-event data distributed among the threads)
  std::mutex fMutex;
  //
  LHCbParticleGun *fGun;

  static const int kMaxThreads = 36;
  static const int kNvolumes     = 4500;
  static const int kNVELOModules = 200;
  static const int kNECALModules = 36;
  static const int kNHCALModules = 112;
  
  bool  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  float fEdepVELO[kNVELOModules][kMaxThreads];  /** Energy deposition in ECAL */
  float fEdepECAL[kNECALModules][kMaxThreads];  /** Energy deposition in ECAL */
  float fEdepHCAL[kNHCALModules][kMaxThreads];  /** Energy deposition in HCAL */
  
  int fVELOid[kNVELOModules];                   /** VELO volume id's */
  int fECALid[kNECALModules];                   /** ECAL volume id's */
  int fHCALid[kNHCALModules];                   /** HCAL volume id's */
  
  std::map<int,int> fVELOMap;                     /** Map of ECAL modules */
  std::map<int,int> fECALMap;                     /** Map of ECAL modules */
  std::map<int,int> fHCALMap;                     /** Map of ECAL modules */

  GeantFactory<MyHit> *fFactory = nullptr;        /** Hits factory */
  ROOT::Experimental::TBufferMerger* fMerger = nullptr;
  int fOutputBlockWrite;
  
};
 
}      // namespace lhcbapp

#endif // LHCbFULLAPP_H
