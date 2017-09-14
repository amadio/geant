#ifndef GEANT_RUN_MANAGER_H
#define GEANT_RUN_MANAGER_H

#include <thread>
#include <atomic>
#include "base/Vector.h"
#include "base/BitSet.h"
#include "Geant/Config.h"
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"
#include "GeantEventServer.h"
#include "GeantConfig.h"
#ifdef USE_ROOT
class TGeoMaterial;
#include "TGeoExtension.h"
#endif


using namespace veccore;
class PhysicsInterface;

GEANT_DEVICE_DECLARE_CONV(Geant,class,GeantPropagator);

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantPropagator;
class TaskBroker;
class PhysicsProcessOld;
class GeantVApplication;
class GeantVDetectorConstruction;
class UserFieldConstruction;
class GeantVTaskMgr;
class GeantEventServer;
class GeantEvent;
class PrimaryGenerator;
class MCTruthMgr;
class EventSet;

class GeantRunManager
{
private:
  bool fInitialized = false;
  int fNpropagators = 0;          /** Number of propagators */
  int fNthreads     = 0;          /** Number of threads per propagator */
  int fNvolumes     = 0;          /** Number of active volumes in the geometry */
  int fNbuff        = 0;          /** Number of event slots per propagator */
  //int fNfeedProp    = 0;          /** Number of propagators with initial feed */
  int fInitialShare = 0;          /** Initial basket share for each propagator */
  GeantConfig *fConfig = nullptr; /** Run configuration */
  TaskBroker *fBroker = nullptr;  /** Task broker */

  GeantVApplication *fApplication = nullptr;    /** User application */
  GeantVApplication *fStdApplication = nullptr; /** Standard application */
  GeantVDetectorConstruction *fDetConstruction = nullptr; /** User detector construction */
  // UserFieldConstruction *fFieldConstruction = nullptr; /** User class to create field */
  
  GeantVTaskMgr     *fTaskMgr = nullptr;  /** GeantV task manager */
  PhysicsProcessOld *fProcess = nullptr;  /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcessOld *fVectorPhysicsProcess = nullptr; /** Interface to vector physics final state sampling */
  PhysicsInterface *fPhysicsInterface = nullptr; /** The new, real physics interface */
  PrimaryGenerator *fPrimaryGenerator = nullptr; /** Primary generator */
  MCTruthMgr *fTruthMgr = nullptr;              /** MCTruth manager */
  GeantEventServer *fEventServer = nullptr;     /** The event server */
  TDManager *fTDManager = nullptr;              /** The task data manager */

  bool fInitialisedRKIntegration=false;  /** Flag: Is RK initialised for tracking in field  */
  float  fBfieldArr[3] = { 0.0, 0.0, 0.0 }; /** Constant Magnetic Field value - if any */

  vector_t<GeantPropagator *> fPropagators;
  vector_t<Volume_t const *> fVolumes;
  vector_t<EventSet *> fEventSets;              /** Registered event sets */
  std::atomic_flag fEventSetsLock;              /** Spinlock for event set access locking */

  // State data
  std::atomic_int fPriorityEvents; /** Number of prioritized events */
  std::atomic_int fTaskId; /** Counter providing unique task id's */
  BitSet *fDoneEvents = nullptr;   /** Array of bits marking done events */
//  int *fNtracks = nullptr;         /** ![fNbuff] Number of tracks per slot */
//  GeantEvent **fEvents = nullptr;  /** ![fNbuff]    Array of events */
//  GeantTaskData **fTaskData = nullptr; /** ![fNthreads] Data private to threads */
  GeantPropagator *fFedPropagator = nullptr; /** Propagator currently being fed */
  std::vector<std::thread> fListThreads; /** Vector of threads */

public:
  GeantRunManager() {}
  GeantRunManager(unsigned int npropagators, unsigned int nthreads, GeantConfig *config);
  ~GeantRunManager();

// Accessors
  void  AddEventSet(EventSet *workload);

  GEANT_FORCE_INLINE
  int  GetNthreads() { return fNthreads; }

  GEANT_FORCE_INLINE
  int  GetNthreadsTotal() { return (fNthreads*fNpropagators); }

  GEANT_FORCE_INLINE
  int  GetNpropagators() { return fNpropagators; }

  GEANT_FORCE_INLINE
  GeantConfig *GetConfig() { return fConfig; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int  GetNvolumes() { return fNvolumes; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNvolumes(int nvol) { fNvolumes = nvol; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  vector_t<Volume_t const *> &GetVolumes() { return fVolumes; }

  GEANT_FORCE_INLINE
  int  GetNprimaries() const { return (fEventServer) ? fEventServer->GetNprimaries() : 0; }

  GEANT_FORCE_INLINE
  int  GetInitialShare() const { return fInitialShare; }

  GEANT_FORCE_INLINE
  GeantVDetectorConstruction *GetDetectorConstruction() const { return fDetConstruction; }

  GEANT_FORCE_INLINE
  void  SetInitialShare(int nbaskets) { fInitialShare = nbaskets; }

  GEANT_FORCE_INLINE
  std::atomic_int &GetPriorityEvents() { return fPriorityEvents; }

  GEANT_FORCE_INLINE
  int GetNpriority() const { return fPriorityEvents.load(); }

  GEANT_FORCE_INLINE
  Volume_t const *GetVolume(int ivol) { return fVolumes[ivol]; }


  GEANT_FORCE_INLINE
  GeantEvent *GetEvent(int slot) { return fEventServer->GetEvent(slot); }

  GEANT_FORCE_INLINE
  GeantEventServer *GetEventServer() const { return fEventServer; }

//  GEANT_FORCE_INLINE
//  GeantTaskData *GetTaskData(int tid) { return fTaskData[tid]; }

  GEANT_FORCE_INLINE
  int  GetTaskId() { return (fTaskId.fetch_add(1)); }

//  GEANT_FORCE_INLINE
//  int GetNtracks(int islot) const { return fNtracks[islot]; }

  GeantPropagator *GetIdlePropagator() const;

  GEANT_FORCE_INLINE
  GeantPropagator *GetFedPropagator() const { return fFedPropagator; }

  GEANT_FORCE_INLINE
  void SetCoprocessorBroker(TaskBroker *broker) { fBroker = broker; }

  GEANT_FORCE_INLINE
  void SetUserApplication(GeantVApplication *app) { fApplication = app; }

  /** @brief Set object to initialize detector, field */
  void SetUserFieldConstruction(UserFieldConstruction* udc);
    
  GEANT_FORCE_INLINE
  GeantVApplication *GetUserApplication() const { return fApplication; }

  GEANT_FORCE_INLINE
  void SetDetectorConstruction(GeantVDetectorConstruction *det) { fDetConstruction = det; }

  GEANT_FORCE_INLINE
  void SetTaskMgr(GeantVTaskMgr *taskmgr) { fTaskMgr = taskmgr; }

  GEANT_FORCE_INLINE
  void SetPhysicsInterface(PhysicsInterface *interface) { fPhysicsInterface = interface; }

  GEANT_FORCE_INLINE
  PhysicsInterface *GetPhysicsInterface() const { return fPhysicsInterface; }

  GEANT_FORCE_INLINE
  void SetPhysicsProcess(PhysicsProcessOld *proc) { fProcess = proc; }

  GEANT_FORCE_INLINE
  void SetVectorPhysicsProcess(PhysicsProcessOld *proc) { fVectorPhysicsProcess = proc; }

  GEANT_FORCE_INLINE
  PrimaryGenerator *GetPrimaryGenerator() const { return fPrimaryGenerator; }

  GEANT_FORCE_INLINE
  void SetPrimaryGenerator(PrimaryGenerator *gen) { fPrimaryGenerator = gen; }

  GEANT_FORCE_INLINE
  void SetMCTruthMgr(MCTruthMgr *mcmgr) { fTruthMgr = mcmgr; }

  GEANT_FORCE_INLINE
  MCTruthMgr *GetMCTruthMgr() const { return fTruthMgr; }

  GEANT_FORCE_INLINE
  TDManager *GetTDManager() const { return fTDManager; }

  GEANT_FORCE_INLINE
  bool IsInitialized() { return fInitialized; }
 
  GeantTaskData *BookTransportTask();

  /** @brief Function checking if transport is completed */
  bool TransportCompleted() const { return ((int)fDoneEvents->FirstNullBit() >= fConfig->fNtotal); }

  /** @brief Check if event is finished */
  bool IsDoneEvent(int ievt) { return fDoneEvents->TestBitNumber(ievt); }

  /** @brief Mark an event as finished */
  void SetDoneEvent(int ievt) { fDoneEvents->SetBitNumber(ievt); }

  /** @brief Initialize classes for RK Integration */
  void PrepareRkIntegration();

  /** @brief Implementation of work stealing */
  int ProvideWorkTo(GeantPropagator *prop);

  void LockEventSets() { while (fEventSetsLock.test_and_set(std::memory_order_acquire)) {}; }
  void UnlockEventSets() { fEventSetsLock.clear(std::memory_order_release); }
 
  EventSet *NotifyEventSets(GeantEvent *finished_event);

  void EventTransported(GeantEvent *event, GeantTaskData *td);
  bool Initialize();
  bool FinishRun();
  bool LoadGeometry(const char *filename);
  void RunSimulation();
  bool RunSimulationTask(EventSet *workload, GeantTaskData *td);
  void StopTransport();

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_RUN_MANAGER_H
