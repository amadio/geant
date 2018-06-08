#ifndef GEANT_RUN_MANAGER_H
#define GEANT_RUN_MANAGER_H

#include <thread>
#include <atomic>
#include "base/Vector.h"
#include "base/BitSet.h"
#include "Geant/Config.h"
#include "Geant/Typedefs.h"
#include "Geant/TaskData.h"
#include "Geant/EventServer.h"
#include "GeantConfig.h"
#ifdef USE_ROOT
class TGeoMaterial;
#include "TGeoExtension.h"
#endif

using namespace veccore;
class PhysicsInterface;

GEANT_DEVICE_DECLARE_CONV(Geant, class, Propagator);

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class Propagator;
class TaskBroker;
class UserApplication;
class UserDetectorConstruction;
class UserFieldConstruction;
class EventServer;
class Event;
class PrimaryGenerator;
class MCTruthMgr;
class EventSet;

class RunManager {
private:
  bool fInitialized    = false;
  int fNpropagators    = 0;       /** Number of propagators */
  int fNthreads        = 0;       /** Number of threads per propagator */
  int fNvolumes        = 0;       /** Number of active volumes in the geometry */
  int fNbuff           = 0;       /** Number of event slots per propagator */
  int fInitialShare    = 0;       /** Initial basket share for each propagator */
  GeantConfig *fConfig = nullptr; /** Run configuration */
  TaskBroker *fBroker  = nullptr; /** Task broker */

  UserApplication *fApplication              = nullptr; /** User application */
  UserApplication *fStdApplication           = nullptr; /** Standard application */
  UserDetectorConstruction *fDetConstruction = nullptr; /** User detector construction */

  PhysicsInterface *fPhysicsInterface = nullptr; /** The new, real physics interface */
  PrimaryGenerator *fPrimaryGenerator = nullptr; /** Primary generator */
  MCTruthMgr *fTruthMgr               = nullptr; /** MCTruth manager */
  EventServer *fEventServer           = nullptr; /** The event server */
  TDManager *fTDManager               = nullptr; /** The task data manager */

  bool fInitialisedRKIntegration = false;           /** Flag: Is RK initialised for tracking in field  */
  float fBfieldArr[3]            = {0.0, 0.0, 0.0}; /** Constant Magnetic Field value - if any */

  vector_t<Propagator *> fPropagators;
  vector_t<Volume_t const *> fVolumes;
  vector_t<EventSet *> fEventSets; /** Registered event sets */
  std::atomic_flag fEventSetsLock; /** Spinlock for event set access locking */
  std::string fProfilingFile;      /** gperftools profiling file name */

  // State data
  std::atomic_int fPriorityEvents;       /** Number of prioritized events */
  std::atomic_int fTaskId;               /** Counter providing unique task id's */
  BitSet *fDoneEvents = nullptr;         /** Array of bits marking done events */
  std::vector<std::thread> fListThreads; /** Vector of threads */

  bool fDoFastSim = false;               /** Flag to swith on/off the fast simulation stage */

public:
  RunManager() {}
  RunManager(unsigned int npropagators, unsigned int nthreads, GeantConfig *config);
  ~RunManager();

  // Accessors
  void AddEventSet(EventSet *workload);

  GEANT_FORCE_INLINE
  int GetNthreads() const { return fNthreads; }

  GEANT_FORCE_INLINE
  int GetNthreadsTotal() const { return (fNthreads * fNpropagators); }

  GEANT_FORCE_INLINE
  int GetNpropagators() const { return fNpropagators; }

  GEANT_FORCE_INLINE
  GeantConfig *GetConfig() const { return fConfig; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNvolumes() const { return fNvolumes; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNvolumes(int nvol) { fNvolumes = nvol; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  vector_t<Volume_t const *> &GetVolumes() { return fVolumes; }

  GEANT_FORCE_INLINE
  int GetNprimaries() const { return (fEventServer) ? fEventServer->GetNprimaries() : 0; }

  GEANT_FORCE_INLINE
  int GetInitialShare() const { return fInitialShare; }

  GEANT_FORCE_INLINE
  UserDetectorConstruction *GetDetectorConstruction() const { return fDetConstruction; }

  GEANT_FORCE_INLINE
  void SetInitialShare(int nbaskets) { fInitialShare = nbaskets; }

  GEANT_FORCE_INLINE
  std::atomic_int &GetPriorityEvents() { return fPriorityEvents; }

  GEANT_FORCE_INLINE
  int GetNpriority() const { return fPriorityEvents.load(); }

  GEANT_FORCE_INLINE
  Volume_t const *GetVolume(int ivol) const { return fVolumes[ivol]; }

  GEANT_FORCE_INLINE
  Event *GetEvent(int slot) const { return fEventServer->GetEvent(slot); }

  GEANT_FORCE_INLINE
  EventServer *GetEventServer() const { return fEventServer; }

  //  GEANT_FORCE_INLINE
  //  TaskData *GetTaskData(int tid) { return fTaskData[tid]; }

  GEANT_FORCE_INLINE
  int GetTaskId() { return (fTaskId.fetch_add(1)); }

  //  GEANT_FORCE_INLINE
  //  int GetNtracks(int islot) const { return fNtracks[islot]; }

  GEANT_FORCE_INLINE
  void SetCoprocessorBroker(TaskBroker *broker) { fBroker = broker; }

  GEANT_FORCE_INLINE
  void SetUserApplication(UserApplication *app) { fApplication = app; }

  /** @brief Set object to initialize detector, field */
  void SetUserFieldConstruction(UserFieldConstruction *udc);

  GEANT_FORCE_INLINE
  UserApplication *GetUserApplication() const { return fApplication; }

  GEANT_FORCE_INLINE
  void SetDetectorConstruction(UserDetectorConstruction *det) { fDetConstruction = det; }

  GEANT_FORCE_INLINE
  void SetPhysicsInterface(PhysicsInterface *interface) { fPhysicsInterface = interface; }

  GEANT_FORCE_INLINE
  PhysicsInterface *GetPhysicsInterface() const { return fPhysicsInterface; }

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
  void SetProfilingFile(const char *fname) { fProfilingFile = fname; }

  GEANT_FORCE_INLINE
  bool IsInitialized() const { return fInitialized; }

  TaskData *BookTransportTask();

  GEANT_FORCE_INLINE
    void SetFastSim( bool dofastsim ) { fDoFastSim = dofastsim; }

  /** @brief Function checking if transport is completed */
  bool TransportCompleted() const { return ((int)fDoneEvents->FirstNullBit() >= fConfig->fNtotal); }

  /** @brief Check if event is finished */
  bool IsDoneEvent(int ievt) const { return fDoneEvents->TestBitNumber(ievt); }

  /** @brief Mark an event as finished */
  void SetDoneEvent(int ievt) { fDoneEvents->SetBitNumber(ievt); }

  /** @brief Initialize classes for RK Integration */
  void PrepareRkIntegration();

  void LockEventSets()
  {
    while (fEventSetsLock.test_and_set(std::memory_order_acquire)) {
    };
  }
  void UnlockEventSets() { fEventSetsLock.clear(std::memory_order_release); }

  EventSet *NotifyEventSets(Event *finished_event);

  void EventTransported(Event *event, TaskData *td);
  bool Initialize();
  void InitializeRKdata(TaskData *td) const;
  bool FinishRun();
  bool LoadGeometry(const char *filename);
  void RunSimulation();
  bool RunSimulationTask(EventSet *workload, TaskData *td);
  void StopTransport();
};

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

#endif // GEANT_RUN_MANAGER_H
