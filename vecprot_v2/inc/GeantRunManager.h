#ifndef GEANT_RUN_MANAGER_H
#define GEANT_RUN_MANAGER_H

#include <atomic>
#include "base/Vector.h"
#include "base/BitSet.h"
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"
#include "GeantConfig.h"

using namespace Geant;
using namespace veccore;

class GeantPropagator;
class GeantEvent;
class TaskBroker;
class PhysicsProcessOld;
class PhysicsInterface;
class GeantVApplication;
class GeantVTaskMgr;
class PrimaryGenerator;
class MCTruthMgr;

// Volume-basket manager connector structure attached to volumes as extension
struct VBconnector {
  int index;                      /** Index of basket manager */
  VBconnector(int i) : index(i) {}
};

class GeantRunManager
{
public:
  template <class T>
  using vector_t = vecgeom::Vector<T>;

private:
  bool fInitialized = false;
  int fNpropagators = 0;          /** Number of propagators */
  int fNthreads     = 0;          /** Number of threads per propagator */
  int fNvolumes     = 0;          /** Number of active volumes in the geometry */
  int fNprimaries   = 0;          /** Total number of primaries in the run */
  GeantConfig *fConfig = nullptr; /** Run configuration */
  TaskBroker *fBroker = nullptr;  /** Task broker */

  GeantVApplication *fApplication = nullptr;    /** User application */
  GeantVApplication *fStdApplication = nullptr; /** Standard application */
  GeantVTaskMgr     *fTaskMgr = nullptr;        /** GeantV task manager */
  PhysicsProcessOld *fProcess = nullptr;           /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcessOld *fVectorPhysicsProcess = nullptr; /** Interface to vector physics final state sampling */
  PhysicsInterface *fPhysicsInterface; /** The new, real physics interface */
  PrimaryGenerator *fPrimaryGenerator = nullptr;   /** Primary generator */
  MCTruthMgr *fTruthMgr = nullptr; /** MCTruth manager */
   
  vector_t<GeantPropagator *> fPropagators;
  vector_t<Volume_t const *> fVolumes;

  // State data
  std::atomic_int fPriorityEvents; /** Number of prioritized events */
  std::atomic_flag fFeederLock = ATOMIC_FLAG_INIT; /** Atomic flag to protect the particle feeder */
  BitSet *fDoneEvents = nullptr;   /** Array of bits marking done events */
  int *fNtracks = nullptr;         /** ![fNbuff] Number of tracks per slot */
  GeantEvent **fEvents = nullptr;  /** ![fNbuff]    Array of events */
  GeantTaskData **fTaskData = nullptr; /** ![fNthreads] Data private to threads */
  
private:
  bool LoadVecGeomGeometry();
  void InitNavigators();

public:
  GeantRunManager() {}
  GeantRunManager(unsigned int npropagators, unsigned int nthreads, GeantConfig *config);
  ~GeantRunManager();

// Accessors
  GEANT_FORCE_INLINE
  int  GetNthreads() { return fNthreads; }

  GEANT_FORCE_INLINE
  int  GetNthreadsTotal() { return (fNthreads*fNpropagators); }

  GEANT_FORCE_INLINE
  GeantConfig *GetConfig() { return fConfig; }

  GEANT_FORCE_INLINE
  int  GetNvolumes() { return fNvolumes; }
  
  GEANT_FORCE_INLINE
  vector_t<Volume_t const *> &GetVolumes() { return fVolumes; }

  GEANT_FORCE_INLINE
  int  GetNprimaries() { return fNprimaries; }

  GEANT_FORCE_INLINE
  std::atomic_int &GetPriorityEvents() { return fPriorityEvents; }

  GEANT_FORCE_INLINE
  int GetNpriority() const { return fPriorityEvents.load(); }
  
  GEANT_FORCE_INLINE
  Volume_t const *GetVolume(int ivol) { return fVolumes[ivol]; }

  GEANT_FORCE_INLINE
  GeantEvent *GetEvent(int i) { return fEvents[i]; }

  GEANT_FORCE_INLINE
  GeantTaskData *GetTaskData(int tid) { return fTaskData[tid]; }

  GEANT_FORCE_INLINE
  int GetNtracks(int islot) { return fNtracks[islot]; }

  GEANT_FORCE_INLINE
  void SetCoprocessorBroker(TaskBroker *broker) { fBroker = broker; }

  GEANT_FORCE_INLINE
  void SetUserApplication(GeantVApplication *app) { fApplication = app; }

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
  void SetPrimaryGenerator(PrimaryGenerator *gen) { fPrimaryGenerator = gen; }

  GEANT_FORCE_INLINE
  void SetMCTruthMgr(MCTruthMgr *mcmgr) { fTruthMgr = mcmgr; } 

  /** @brief Function checking if transport is completed */
  bool TransportCompleted() const { return ((int)fDoneEvents->FirstNullBit() >= fConfig->fNtotal); }

  /** @brief Feeder for importing tracks */
  int Feeder(GeantTaskData *td);

  /** @brief Check if transport is feeding with new tracks. */
  GEANT_FORCE_INLINE
  bool IsFeeding() {
    bool feeding = fFeederLock.test_and_set(std::memory_order_acquire);
    if (feeding)
      return true;
    fFeederLock.clear(std::memory_order_release);
    return false;
  }

  /** @brief Try to acquire the lock */
  bool TryLock() { return (fFeederLock.test_and_set(std::memory_order_acquire)); }

  /** @brief Release the lock */
  void ReleaseLock() { fFeederLock.clear(std::memory_order_release); }

  /** @brief Function for importing tracks */
  int ImportTracks(int nevents, int startevent, int startslot, GeantTaskData *td);

  /** @brief Initialize classes for RK Integration */
  void PrepareRkIntegration();

  bool Initialize();
  bool FinishRun();
  bool LoadGeometry(const char *filename);
  void RunSimulation();

};

#endif // GEANT_RUN_MANAGER_H
