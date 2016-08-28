#ifndef GEANT_RUN_MANAGER_H
#define GEANT_RUN_MANAGER_H

#include "base/Vector.h"
#include "Geant/Typedefs.h"

class GeantConfig;
class GeantPropagator;
class TaskBroker;
class PhysicsProcess;
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
  PhysicsProcess *fProcess = nullptr;           /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcess *fVectorPhysicsProcess = nullptr; /** Interface to vector physics final state sampling */
  PrimaryGenerator *fPrimaryGenerator = nullptr;   /** Primary generator */
  MCTruthMgr *fTruthMgr = nullptr; /** MCTruth manager */
   
  vector_t<GeantPropagator *> fPropagators;
  vector_t<Volume_t const *> fVolumes;

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
  Volume_t const *GetVolume(int ivol) { return fVolumes[ivol]; }

  GEANT_FORCE_INLINE
  void SetCoprocessorBroker(TaskBroker *broker) { fBroker = broker; }

  GEANT_FORCE_INLINE
  void SetUserApplication(GeantVApplication *app) { fApplication = app; }

  GEANT_FORCE_INLINE
  void SetTaskMgr(GeantVTaskMgr *taskmgr) { fTaskMgr = taskmgr; }

  GEANT_FORCE_INLINE
  void SetPhysicsProcess(PhysicsProcess *proc) { fProcess = proc; }

  GEANT_FORCE_INLINE
  void SetVectorPhysicsProcess(PhysicsProcess *proc) { fVectorPhysicsProcess = proc; }

  GEANT_FORCE_INLINE
  void SetPrimaryGenerator(PrimaryGenerator *gen) { fPrimaryGenerator = gen; }

  GEANT_FORCE_INLINE
  void SetMCTruthMgr(MCTruthMgr *mcmgr) { fTruthMgr = mcmgr; } 

  bool Initialize();
  bool FinishRun();
  bool LoadGeometry(const char *filename);
  void RunSimulation();

};

#endif // GEANT_RUN_MANAGER_H
