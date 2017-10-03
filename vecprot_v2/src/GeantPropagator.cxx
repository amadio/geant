// A simple propagator taking as input a set of particles located in a given
// volume AND the global matrix of the volume.
// The ProposeStep() method choses between a "scattering" process with no eloss
// and a "ionization" process and generates a random "physical" step. In this simple
// model all particlea undertake the same list of processes
// The ScatteringProcess() method emulates scattering and changes the particle
// direction randomly in a forward cone with opening angle proportional with 1/p
// The IonizationProcess() method simulates an energy deposition with an amount
// epsil*int(1+K*rnd) (epsil, 2*epsil, ..., K*epsil)
// In future we can add absorption and decay to check the stack model...
// The PropagateInField(step) method propagates in a uniform magnetic field with
// an amount equal to step
// The Transport() method computes the safety and snext and compares with the
// physical step. If (safety>pstep), PropagateInField(pstep) is called and the
// physics process is simulated. Otherwise, PropagateInField(safety) is called
// and safety subtracted from the physical step, then the procedure is repeated
// until C*snext/4 < 1E-6 (tangent of angle with sagita, C=1/R is the curvature)
//
#include "GeantPropagator.h"

#ifdef USE_ROOT
#include "TTimer.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include <fenv.h>
#endif


#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/VNavigator.h"
#ifdef USE_ROOT
#include "management/RootGeoManager.h"
#endif
#include "volumes/PlacedVolume.h"
#else
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#endif

#include "Geant/Error.h"
#include "LocalityManager.h"
#include "TrackManager.h"
#include "GeantRunManager.h"
#include "GeantTrackVec.h"
#include "PhysicsInterface.h"
#include "PhysicsProcessOld.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "GeantVTaskMgr.h"
#include "StdApplication.h"
#include "GeantFactoryStore.h"
#include "GeantEvent.h"
#include "GeantScheduler.h"
#include "PrimaryGenerator.h"
#include "MCTruthMgr.h"

#include "PreStepStage.h"
#include "XSecSamplingStage.h"
#include "GeomQueryStage.h"
#include "PropagationStage.h"
#include "ContinuousProcStage.h"
#include "DiscreteProcStage.h"
#include "SteppingActionsStage.h"

#ifdef USE_CALLGRIND_CONTROL
#include <valgrind/callgrind.h>
#endif

using namespace vecgeom;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeantPropagator::GeantPropagator(int nthreads)
    : fNthreads(nthreads), fNtransported(0), fNsteps(0), fNsnext(0),
      fNphys(0), fNmag(0), fNsmall(0), fNcross(0), fNpushed(0), fNkilled(0), fNidle(0), fNbfeed(0) {
  // Constructor
  // Single instance of the propagator

  // Initialize workload manager
#ifndef VECCORE_CUDA
  fWMgr = WorkloadManager::NewInstance(this, nthreads);
#endif
}

//______________________________________________________________________________
GeantPropagator::GeantPropagator(const GeantPropagator &orig) : GeantPropagator(orig.fNthreads) {

  fRunMgr = orig.fRunMgr;
  SetConfig(orig.fConfig);
  fApplication = orig.fApplication;
  fStdApplication = orig.fStdApplication;
  fTaskMgr = orig.fTaskMgr;
  fProcess = orig.fProcess;
  fPhysicsInterface = orig.fPhysicsInterface;
  fVectorPhysicsProcess = orig.fVectorPhysicsProcess;
  fPrimaryGenerator = orig.fPrimaryGenerator;
  fTruthMgr = orig.fTruthMgr;

}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator() {
  // Destructor
  delete fTimer;
  delete fWMgr;
}

//______________________________________________________________________________
int GeantPropagator::AddTrack(GeantTrack &track) {
#ifdef VECCORE_CUDA
  assert(0 && "DispatchTrack not implemented yet for CUDA host/device code.");
  return 0;
#else
    // Add a new track in the system. returns track number within the event.
  track.fParticle = fRunMgr->GetEvent(track.fEvslot)->AddTrack();

  // call MCTruth manager if it has been instantiated
  if(fTruthMgr) fTruthMgr->AddTrack(track);

  fNtransported++;
  return track.fParticle;
#endif
}

//______________________________________________________________________________
int GeantPropagator::DispatchTrack(GeantTrack &track, GeantTaskData *td) {
  // Dispatch a registered track produced by the generator.
#ifdef VECCORE_CUDA
  assert(0 && "DispatchTrack not implemented yet for CUDA host/device code.");
  return 0;
#else  // stoping track in MCTruthManager
  return fWMgr->GetScheduler()->AddTrack(track, td);
#endif
}

//______________________________________________________________________________
int GeantPropagator::GetNpending() const {
  // Returns number of baskets pending in the queue
  return fWMgr->GetNpending();
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(const GeantTrack_v &tracks, int itr, GeantTaskData *td) {
  // Mark track as stopped for tracking.
  //   Printf("Stopping track %d", track->particle);

#ifdef VECCORE_CUDA
  assert(0 && "StopTrack not implemented yet for CUDA host/device code.");
#else  // stoping track in MCTruthManager
  if(fTruthMgr)
    {
      if(tracks.fStatusV[itr] == kKilled) fTruthMgr->EndTrack(tracks, itr);
    }

  if (fRunMgr->GetEvent(tracks.fEvslotV[itr])->StopTrack(fRunMgr, td)) {
    std::atomic_int &priority_events = fRunMgr->GetPriorityEvents();
    priority_events++;
  }
#endif
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(GeantTrack *track, GeantTaskData *td) {
  // Mark track as stopped for tracking.
  //   Printf("Stopping track %d", track->particle);

#ifdef VECCORE_CUDA
  assert(0 && "StopTrack not implemented yet for CUDA host/device code.");
#else  // stoping track in MCTruthManager
  // stoping track in MCTruthManager
  if(fTruthMgr)
    {
      if(track->fStatus == kKilled) fTruthMgr->EndTrack(track);
    }

  if (fRunMgr->GetEvent(track->fEvslot)->StopTrack(fRunMgr, td)) {
    std::atomic_int &priority_events = fRunMgr->GetPriorityEvents();
    priority_events++;
  }
#endif
}

//______________________________________________________________________________
bool GeantPropagator::IsIdle() const {
  // Check if work queue is empty and all used threads are waiting
  return (!fCompleted && GetNworking()==0 && GetNpending()==0);
}

//______________________________________________________________________________
void GeantPropagator::Initialize() {
  // Initialize the propagator.
#ifndef VECCORE_CUDA
  LocalityManager *mgr = LocalityManager::Instance();
  int numa = (fNuma >= 0) ? fNuma : 0;
  fTrackMgr = &mgr->GetTrackManager(numa);
#endif

  // Add some empty baskets in the queue
#ifdef VECCORE_CUDA
  // assert(0 && "Initialize not implemented yet for CUDA host/device code.");
#else
  if (!fConfig->fUseV3)
    fWMgr->CreateBaskets(this);
#endif
  CreateSimulationStages();
}

//______________________________________________________________________________
void GeantPropagator::SetNuma(int numa)
{
// Set locality for the propagator
#ifndef VECCORE_CUDA
  LocalityManager *mgr = LocalityManager::Instance();
  fNuma  = numa;
  if (!mgr->IsInitialized()) return;
  if (numa >= 0)
    fTrackMgr = &mgr->GetTrackManager(numa);
  else
    fTrackMgr = &mgr->GetTrackManager(0);
#endif
}

//NOTE: We don't do anything here so it's not called from the WorkloadManager anymore
//______________________________________________________________________________
void GeantPropagator::ApplyMsc(int /*ntracks*/, GeantTrack_v & /*tracks*/, GeantTaskData * /*td*/) {
/*
  // Apply multiple scattering for charged particles.
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = (Material_t *)td->fVolume->GetMaterialPtr();
#else
    mat = td->fVolume->GetMaterial();
#endif

#ifdef USE_REAL_PHYSICS
  tracks;
#else
  // actually nothing happens in this call in the TTabPhysProcess
  fProcess->ApplyMsc(mat, ntracks, tracks, td);
#endif
*/
}

//______________________________________________________________________________
void GeantPropagator::ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Generate all physics steps for the tracks in trackin.
  // Reset the current step length to 0
  for (int i = 0; i < ntracks; ++i) {
    tracks.fStepV[i] = 0.;
    tracks.fEdepV[i] = 0.;
  }
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = (Material_t *)td->fVolume->GetMaterialPtr();
#else
    mat = td->fVolume->GetMaterial();
#endif

#ifdef USE_REAL_PHYSICS
  fPhysicsInterface->ComputeIntLen(mat, ntracks, tracks, 0, td);
#else
  fProcess->ComputeIntLen(mat, ntracks, tracks, td);
#endif
}

//______________________________________________________________________________
void GeantPropagator::RunSimulation(GeantPropagator *prop, int nthreads)
{
// Static thread running the simulation for one propagator
  Info("RunSimulation", "Starting propagator %p with %d threads", prop, nthreads);
  prop->PropagatorGeom(nthreads);
}

//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(int nthreads) {
#ifdef VECCORE_CUDA
  assert(0 && "PropagatorGeom not implemented yet for CUDA host/device code.");
#else
// Steering propagation method
  fNthreads = nthreads;

  // Import the input events. This will start also populating the main queue

  // Loop baskets and transport particles until there is nothing to transport anymore
  fTransportOngoing = true;

#ifdef USE_CALLGRIND_CONTROL
  CALLGRIND_START_INSTRUMENTATION;
#endif

  // Start system tasks
  if (!fWMgr->StartTasks(fTaskMgr)) {
    Fatal("PropagatorGeom", "%s", "Cannot start tasks.");
    return;
  }

  // Wait all workers to finish, then join all threads
  fWMgr->WaitWorkers();
  fCompleted = true; // Make sure someone is not just feeding work...
  fWMgr->JoinThreads();
#ifdef USE_CALLGRIND_CONTROL
  CALLGRIND_STOP_INSTRUMENTATION;
  CALLGRIND_DUMP_STATS;
#endif
#endif
}

//______________________________________________________________________________
void GeantPropagator::StopTransport()
{
// Stop the transport threads. Needed only when controlling the transport
// from the transport manager
#ifdef VECCORE_CUDA
  assert(0 && "StopTransport not implemented yet for CUDA host/device code.");
#else
  if (fCompleted) return;
  std::unique_lock<std::mutex> lk(fStopperLock);
  fCompleted = true;
  Printf("+++ Stopping propagator %p", this);
  fWMgr->StopTransportThreads();
#endif
}

//______________________________________________________________________________
void GeantPropagator::SetTaskBroker(TaskBroker *broker) {
  // Setter for task broker
#ifdef VECCORE_CUDA
  assert(0 && "SetTaskBroker not implemented yet for CUDA host/device code.");
#else
  fWMgr->SetTaskBroker(broker);
#endif
}

//______________________________________________________________________________
TaskBroker *GeantPropagator::GetTaskBroker() {
  // Getter for task broker
#ifdef VECCORE_CUDA
  assert(0 && "GetTaskBroker not implemented yet for CUDA host/device code.");
  return nullptr;
#else
  return fWMgr->GetTaskBroker();
#endif
}

//______________________________________________________________________________
void GeantPropagator::SetConfig(GeantConfig *config)
{
// Set run configuration.
  fConfig = config;
  fNbuff = config->fNbuff;
  fNtotal = config->fNtotal;
  if (fNtotal <= 0 || fNbuff <= 0) {
    Fatal("GeantPropagator::SetConfig", "%s", "Number of transported/buffered events should be positive");
    return;
  }
  if (fNbuff > fNtotal) {
    Info("GeantPropagator::SetCofig", "%s", "Number of buffered events changed to %d", fNtotal);
    fNbuff = fNtotal;
  }
  // Instantiate factory store
  GeantFactoryStore::Instance(fNbuff);
}

//______________________________________________________________________________
int GeantPropagator::ShareWork(GeantPropagator &other)
{
#ifdef VECCORE_CUDA
  assert(0 && "ShareWork not implemented yet for CUDA host/device code.");
  return 0;
#else
  return ( fWMgr->ShareBaskets(other.fWMgr) );
#endif
}

#ifdef USE_REAL_PHYSICS
//______________________________________________________________________________
int GeantPropagator::CreateSimulationStages()
{
  // Create stages in the same order as the enumeration ESimulationStage
  SimulationStage *stage = nullptr;
  (void)stage;
  // kPreStepStage
  stage = new PreStepStage(this);
  if (stage->GetId() != int(kPreStepStage))
    Fatal("CreateSimulationStages", "Wrong stages start index");
  // kComputeIntLStage: physics step limit
  stage = fPhysicsInterface->CreateComputeIntLStage(this);
  assert(stage->GetId() == int(kComputeIntLStage));
  // kGeometryStepStage
  stage = new GeomQueryStage(this);
  assert(stage->GetId() == int(kGeometryStepStage));
  // kPrePropagationStage: possible msc step limit or msc
  stage = fPhysicsInterface->CreatePrePropagationStage(this);
  assert(stage->GetId() == int(kPrePropagationStage));

//  stage = new GeomQueryStage(this);
//  assert(stage->GetId() == int(kGeometryStepStage));


  // kPropagationStage
  stage = new PropagationStage(this);
  assert(stage->GetId() == int(kPropagationStage));
  // kPostPropagationStage
  stage = fPhysicsInterface->CreatePostPropagationStage(this);
  assert(stage->GetId() == int(kPostPropagationStage));

  // kContinuousProcStage
  stage = fPhysicsInterface->CreateAlongStepActionStage(this);
  assert(stage->GetId() == int(kAlongStepActionStage));
  // kDiscreteProcStage
  stage = fPhysicsInterface->CreatePostStepActionStage(this);
  assert(stage->GetId() == int(kPostStepActionStage));
  // kSteppingActionsStage
  stage = new SteppingActionsStage(this);
  assert(stage->GetId() == int(kSteppingActionsStage));

  /**************************************
   *  Define connections between stages *
   **************************************/
  GetStage(kPreStepStage)->SetFollowUpStage(kComputeIntLStage, false);
  // Follow-up not unique: new tracks may be killed by the user -> SteppingActions
  GetStage(kPreStepStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kComputeIntLStage)->SetFollowUpStage(kGeometryStepStage, true);
//  GetStage(kComputeIntLStage)->SetFollowUpStage(kPrePropagationStage, true);
  GetStage(kComputeIntLStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kGeometryStepStage)->SetFollowUpStage(kPrePropagationStage, true);
  GetStage(kGeometryStepStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kPrePropagationStage)->SetFollowUpStage(kPropagationStage, true);
//  GetStage(kPrePropagationStage)->SetFollowUpStage(kGeometryStepStage, true);
  GetStage(kPrePropagationStage)->ActivateBasketizing(false);
  //        V

//  GetStage(kGeometryStepStage)->SetFollowUpStage(kPropagationStage, true);
//  GetStage(kGeometryStepStage)->ActivateBasketizing(false);

  //        V
  //        V
  GetStage(kPropagationStage)->SetFollowUpStage(kPostPropagationStage, false);
  // Follow-up not unique: stuck tracks are killed -> SteppingActions
  GetStage(kPropagationStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kPostPropagationStage)->SetFollowUpStage(kAlongStepActionStage, true);
  GetStage(kPrePropagationStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kAlongStepActionStage)->SetFollowUpStage(kPostStepActionStage, false);
  // Follow-up not unique: particle can be become at-rest/killed and new particles can be produced
  //   - (at the moment we don't have at-rest stage because we don't have at rest process)
  GetStage(kAlongStepActionStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kPostStepActionStage)->SetFollowUpStage(kSteppingActionsStage, false);
  // Follow-up not unique: particle can be become at-rest/killed and new particles can be produced
  //   - (at the moment we don't have at-rest stage because we don't have at rest process)
  GetStage(kPostStepActionStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kSteppingActionsStage)->SetEndStage();
  GetStage(kSteppingActionsStage)->ActivateBasketizing(false);

  for (auto stage : fStages) {
    int nhandlers = stage->CreateHandlers();
    (void)nhandlers;
    assert((nhandlers > 0) && "Number of handlers for a simulation stage cannot be 0");
  }
  return fStages.size();
}
#else
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantPropagator::CreateSimulationStages()
{
  // Create stages in the same order as the enumeration ESimulationStage
  SimulationStage *stage = nullptr;
  (void)stage;
  // kPreStepStage
  stage = new PreStepStage(this);
  if (stage->GetId() != int(kPreStepStage))
    Fatal("CreateSimulationStages", "Wrong stages start index");
  // kXSecSamplingStage
  stage = new XSecSamplingStage(this);
  assert(stage->GetId() == int(kXSecSamplingStage));
  // kGeometryStepStage
  stage = new GeomQueryStage(this);
  assert(stage->GetId() == int(kGeometryStepStage));
  // kPropagationStage
  stage = new PropagationStage(this);
  assert(stage->GetId() == int(kPropagationStage));
  // kMSCStage
  // stage = new MSCStage(this);
  //assert(stage->GetId() == int(kMSCStage));
  // kContinuousProcStage
  stage = new ContinuousProcStage(this);
  assert(stage->GetId() == int(kContinuousProcStage));
  // kDiscreteProcStage
  stage = new DiscreteProcStage(this);
  assert(stage->GetId() == int(kDiscreteProcStage));
  // kSteppingActionsStage
  stage = new SteppingActionsStage(this);
  assert(stage->GetId() == int(kSteppingActionsStage));

  /**************************************
   *  Define connections between stages *
   **************************************/
  GetStage(kPreStepStage)->SetFollowUpStage(kXSecSamplingStage, false);
  // Follow-up not unique: new tracks may be killed by the user -> SteppingActions
  GetStage(kPreStepStage)->ActivateBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kXSecSamplingStage)->SetFollowUpStage(kGeometryStepStage, true);
  GetStage(kXSecSamplingStage)->ActivateBasketizing(true);
  //        V
  //        V
  //        V
  GetStage(kGeometryStepStage)->SetFollowUpStage(kPropagationStage, true);
  GetStage(kGeometryStepStage)->ActivateBasketizing(true);
  //        V
  //        V
  //        V
  GetStage(kPropagationStage)->SetFollowUpStage(kContinuousProcStage, false);
  // Follow-up not unique: stuck tracks are killed -> SteppingActions
  GetStage(kPropagationStage)->ActivateBasketizing(true);
  //        V
  //        V
  //        V
  GetStage(kContinuousProcStage)->SetFollowUpStage(kDiscreteProcStage, true);
  GetStage(kContinuousProcStage)->ActivateBasketizing(true);
  //        V
  //        V
  //        V
  GetStage(kDiscreteProcStage)->SetFollowUpStage(kSteppingActionsStage, true);
  GetStage(kDiscreteProcStage)->ActivateBasketizing(true);
  //        V
  //        V
  //        V
  GetStage(kSteppingActionsStage)->SetEndStage();
  GetStage(kSteppingActionsStage)->ActivateBasketizing(false);

  for (auto stage : fStages) {
    int nhandlers = stage->CreateHandlers();
    (void)nhandlers;
    assert((nhandlers > 0) && "Number of handlers for a simulation stage cannot be 0");
  }
  return fStages.size();
}
#endif

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantPropagator::GetNextStage(GeantTrack &/*track*/, int /*current*/)
{
// Get the next simulation stage for a track
//  0 - Sample X-sec for discrete processes to propose the physics step
//  1 - Compute geometry step length to boundary with a physics limitation
//  2 - Propagate track with selected step, performing relocation if needed.
//        - Follow-up to 1 for charged tracks if neither the geometry nor
//          physics steps are completed
//        - In case MSC is available, apply it for charged tracks
//        - Detect loopers and send them to RIP stage
//  2'- Apply multiple scattering and change track position/direction
//        - Follow-up to 1 after
//        - Handle boundary crossing by MSC
//  3 - Apply along-step continuous processes
//        - Send particles killed by energy threshold to graveyard
//  4 - Do post step actions for particles suffering a physics process
//        - Follow-up to stage 0 after running user actions stage
//  5 - RIP stage - execute user actions then terminate tracks
//  6 - Stepping actions is invoked as a stage, but the follow-up stage is backed-up
//      beforehand in the track state.
  return -1;
}


} // GEANT_IMPL_NAMESPACE
} // Geant
