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
#include "Geant/Propagator.h"

#ifdef USE_ROOT
#include "TTimer.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include <fenv.h>
#endif

#include "navigation/VNavigator.h"
#ifdef USE_ROOT
#include "management/RootGeoManager.h"
#endif
#include "volumes/PlacedVolume.h"

#include "Geant/Error.h"
#include "Geant/LocalityManager.h"
#include "Geant/TrackManager.h"
#include "Geant/RunManager.h"
#include "Geant/PhysicsInterface.h"
#include "Geant/WorkloadManager.h"
#include "Geant/TaskData.h"
#include "Geant/UserApplication.h"
#include "Geant/StdApplication.h"
#include "Geant/FactoryStore.h"
#include "Geant/Event.h"
#include "Geant/PrimaryGenerator.h"
#include "Geant/MCTruthMgr.h"

#include "Geant/PreStepStage.h"
#include "Geant/GeomQueryStage.h"
#include "Geant/PropagationStage.h"
#include "Geant/SteppingActionsStage.h"

#ifdef USE_CALLGRIND_CONTROL
#include <valgrind/callgrind.h>
#endif

using namespace vecgeom;

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Propagator::Propagator(int nthreads)
    : fNthreads(nthreads), fNtransported(0), fNsteps(0), fNsnext(0), fNphys(0), fNmag(0), fNsmall(0), fNcross(0),
      fNpushed(0), fNkilled(0), fNbfeed(0)
{
// Constructor
// Single instance of the propagator

// Initialize workload manager
#ifndef VECCORE_CUDA
  fWMgr = WorkloadManager::NewInstance(this, nthreads);
#endif
}

//______________________________________________________________________________
Propagator::Propagator(const Propagator &orig) : Propagator(orig.fNthreads)
{

  fRunMgr = orig.fRunMgr;
  SetConfig(orig.fConfig);
  fApplication      = orig.fApplication;
  fStdApplication   = orig.fStdApplication;
  fPhysicsInterface = orig.fPhysicsInterface;
  fPrimaryGenerator = orig.fPrimaryGenerator;
  fTruthMgr         = orig.fTruthMgr;
}

//______________________________________________________________________________
Propagator::~Propagator()
{
  // Destructor
  delete fTimer;
  delete fWMgr;
}

//______________________________________________________________________________
int Propagator::AddTrack(Track &track)
{
#ifdef VECCORE_CUDA
  assert(0 && "DispatchTrack not implemented yet for CUDA host/device code.");
  return 0;
#else
  // Add a new track in the system. returns track number within the event.
  track.SetParticle(fRunMgr->GetEvent(track.EventSlot())->AddTrack());

  // call MCTruth manager if it has been instantiated
  if (fTruthMgr) fTruthMgr->AddTrack(track);

  fNtransported++;
  return track.Particle();
#endif
}

//______________________________________________________________________________
void Propagator::StopTrack(Track *track, TaskData *td)
{
  // Mark track as stopped for tracking.
  //   Printf("Stopping track %d", track->particle);

#ifdef VECCORE_CUDA
  assert(0 && "StopTrack not implemented yet for CUDA host/device code.");
#else // stoping track in MCTruthManager
  // stoping track in MCTruthManager
  if (fTruthMgr) {
    if (track->Status() == kKilled) fTruthMgr->EndTrack(track);
  }

  if (fRunMgr->GetEvent(track->EventSlot())->StopTrack(fRunMgr, td)) {
    std::atomic_int &priority_events = fRunMgr->GetPriorityEvents();
    priority_events++;
  }
#endif
}

//______________________________________________________________________________
void Propagator::Initialize()
{
// Initialize the propagator.
#ifndef VECCORE_CUDA
  LocalityManager *mgr = LocalityManager::Instance();
  int numa             = (fNuma >= 0) ? fNuma : 0;
  fTrackMgr            = &mgr->GetTrackManager(numa);
#endif

  // Create the simulation stages
  CreateSimulationStages();
}

//______________________________________________________________________________
void Propagator::SetNuma(int numa)
{
// Set locality for the propagator
#ifndef VECCORE_CUDA
  LocalityManager *mgr = LocalityManager::Instance();
  fNuma                = numa;
  if (!mgr->IsInitialized()) return;
  if (numa >= 0)
    fTrackMgr = &mgr->GetTrackManager(numa);
  else
    fTrackMgr = &mgr->GetTrackManager(0);
#endif
}

//______________________________________________________________________________
void Propagator::RunSimulation(Propagator *prop, int nthreads)
{
  // Static thread running the simulation for one propagator
  Info("RunSimulation", "Starting propagator %p with %d threads", prop, nthreads);
  prop->PropagatorGeom(nthreads);
}

//______________________________________________________________________________
void Propagator::PropagatorGeom(int nthreads)
{
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
  if (!fWMgr->StartTasks()) {
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
void Propagator::StopTransport()
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
void Propagator::SetTaskBroker(TaskBroker *broker)
{
// Setter for task broker
#ifdef VECCORE_CUDA
  assert(0 && "SetTaskBroker not implemented yet for CUDA host/device code.");
#else
  fWMgr->SetTaskBroker(broker);
#endif
}

//______________________________________________________________________________
TaskBroker *Propagator::GetTaskBroker()
{
// Getter for task broker
#ifdef VECCORE_CUDA
  assert(0 && "GetTaskBroker not implemented yet for CUDA host/device code.");
  return nullptr;
#else
  return fWMgr->GetTaskBroker();
#endif
}

//______________________________________________________________________________
void Propagator::SetConfig(GeantConfig *config)
{
  // Set run configuration.
  fConfig = config;
  fNbuff  = config->fNbuff;
  fNtotal = config->fNtotal;
  if (fNtotal <= 0 || fNbuff <= 0) {
    Fatal("Propagator::SetConfig", "%s", "Number of transported/buffered events should be positive");
    return;
  }
  if (fNbuff > fNtotal) {
    Info("Propagator::SetCofig", "%s", "Number of buffered events changed to %d", fNtotal);
    fNbuff = fNtotal;
  }
  // Instantiate factory store
  FactoryStore::Instance(fNbuff);
}

//______________________________________________________________________________
int Propagator::CreateSimulationStages()
{
  // Create stages in the same order as the enumeration ESimulationStage
  SimulationStage *stage = nullptr;
  (void)stage;
  // kPreStepStage
  stage = new PreStepStage(this);
  if (stage->GetId() != int(kPreStepStage)) Fatal("CreateSimulationStages", "Wrong stages start index");
  // kFastSimStage: optional fast simulation for selected volumes/particles
  stage = fPhysicsInterface->CreateFastSimStage(this);
  assert(stage->GetId() == int(kFastSimStage));
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

  // stage = fPhysicsInterface->CreatePostStepActionStage(this);
  // assert(stage->GetId() == int(kPostStepActionStage));

  stage = fPhysicsInterface->CreatePostStepActionPhysProcessStage(this);
  assert(stage->GetId() == int(kPostStepActionStage));

  // kAtRestActionStage
  stage = fPhysicsInterface->CreateAtRestActionStage(this);
  assert(stage->GetId() == int(kAtRestActionStage));
  // kSteppingActionsStage
  stage = new SteppingActionsStage(this);
  assert(stage->GetId() == int(kSteppingActionsStage));

  // Create the handlers
  for (auto stage : fStages) {
    int nhandlers = stage->CreateHandlers();
    stage->CountBasketized();
    (void)nhandlers;
    assert((nhandlers > 0) && "Number of handlers for a simulation stage cannot be 0");
  }

  /**************************************
   *  Define connections between stages *
   **************************************/
  if (fDoFastSim) {
    GetStage(kPreStepStage)->SetFollowUpStage(kFastSimStage, false);
    // Follow-up not unique: new tracks may be killed by the user -> SteppingActions
    GetStage(kPreStepStage)->SetBasketizing(false);
    //        V
    //        V
    //        V
    GetStage(kFastSimStage)->SetFollowUpStage(kComputeIntLStage, false);
    // Follow-up not unique: new tracks may be killed by the user -> SteppingActions
    GetStage(kFastSimStage)->SetBasketizing(false);
  } else {
    GetStage(kPreStepStage)->SetFollowUpStage(kComputeIntLStage, false);
    // Follow-up not unique: new tracks may be killed by the user -> SteppingActions
    GetStage(kPreStepStage)->SetBasketizing(false);
  }
  //        V
  //        V
  //        V
  GetStage(kComputeIntLStage)->SetFollowUpStage(kGeometryStepStage, true);
  //  GetStage(kComputeIntLStage)->SetFollowUpStage(kPrePropagationStage, true);
  GetStage(kComputeIntLStage)->SetBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kGeometryStepStage)->SetFollowUpStage(kPrePropagationStage, false);
  // Follow-up not unique: tracks may have done pre-propagation/propagation returning to
  // geometry due to partial field propagation
  GetStage(kGeometryStepStage)->SetBasketizing(fConfig->fUseVectorizedGeom);
  //        V
  //        V
  //        V
  GetStage(kPrePropagationStage)->SetFollowUpStage(kPropagationStage, true);
  //  GetStage(kPrePropagationStage)->SetFollowUpStage(kGeometryStepStage, true);
  GetStage(kPrePropagationStage)->SetBasketizing(false);
  //        V

  //  GetStage(kGeometryStepStage)->SetFollowUpStage(kPropagationStage, true);
  //  GetStage(kGeometryStepStage)->SetBasketizing(false);

  //        V
  //        V
  GetStage(kPropagationStage)->SetFollowUpStage(kPostPropagationStage, false);
  // Follow-up not unique: stuck tracks are killed -> SteppingActions
  GetStage(kPropagationStage)->SetBasketizing(fConfig->fUseVectorizedField);
  //        V
  //        V
  //        V
  GetStage(kPostPropagationStage)->SetFollowUpStage(kAlongStepActionStage, true);
  GetStage(kPrePropagationStage)->SetBasketizing(false);
  //        V
  //        V
  //        V
  GetStage(kAlongStepActionStage)->SetFollowUpStage(kPostStepActionStage, false);
  // Follow-up not unique: particle can be become at-rest/killed and new particles can be produced
  //   - (at the moment we don't have at-rest stage because we don't have at rest process)
  GetStage(kAlongStepActionStage)->SetBasketizing(false);
  //        V
  //        V
  //        V
  // GetStage(kPostStepActionStage)->SetFollowUpStage(kSteppingActionsStage, false);
  GetStage(kPostStepActionStage)->SetFollowUpStage(kSteppingActionsStage, false);
  // Follow-up not unique: particle can be become at-rest/killed and new particles can be produced
  //   - (at the moment we don't have at-rest stage because we don't have at rest process)
  GetStage(kPostStepActionStage)->SetBasketizing(fConfig->fUseVectorizedPhysics);
  //        V
  //        V
  //        V
  GetStage(kSteppingActionsStage)->SetEndStage();
  GetStage(kSteppingActionsStage)->SetBasketizing(false);

  return fStages.size();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int Propagator::GetNextStage(Track & /*track*/, int /*current*/) const
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

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
