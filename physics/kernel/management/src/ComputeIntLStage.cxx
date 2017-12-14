
#include "ComputeIntLStage.h"

// from geantV
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "Handler.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"

// handler(s)
#include "ComputeIntLHandler.h"


namespace geantphysics {


ComputeIntLStage::ComputeIntLStage(Geant::GeantPropagator *prop)
: SimulationStage(Geant::kComputeIntLStage, prop) {}


// base class will delete the created handlers
ComputeIntLStage::~ComputeIntLStage() {}


int ComputeIntLStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new ComputeIntLHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}


// Selects tracks that have any physics processes and set some members to others
Geant::Handler* ComputeIntLStage::Select(Geant::GeantTrack *track, Geant::GeantTaskData * /*td*/) {
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  if (!pManager || (pManager->GetListAlongStepProcesses().size()+pManager->GetListPostStepCandidateProcesses().size()==0)) {
    track->SetPstep(PhysicsProcess::GetAVeryLargeValue());
    // no physics limit nothing to do => no handler action
    return nullptr;
  }
  // give back the only one handler that will compute the physics step limit for all particles that has any along/post
  // step interactions
  return fHandlers[0];
}


}  // namespace geantphysics
