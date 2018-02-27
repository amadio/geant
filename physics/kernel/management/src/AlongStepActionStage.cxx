
#include "Geant/AlongStepActionStage.h"

// from geantV
#include "Geant/Propagator.h"
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/Handler.h"

// from realphysics
#include "Geant/Material.h"
#include "Geant/MaterialCuts.h"
#include "Geant/Region.h"
#include "Geant/Particle.h"

#include "Geant/PhysicsProcess.h"
#include "Geant/PhysicsManagerPerParticle.h"
#include "Geant/LightTrack.h"

// handler(s)
#include "Geant/AlongStepActionHandler.h"


namespace geantphysics {

AlongStepActionStage::AlongStepActionStage(geant::Propagator *prop)
: SimulationStage(geant::kAlongStepActionStage, prop) {}

// base class will delete the created handlers
AlongStepActionStage::~AlongStepActionStage() {}


int AlongStepActionStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new AlongStepActionHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have any along step processes i.e. continuous part
geant::Handler* AlongStepActionStage::Select(geant::Track *track, geant::TaskData * /*td*/) {
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  if (!pManager || (pManager->GetListAlongStepProcesses().size()==0)) {
    // no any physics process for the particle or no any continuous process i.e. no alongstepaction => no handler action
    return nullptr;
  }
  // give back the only one handler for all other particles that will compute the along-step-action(s)
  return fHandlers[0];
}


}  // namespace geantphysics
