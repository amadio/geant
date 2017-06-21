
#include "PrePropagationStage.h"

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

#include "PrePropagationHandler.h"


namespace geantphysics {

PrePropagationStage::PrePropagationStage(Geant::GeantPropagator *prop)
: SimulationStage(Geant::kPrePropagationStage, prop) { }

// base class will delete the created handlers
PrePropagationStage::~PrePropagationStage() {}


int PrePropagationStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new PrePropagationHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have msc process
Geant::Handler* PrePropagationStage::Select(Geant::GeantTrack *track, Geant::GeantTaskData * /*td*/) {
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  if (!pManager || !(pManager->HasMSCProcess())) {
    // no any physics process for the particle or no msc process
    return nullptr;
  }
  // give back the only one handler that will perform the msc step limit phase actions
  return fHandlers[0];
}


}  // namespace geantphysics
