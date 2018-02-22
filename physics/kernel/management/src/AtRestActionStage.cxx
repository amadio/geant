
#include "AtRestActionStage.h"

// from geantV
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "Handler.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"
#include "Particle.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"

// handler(s)
#include "AtRestActionHandler.h"


namespace geantphysics {

AtRestActionStage::AtRestActionStage(geant::GeantPropagator *prop)
: SimulationStage(geant::kAtRestActionStage, prop) {}


// base class will delete the created handlers
AtRestActionStage::~AtRestActionStage() {}


int AtRestActionStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new AtRestActionHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}


// Selects tracks that have any at-rest processes
geant::Handler* AtRestActionStage::Select(geant::GeantTrack *track, geant::GeantTaskData * /*td*/) {
  if (track->T()<=0.) {
    // here we will get the MaterialCuts from the LogicalVolume
    const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
    // get the internal code of the particle
    int   particleCode         = track->GVcode();
    const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
    // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager || (pManager->GetListAtRestCandidateProcesses().size()==0)) {
      // no any physics process for the particle or no any at-rest process
      return nullptr;
    }
    // at-rest process needs to be invoked
    return fHandlers[0];
  }
  return nullptr;
}



}  // namespace geantphysics
