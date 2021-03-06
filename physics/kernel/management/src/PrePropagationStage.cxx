
#include "Geant/PrePropagationStage.h"

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

#include "Geant/PrePropagationHandler.h"

namespace geantphysics {

PrePropagationStage::PrePropagationStage(geant::Propagator *prop)
    : geant::SimulationStage(geant::kPrePropagationStage, prop)
{
}

// base class will delete the created handlers
PrePropagationStage::~PrePropagationStage() {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PrePropagationStage::PrePropagationStage(const PrePropagationStage &other) : geant::SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PrePropagationStage &PrePropagationStage::operator=(const PrePropagationStage &other)
{
  geant::SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
geant::SimulationStage *PrePropagationStage::Clone() const
{
  PrePropagationStage *stage = new PrePropagationStage(*this);
  return stage;
}

int PrePropagationStage::CreateHandlers()
{
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new PrePropagationHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have msc process
geant::Handler *PrePropagationStage::Select(geant::Track *track, geant::TaskData * /*td*/)
{
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int particleCode         = track->GVcode();
  const Particle *particle = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  if (!pManager || !(pManager->HasMSCProcess())) {
    // no any physics process for the particle or no msc process
    return nullptr;
  }
  // give back the only one handler that will perform the msc step limit phase actions
  return fHandlers[0];
}

} // namespace geantphysics
