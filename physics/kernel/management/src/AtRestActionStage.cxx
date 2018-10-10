
#include "Geant/AtRestActionStage.h"

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
#include "Geant/AtRestActionHandler.h"

namespace geantphysics {

AtRestActionStage::AtRestActionStage(geant::Propagator *prop) : geant::SimulationStage(geant::kAtRestActionStage, prop)
{
}

// base class will delete the created handlers
AtRestActionStage::~AtRestActionStage() {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
AtRestActionStage::AtRestActionStage(const AtRestActionStage &other) : geant::SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
AtRestActionStage &AtRestActionStage::operator=(const AtRestActionStage &other)
{
  geant::SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
geant::SimulationStage *AtRestActionStage::Clone() const
{
  AtRestActionStage *stage = new AtRestActionStage(*this);
  return stage;
}

int AtRestActionStage::CreateHandlers()
{
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new AtRestActionHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have any at-rest processes
geant::Handler *AtRestActionStage::Select(geant::Track *track, geant::TaskData * /*td*/)
{
  if (track->Ekin() <= 0.) {
    // here we will get the MaterialCuts from the LogicalVolume
    const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
        (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
    // get the internal code of the particle
    int particleCode         = track->GVcode();
    const Particle *particle = Particle::GetParticleByInternalCode(particleCode);
    // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager || (pManager->GetListAtRestCandidateProcesses().size() == 0)) {
      // no any physics process for the particle or no any at-rest process
      return nullptr;
    }
    // at-rest process needs to be invoked
    return fHandlers[0];
  }
  return nullptr;
}

} // namespace geantphysics
