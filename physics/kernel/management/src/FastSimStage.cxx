
#include "Geant/FastSimStage.h"

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
#include "Geant/FastSimHandler.h"

namespace geantphysics {

FastSimStage::FastSimStage(geant::Propagator *prop) : geant::SimulationStage(geant::kFastSimStage, prop) {}

// base class will delete the created handlers
FastSimStage::~FastSimStage() {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FastSimStage::FastSimStage(const FastSimStage &other) : geant::SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FastSimStage &FastSimStage::operator=(const FastSimStage &other)
{
  geant::SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
geant::SimulationStage *FastSimStage::Clone() const
{
  FastSimStage *stage = new FastSimStage(*this);
  return stage;
}

int FastSimStage::CreateHandlers()
{
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new FastSimHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have any fast sim processes
geant::Handler *FastSimStage::Select(geant::Track *track, geant::TaskData * /*td*/)
{
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int particleCode         = track->GVcode();
  const Particle *particle = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());

  if (!pManager || (pManager->GetListFastSimProcesses().size() == 0)) {
    // no fast sim process, nothing to do => no handler action
    return nullptr;
  }

  for (unsigned long i = 0; i < pManager->GetListFastSimProcesses().size(); ++i) {

    // give back the only one handler that will call the fast sim process
    if (pManager->GetListFastSimProcesses()[i]->IsApplicable(track)) {
      track->SetPhysicsProcessIndex(i);
      return fHandlers[0];
    }
  }

  return nullptr;
}

} // namespace geantphysics
