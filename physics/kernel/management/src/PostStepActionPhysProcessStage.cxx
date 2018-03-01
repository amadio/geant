
#include "Geant/PhysicsListManager.h"
#include "Geant/Particle.h"
#include "Geant/PostStepActionPhysProcessStage.h"

// from geantV
#include "Geant/Propagator.h"
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/Handler.h"

// from realphysics
#include "Geant/Material.h"
#include "Geant/MaterialCuts.h"
#include "Geant/Region.h"

#include "Geant/PhysicsProcess.h"
#include "Geant/PhysicsManagerPerParticle.h"
#include "Geant/LightTrack.h"

// handler(s)
#include "Geant/PostStepActionPhysProcessHandler.h"

namespace geantphysics {

PostStepActionPhysProcessStage::PostStepActionPhysProcessStage(geant::Propagator *prop)
    : SimulationStage(geant::kPostStepActionPhysProcessStage, prop)
{
}

// base class will delete the created handlers
PostStepActionPhysProcessStage::~PostStepActionPhysProcessStage()
{
}

int PostStepActionPhysProcessStage::CreateHandlers()
{
  int threshold = fPropagator->fConfig->fNperBasket;

  auto cutsTable     = MaterialCuts::GetTheMaterialCutsTable();
  auto particleTable = Particle::GetInternalParticleTable();
  int numOfRegions   = PhysicsListManager::Instance().GetNumberOfRegions();

  for (int r = 0; r < numOfRegions; ++r) {
    fHandlersPerRegionPerParticlePerProcess.emplace_back();
    for (size_t pa = 0; pa < particleTable.size(); ++pa) {
      fHandlersPerRegionPerParticlePerProcess[r].emplace_back();
      if (particleTable[pa] == nullptr) continue;

      for (size_t p = 0; p < PhysicsProcess::GetProcessTableSize(); ++p) {
        bool particleAssigned = false;

        auto &particleProcList = PhysicsProcess::GetProcessByGlobalIndex(p)->GetListParticlesAssignedTo();
        for (auto particle : particleProcList) {
          if (particle->GetInternalCode() == particleTable[pa]->GetInternalCode()) {
            particleAssigned = true;
            break;
          }
        }
        if (particleAssigned) {
          auto handler = new PostStepActionPhysProcessHandler(threshold, fPropagator);
          fHandlersPerRegionPerParticlePerProcess[r][pa].push_back(handler);
          AddHandler(handler);
        } else {
          fHandlersPerRegionPerParticlePerProcess[r][pa].push_back(nullptr);
        }
      }
    }
  }
  return PhysicsProcess::GetProcessTableSize();
}

// Selects tracks that have any processes, any post step processes i.e. discrete part and that limited the step
geant::Handler *PostStepActionPhysProcessStage::Select(geant::Track *track, geant::TaskData * /*td*/)
{
  if (track->Status() == geant::TrackStatus_t::kPhysics && track->EIndex() == 1000) {
    // these tracks should always have psorcesses active in the given region moreover should always have discrete
    // processes that limited the step (fEindex==1000)
    assert(track->Process() >= 0 && track->Process() < (int)fHandlers.size());

    const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
        (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
    int regionIndex = matCut->GetRegionIndex();
    auto handler    = fHandlersPerRegionPerParticlePerProcess[regionIndex][track->GVcode()][track->Process()];
    assert(handler != nullptr);
    return handler;
  }
  // not physics or not discrete part of limited the step
  return nullptr;
}

} // namespace geantphysics
