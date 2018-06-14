
#include "Geant/PhysicsListManager.h"
#include "Geant/Particle.h"
#include "Geant/PostStepActionPhysModelStage.h"

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
#include "Geant/EMModel.h"
#include "Geant/EMPhysicsProcess.h"
#include "Geant/EMModelManager.h"

// handler(s)
#include "Geant/PostStepActionPhysModelHandler.h"

namespace geantphysics {

PostStepActionPhysModelStage::PostStepActionPhysModelStage(geant::Propagator *prop)
    : SimulationStage(geant::kPostStepActionStage, prop)
{
}

// base class will delete the created handlers
PostStepActionPhysModelStage::~PostStepActionPhysModelStage() {}

int PostStepActionPhysModelStage::CreateHandlers()
{
  int threshold = fPropagator->fConfig->fNperBasket;

  auto &modelTable = geantphysics::EMModel::GetGlobalTable();

  for (size_t m = 0; m < modelTable.size(); ++m) {
    auto handler = new PostStepActionPhysModelHandler(threshold, fPropagator, m);
    handler->SetMayBasketize(modelTable[m]->IsBasketizable());
    fHandlersPerModel.push_back(handler);
    AddHandler(handler);
  }

  return modelTable.size();
}

// Selects tracks that have any processes, any post step processes i.e. discrete part and that limited the step
geant::Handler *PostStepActionPhysModelStage::Select(geant::Track *track, geant::TaskData *td)
{
  if (track->Status() == geant::TrackStatus_t::kPhysics && track->EIndex() == 1000) {
    // these tracks should always have psorcesses active in the given region moreover should always have discrete
    // processes that limited the step (fEindex==1000)
    // assert(track->Process() >= 0 && track->Process() < (int)fHandlers.size());

    const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
        (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));

    int regionIndex = matCut->GetRegionIndex();

    int particleCode         = track->GVcode();
    const Particle *particle = Particle::GetParticleByInternalCode(particleCode);

    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());

    auto pProc = pManager->PostStepSelectProcess(track, td);
    if (!pProc) {
      // Delta interaction
      return nullptr;
    }

    auto emProc  = static_cast<geantphysics::EMPhysicsProcess *>(pProc);
    auto emModel = emProc->PostStepSelectModel(track->Ekin(), regionIndex);
    if (!emModel) {
      // Energy < min energy for process
      return nullptr;
    }
    if (!emModel->IsModelUsable(matCut, track->Ekin())) return nullptr;

    int modelIndex = emModel->GetGlobalIndex();
    return fHandlersPerModel[modelIndex];
  }
  // not physics or not discrete part of limited the step
  return nullptr;
}

} // namespace geantphysics
