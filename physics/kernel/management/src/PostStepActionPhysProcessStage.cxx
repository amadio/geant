
#include <Geant/PostStepActionDeltaIntHandler.h>
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
#include "Geant/EMModel.h"
#include "Geant/EMPhysicsProcess.h"
#include "Geant/EMModelManager.h"

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

  auto &modelTable = geantphysics::EMModel::GetGlobalTable();

  for (size_t m = 0; m < modelTable.size(); ++m) {
    auto handler = new PostStepActionPhysProcessHandler(threshold, fPropagator, m);
    fHandlersPerModel.push_back(handler);
    AddHandler(handler);
  }

  fDeltaIntHandler = new PostStepActionDeltaIntHandler(threshold, fPropagator);
  AddHandler(fDeltaIntHandler);

  return modelTable.size();
}

// Selects tracks that have any processes, any post step processes i.e. discrete part and that limited the step
geant::Handler *PostStepActionPhysProcessStage::Select(geant::Track *track, geant::TaskData *td)
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
    auto pProc                          = pManager->PostStepSelectProcess(track, td);
    if (!pProc) {
      return fDeltaIntHandler;
    }

    auto emProc    = (geantphysics::EMPhysicsProcess *)pProc;
    auto emModel   = emProc->GetModelManager()->SelectModel(track->E() - track->Mass(), regionIndex);
    int modelIndex = emModel->GetGlobalIndex();

    return fHandlersPerModel[modelIndex];
  }
  // not physics or not discrete part of limited the step
  return nullptr;
}

} // namespace geantphysics
