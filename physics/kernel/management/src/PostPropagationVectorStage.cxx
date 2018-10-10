
#include "Geant/PostPropagationVectorStage.h"

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
#include "Geant/MSCProcess.h"
#include "Geant/MSCModel.h"
#include "Geant/PhysicsManagerPerParticle.h"

#include "Geant/PostPropagationVectorHandler.h"

namespace geantphysics {

PostPropagationVectorStage::PostPropagationVectorStage(geant::Propagator *prop)
    : SimulationStage(geant::kPostPropagationStage, prop)
{
}

// base class will delete the created handlers
PostPropagationVectorStage::~PostPropagationVectorStage() {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PostPropagationVectorStage::PostPropagationVectorStage(const PostPropagationVectorStage &other)
    : geant::SimulationStage(other)
{
  for (auto handler : other.fHandlersPerModel)
    fHandlersPerModel.push_back(handler);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PostPropagationVectorStage &PostPropagationVectorStage::operator=(const PostPropagationVectorStage &other)
{
  geant::SimulationStage::operator=(other);
  for (auto handler : other.fHandlersPerModel)
    fHandlersPerModel.push_back(handler);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
geant::SimulationStage *PostPropagationVectorStage::Clone() const
{
  PostPropagationVectorStage *stage = new PostPropagationVectorStage(*this);
  return stage;
}

int PostPropagationVectorStage::CreateHandlers()
{
  int threshold = fPropagator->fConfig->fNperBasket;

  auto &modelTable = geantphysics::EMModel::GetGlobalTable();
  int modelsAdded  = 0;

  for (size_t m = 0; m < modelTable.size(); ++m) {
    if (!modelTable[m]->IsMSCModel()) {
      fHandlersPerModel.push_back(nullptr);
      continue;
    }
    auto handler = new PostPropagationVectorHandler(threshold, fPropagator, m);
    handler->SetMayBasketize(modelTable[m]->IsBasketizable());
    handler->SetLocal(bool(fLocalHandlers & geantphysics::kMSC));
    fHandlersPerModel.push_back(handler);
    AddHandler(handler);
    modelsAdded++;
  }

  return modelsAdded;
}

void PostPropagationVectorStage::ReplaceLocalHandlers()
{
  geant::SimulationStage::ReplaceLocalHandlers();
  auto &modelTable = geantphysics::EMModel::GetGlobalTable();
  size_t i         = 0;
  for (size_t m = 0; m < modelTable.size(); ++m) {
    if (!modelTable[m]->IsMSCModel()) continue;
    if (fHandlers[i]->IsLocal()) {
      fHandlersPerModel[m] = fHandlers[i];
    }
    i++;
  }
}

// Selects tracks that have msc process
geant::Handler *PostPropagationVectorStage::Select(geant::Track *track, geant::TaskData *td)
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
    // No any physics process for the particle or no msc process
    //
    // True path length is identical to the geometrical one so update time and number of interaction left
    // here based on fSnext
    track->IncreaseTime(track->TimeStep(track->GetStep()));
    // NOTE: we should save the previous step length in the Track and update this in the process or
    // in the process manager per particle only for the discrete processes BUT FOR THAT WE NEED TO SAVE the previous
    // step and we do it in the next step
    // track->fNintLen -= track->fStep/track->fIntLen;
    for (size_t i = 0; i < geant::kNumPhysicsProcess; ++i) {
      track->DecreasePhysicsNumOfInteractLengthLeft(i, track->GetStep() / track->GetPhysicsInteractLength(i));
    }
    return nullptr;
  }

  const MSCProcess *mscProc = static_cast<const MSCProcess *>(pManager->GetMSCProcess());
  MSCModel *model           = mscProc->SelectModel(track);
  if (!model->SamplingNeeded(track, td)) return nullptr;

  // give back the only one handler that will perform the msc step limit phase actions
  return fHandlersPerModel[model->GetGlobalIndex()];
}

} // namespace geantphysics
