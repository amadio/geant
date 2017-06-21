
#include "ComputeIntLHandler.h"

// from geantV
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "Basket.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"

namespace geantphysics {

ComputeIntLHandler::ComputeIntLHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


ComputeIntLHandler::~ComputeIntLHandler() {}

//
// The select method of the ComputeIntLStage selected only those tracks that (1) that has any process assigned i.e.
// their PhysicsManagerPerParticle object is not null in the given region and (2) they has AlongStep and/or PostStep
// processes that can limit the step.
void ComputeIntLHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td) {
  // reset step length and energy deposit
  track->fStep=0.; // no setter for this member in GeantTrack
  track->SetEdep(0.);
  // ---
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // put some asserts here to make sure (1) that the partcile has any processes, (2) the particle has at least one
  // process with continuous or discrete parts.
  assert(pManager!=nullptr); // (1)
  assert(pManager->GetListAlongStepProcesses().size()+pManager->GetListPostStepCandidateProcesses().size()!=0);
  //
  // compute the intercation length:
  LightTrack primaryLT;
  //LightTrack *lt = &(td->fPhysicsData->GetListOfSecondaries()[0]);
  // we will use members:
  //  fNintLen      <==>  fNintLen  // number of interaction left
  //  fTargetZ      <==>  fEindex   // will be set to flag if disc. or cont. step limit won
  //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
  //  fKinE         <==>  fE-fMass  // kinetic energy
  //  fGVcode       <==>  fGVcode   // internal particle code
  //  fIntLen       <==>  fIntLen   // will be set to store the current inverse total lambda
  //  fStepLength   <==>  fPstep    // will be set to store the physics step limit
  primaryLT.SetNumOfInteractionLegthLeft(track->GetNintLen());
  primaryLT.SetKinE(track->E()-track->Mass());
  primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
  primaryLT.SetGVcode(track->GVcode());
  pManager->ComputeIntLen(primaryLT,td);
  // update GeantTrack
  track->SetPstep(primaryLT.GetStepLength());
  track->SetIntLen(primaryLT.GetTotalMFP());
  track->SetEindex(primaryLT.GetTargetZ()); // just indicates if along or post step limit happened
  if (track->GetNintLen()<=0.0) { // was resampled in pManager->ComputeIntLen()
    track->SetNintLen(primaryLT.GetNumOfInteractionLegthLeft());
  }
  // ---
  // copy input track to the output
  output.AddTrack(track);
}


} // namespace geantphysics
