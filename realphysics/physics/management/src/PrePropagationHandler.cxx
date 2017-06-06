
#include "PrePropagationHandler.h"

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
#include "MSCProcess.h"

#include "Geant/NavigationInterface.h"


namespace geantphysics {

PrePropagationHandler::PrePropagationHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


PrePropagationHandler::~PrePropagationHandler() {}


// The PrePropagationStage will select only tracks with particles that (1) has any physics processes
// active in the given region and (2) has msc process
void PrePropagationHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td) {
  // ---
  int numSecondaries = 0;
  // here we will get the MaterialCuts from the LogicalVolume later
  int   matIndx              = track->GetMaterial()->GetIndex();
  int   regIndx              = const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetRegion()->GetIndex();
  const MaterialCuts *matCut =  MaterialCuts::GetMaterialCut(regIndx,matIndx);
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // put some asserts here to make sure (1) that the partcile has any processes, (2) the particle has at least one
  // process with continuous parts.
  assert(pManager!=nullptr); // (1)
  assert(pManager->HasMSCProcess()); // (2)
  //
  // invoke the step limit part of the msc process
  const MSCProcess *mscProc = static_cast<const MSCProcess*>(pManager->GetMSCProcess());
  assert(mscProc);  // make sure that it is not nullptr
  // invoke the step limit method
  Geant::ScalarNavInterfaceVGM::NavFindNextBoundaryForMSC(*track, track->fPstep);
  mscProc->AlongStepLimitationLength(track, td);
  // --
  // copy the input track to the output
  // (secondaries should be copied as well but there is no secondary production in the msc along-step limit)
  output.AddTrack(track);
}


}  // namespace geantphysics
