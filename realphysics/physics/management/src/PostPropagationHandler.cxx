
#include "PostPropagationHandler.h"

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


namespace geantphysics {

PostPropagationHandler::PostPropagationHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


PostPropagationHandler::~PostPropagationHandler() {}


// The PostPropagationStage will select only tracks with particles that (1) has any physics processes
// active in the given region and (2) has msc process
void PostPropagationHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td) {
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
  // invoke the along step doit method
  mscProc->AlongStepDoIt(track, td);
  // Particles have been propagated by the geometrical step and after the msc AlongStepDoIt this geometrical step
  // is converted back to true step length. Time and number of interaction left must be updated by using this true
  // path length (that have been written into fStep by msc).  
  track->fTime += track->TimeStep(track->fStep);
  track->fNintLen -= track->fStep/track->fIntLen;

  // --
  // copy the input track to the output
  // (secondaries should be copied as well but there is no secondary production in the msc along-step-action)
  output.AddTrack(track);
}


}  // namespace geantphysics
