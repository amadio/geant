#include "GXeBremsstrahlung.h"

FQUALIFIER
GXeBremsstrahlung::GXeBremsstrahlung(curandState* devStates,
                                     G4int threadId,
				     GPPhysicsTable* lamdaTable)
{
  // curand
  fThreadId = threadId;
  fDevStates = devStates;

  // Geant4 lambda table and variables from G4VProcess
  fLambdaTable = lamdaTable;
  fPreStepLambda = 0.0;

  // G4eBremsstrahlung Model
  isInitialised = false;
  currentModel = 0;

  // G4 variables : only for testing purpose
  theNumberOfInteractionLengthLeft = -1;
  currentInteractionLength = -1;
  theInitialNumberOfInteractionLength = -1.0;
  preStepLambda = 0.0;
  preStepKinEnergy = 0.0;
  preStepScaledEnergy = 0.0;
  mfpKinEnergy  = DBL_MAX;
  massRatio = 1.0;
  fFactor = 1.0;
  lowestKinEnergy  = 1.*eV;

}

FQUALIFIER
GXeBremsstrahlung::GXeBremsstrahlung(curandState* devStates,
				     G4int threadId,
				     G4double* xsecTable)
{
  // curand
  fThreadId = threadId;
  fDevStates = devStates;

  // tabulated cross section data
  fXsecTable = xsecTable;
  fPreStepLambda = 0.0;

  // G4eBremsstrahlung Model
  isInitialised = false;
  currentModel = 0;

  // G4 variables : only for testing purpose
  theNumberOfInteractionLengthLeft = -1;
  currentInteractionLength = -1;
  theInitialNumberOfInteractionLength = -1.0;
  preStepLambda = 0.0;
  preStepKinEnergy = 0.0;
  preStepScaledEnergy = 0.0;
  mfpKinEnergy  = DBL_MAX;
  massRatio = 1.0;

}

FQUALIFIER 
void GXeBremsstrahlung::InitialiseProcess(GXSeltzerBergerModel* model)
{
  if(!isInitialised) {

    currentModel = model;

    /*
    // energyLimit between SeltzerBerger eBremsstrahlungRel = 1*GeV;
    currentModel->SetEnergyLimitModes(1.0*GeV);

    //minKinEnergy and maxKinEnergy of G4LossTableManager
    currentModel->SetLowEnergyLimit(0.1*keV);
    currentModel->SetHighEnergyLimit(10.0*TeV);
    */

    //set ParticleChange
    currentModel->SetParticleChange(&fParticleChange);
    isInitialised = true;
  }

  G4double eth = DBL_MAX; //from G4LossTableManager
  currentModel->SetSecondaryThreshold(eth);

  // Only high energy model LMP flag is ON/OFF
  //  currentModel->SetLPMFlag(false); //for G4SeltzerBergerModel
  //  currentModel->SetLPMFlag(true); //for G4eBremsstrahlungRelModel
}

FQUALIFIER 
void GXeBremsstrahlung::ResetNumberOfInteractionLengthLeft_G4()
{
  theNumberOfInteractionLengthLeft = -log(GPUniformRand(fDevStates,fThreadId));
  theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 
}

FQUALIFIER
void GXeBremsstrahlung::StartTracking_G4()
{
  // reset parameters for the new track
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 
}

FQUALIFIER
void GXeBremsstrahlung::StartTracking_s(GXTrack* track)
{
  // initialize
  track->ns = -1.0;
}

FQUALIFIER
void GXeBremsstrahlung::StartTracking_v(G4int ntracks,
				       GXTrack* track)
{
  // initialize 
  for(G4int it = 0 ; it < ntracks ; ++it) {
    track[it].ns = -1.0;
  }
}

FQUALIFIER
void GXeBremsstrahlung::StartTracking_soa(G4int ntracks,
					  GXSoATrack track)
{
  // initialize 
  for(G4int it = 0 ; it < ntracks ; ++it) {
    (track.ns)[it] = -1.0;
  }
}

FQUALIFIER 
void GXeBremsstrahlung::DefineMaterial(GPMaterial* material)
{
  //get/set material dependencies
}

FQUALIFIER
G4double GXeBremsstrahlung::AlongStepGPIL(GPGPILSelection* selection)
{
  // condition is set to "Not Forced" 
  *selection = NotCandidateForSelection;

  // get GPIL
  return DBL_MAX;
}

FQUALIFIER
G4double GXeBremsstrahlung::PostStepGPIL_G4(GXTrack* track,
					    GPMaterial* material,
					    G4double previousStepSize,
					    GPForceCondition* condition)
{
  // condition is set to "Not Forced" 
  *condition = NotForced;
  G4double x = DBL_MAX;

  DefineMaterial(material);
  preStepKinEnergy    = track->E ; //track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;

  //@@@ we select SeltzerBergerModel and eBremsstrahlungRelModel
  //  SelectModel(preStepScaledEnergy);
  //  if(!currentModel->IsActive(preStepScaledEnergy)) { return x; }

  // compute mean free path
  if(preStepScaledEnergy < mfpKinEnergy) {
    //@@@ we use lambda table
    //    if (integral) { ComputeLambdaForScaledEnergy(preStepScaledEnergy); }
    //    else  { 
    preStepLambda = GetLambdaForScaledEnergy_G4(preStepScaledEnergy); 

    // zero cross section
    if(preStepLambda <= 0.0) { 
      theNumberOfInteractionLengthLeft = -1.0;
      currentInteractionLength = DBL_MAX;
    }
  }

  // non-zero cross section
  if(preStepLambda > 0.0) { 
    if (theNumberOfInteractionLengthLeft < 0.0) {

      // beggining of tracking (or just after DoIt of this process)
      ResetNumberOfInteractionLengthLeft_G4();

    } else if(currentInteractionLength < DBL_MAX) {
      theNumberOfInteractionLengthLeft -= 
        previousStepSize/currentInteractionLength;
      if(theNumberOfInteractionLengthLeft < 0.) {
        theNumberOfInteractionLengthLeft = 0.0;
      }
    }
    // new mean free path and step limit
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;

  }
  return x;
}

FQUALIFIER
G4double GXeBremsstrahlung::PostStepGPIL_s(GXTrack* track,
					   GPMaterial* material,
					   G4double previousStepSize,
					   GPForceCondition* condition)
{
  // condition is set to "Not Forced" 
  *condition = NotForced;

  // get/set material dependencies
  DefineMaterial(material);

  // get GPIL
  fPreStepLambda = GetLambda_s(track);

 // update the number of step length left
  if(track->ns < 0.0) {
    track->ns = -log(GPUniformRand(fDevStates,fThreadId));
  }
  else {
    track->ns -= previousStepSize/track->ci;
  }
  // new mean free path and step limit
  track->ci = 1.0/fPreStepLambda;
  return track->ns*track->ci;
}

FQUALIFIER
G4double GXeBremsstrahlung::PostStepGPIL_t(GXTrack* track,
					   GPMaterial* material,
					   G4double previousStepSize,
					   GPForceCondition* condition)
{
  // condition is set to "Not Forced" 
  *condition = NotForced;

  // get/set material dependencies
  DefineMaterial(material);

  // get GPIL
  fPreStepLambda = GetLambda_t(track);

 // update the number of step length left
  if(track->ns < 0.0) {
    track->ns = -log(GPUniformRand(fDevStates,fThreadId));
  }
  else {
    track->ns -= previousStepSize/track->ci;
  }
  // new mean free path and step limit
  track->ci = 1.0/fPreStepLambda;
  return track->ns*track->ci;
}

FQUALIFIER
void GXeBremsstrahlung::PostStepGPIL_v(G4int ntracks,
				       GXTrack* track,
				       GPMaterial* material,
				       GPForceCondition* condition)
{
  // assumptions: 
  // 1) all tracks are in the same logical volume
  // 2) energies of tracks are in the valid range of the current model
  // 3) lambda is a positive definition
  // 4) track[it].ns = -log(rand()) at the start of tracking or after post 
  //    stepping if this process is selected as the process governing the step 
  // 5) track was killed if track[it].ns < 0.0

  // condition is set to "Not Forced" for this process
  *condition = NotForced;

  // material is sharable for the data locality - not used now for one material
  DefineMaterial(material);

  GetLambda_v(ntracks,track);
  // vectorized loop if track data are properly aligned
  // simd instruction goes here
  for(G4int it = 0 ; it < ntracks ; ++it) {
    //    fPreStepLambda = GetLambda_t(track);
    //    fPreStepLambda = GetLambda(track[it].E);
    // update the number of step length left
    if(track[it].ns < 0.0) {
      track[it].ns = -log(GPUniformRand(fDevStates,fThreadId));
    }
    else {
      track[it].ns -= track[it].ps/track[it].ci;
    }
    // proposed step length and save the current interaction length
    //    track[it].ci = 1.0/fPreStepLambda;
    track[it].ci = 1.0/track[it].tv;
    track[it].s  = track[it].ns*track[it].ci; 
  }
}

FQUALIFIER
void GXeBremsstrahlung::PostStepGPIL_soa(G4int ntracks,
					 GXSoATrack track,
					 GPMaterial* material,
					 GPForceCondition* condition)
{
  // assumptions: 
  // 1) all tracks are in the same logical volume
  // 2) energies of tracks are in the valid range of the current model
  // 3) lambda is a positive definition
  // 4) track[it].ns = -log(rand()) at the start of tracking or after post 
  //    stepping if this process is selected as the process governing the step 
  // 5) track was killed if track[it].ns < 0.0

  // condition is set to "Not Forced" for this process
  *condition = NotForced;

  // material is sharable for the data locality - not used now for one material
  DefineMaterial(material);

  GetLambda_soa(ntracks,track);
  // vectorized loop if track data are properly aligned
  // simd instruction goes here
  for(G4int it = 0 ; it < ntracks ; ++it) {
    // update the number of step length left
    //    fPreStepLambda = GetLambda((track.E)[it]);
 
    if((track.ns)[it] < 0.0) {
      (track.ns)[it] = -log(GPUniformRand(fDevStates,fThreadId));
    }
    else {
      (track.ns)[it] -= (track.ps)[it]/(track.ci)[it];
    }
    // proposed step length and save the current interaction length
    //    (track.ci)[it] = 1.0/fPreStepLambda;
    (track.ci)[it] = 1.0/(track.tv)[it];
    (track.s)[it]  = (track.ns)[it]*(track.ci)[it]; 
  }
}

FQUALIFIER
GPVParticleChange& GXeBremsstrahlung::AlongStepDoIt(GXTrack* track, 
						    GPMaterial* material) 
{
  fParticleChange.ParticleChangeForLoss_InitializeForAlongStep(track);
  return fParticleChange;
}

FQUALIFIER
GPVParticleChange& GXeBremsstrahlung::PostStepDoIt_G4(GXTrack* track,
						      GPMaterial* material)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 


  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);

  G4double finalT = track->E ; //track.GetKineticEnergy();
  if(finalT <= lowestKinEnergy) { return fParticleChange; }

  G4double postStepScaledEnergy = finalT*massRatio;

  //  if(!currentModel->IsActive(postStepScaledEnergy)) { 
  //    return fParticleChange; 
  //  }

  //  G4double tcut = (*theCuts)[currentCoupleIndex]; 
  G4double tcut = 0.0893347; //MeV

  // sample secondaries
  //  secParticles.clear();

  //@@@ use SampleSecondaries of SeltzerBergerModel or eBremsstrahlungRelModel
  //    currentModel->SampleSecondaries(&secParticles, currentCouple, 
  //                                dynParticle, tcut);

  currentModel->SampleSecondaries(track, material,tcut,100*GeV);  

  // save secondaries
  G4int num = 1; // secParticles.size();
  if(num > 0) {

    //    fParticleChange.SetNumberOfSecondaries(num);

    //    for (G4int i=0; i<num; ++i) {
    //      if(secParticles[i]) {
    //        G4Track* t = new G4Track(secParticles[i], track.GetGlobalTime(), 
    //                                 track.GetPosition());
    GXTrack t = currentModel->GetSecondary();
    //        t->SetTouchableHandle(track.GetTouchableHandle());
    //        t->SetWeight(weight); 
    //    fParticleChange.AddSecondary(t);
    //      }
    //    }
  }

  if(0.0 == fParticleChange.GetProposedKineticEnergy() &&
     fAlive == fParticleChange.GetTrackStatus()) {
    //if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0){ 
    //  fParticleChange.ProposeTrackStatus(fStopButAlive); 
    //}
    //else { 
    fParticleChange.ProposeTrackStatus(fStopAndKill); 
    //}
  }

  return fParticleChange;

}

FQUALIFIER
GPVParticleChange& GXeBremsstrahlung::PostStepDoIt_s(GXTrack* track,
						     GPMaterial* material)
{
  // clear number of interaction lengths and initialize the particle change
  track->ns  = -1.0;
  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);

  //assumption
  //1) neither the bias nor the integral method are used
  //2) G4eBremsstrahlungRelModel and G4SeltzerBergerModel are selected
  //3) the particle (track) is electron
  //4) use a temporary cut - G4double tcut = (*theCuts)[currentCoupleIndex]; 
  G4double tcut = 0.0893347; //MeV

  currentModel->SampleSecondaries(track, 
				  material,
				  tcut,
				  100*GeV);  
  // save secondaries
  //  G4int num = 1; // secParticles.size();
  //  fParticleChange.SetNumberOfSecondaries(num);

  GXTrack t = currentModel->GetSecondary();
  fParticleChange.AddSecondary(t);

  return fParticleChange;
}

FQUALIFIER
GPVParticleChange& GXeBremsstrahlung::PostStepDoIt_t(GXTrack* track,
						     GPMaterial* material)
{
  // clear number of interaction lengths and initialize the particle change
  track->ns  = -1.0;
  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);

  //assumption
  //1) neither the bias nor the integral method are used
  //2) G4eBremsstrahlungRelModel and G4SeltzerBergerModel are selected
  //3) the particle (track) is electron
  //4) use a temporary cut - G4double tcut = (*theCuts)[currentCoupleIndex]; 
  G4double tcut = 0.0893347; //MeV - temporary

  currentModel->SampleSecondaries_t(track, material, tcut, 100*GeV);  

  // save secondaries

  //  G4int num = 1; // secParticles.size();
  //  fParticleChange.SetNumberOfSecondaries(num);

  GXTrack t = currentModel->GetSecondary();
  fParticleChange.AddSecondary(t);

  return fParticleChange;
}

FQUALIFIER
void GXeBremsstrahlung::PostStepDoIt_v(G4int ntracks,
				       GXTrack* track,
				       GXTrack* secTracks,
				       G4int* stackSize,
				       GPMaterial* material)
{
  // clear number of interaction lengths and initialize the particle change
  for (int it = 0; it < ntracks; ++it) {
    track[it].ns  = -1.0;
  }

  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);

  //assumption
  //1) neither the bias nor the integral method are used
  //2) G4eBremsstrahlungRelModel and G4SeltzerBergerModel are selected
  //3) the particle (track) is electron
  //4) use a temporary cut - G4double tcut = (*theCuts)[currentCoupleIndex]; 
  G4double tcut = 0.0893347; //MeV - temporary

  currentModel->SampleSecondaries_v(ntracks, track, secTracks, stackSize, 
				    material, tcut, 100*GeV);  
}

FQUALIFIER
G4double GXeBremsstrahlung::GetLambdaForScaledEnergy_G4(G4double e)
{
  //  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e);
  return fFactor*fLambdaTable->physicsVectors[1].Value(e);
}

FQUALIFIER 
G4double GXeBremsstrahlung::GetLambda_s(GXTrack* track)
{
  return fLambdaTable->physicsVectors[1].Value(track->E);
}

FQUALIFIER 
G4double GXeBremsstrahlung::GetLambda(G4double energy)
{
  size_t ibin = (log(energy/GeV)-TLEmin)*(TNBins/(TLEmax-TLEmin));
  if(ibin<0)      ibin = 0;
  if(ibin>TNBins) ibin = TNBins -1;
  return fXsecTable[ibin];
}

FQUALIFIER 
G4double GXeBremsstrahlung::GetLambda_t(GXTrack* track)
{
  size_t ibin = (log(track->E/GeV)-TLEmin)*(TNBins/(TLEmax-TLEmin));
  if(ibin<0)      ibin = 0;
  if(ibin>TNBins) ibin = TNBins -1;
  return fXsecTable[ibin];
}

FQUALIFIER 
void GXeBremsstrahlung::GetLambda_v(G4int ntracks,
				    GXTrack* track)
{
  // loop that may not be auto-vectorized
  // track[].tv is a place-holder for any transit variable which can be utilized
  // for vector instructions if necessary

  size_t ibin = 0;
  for(G4int it = 0 ; it < ntracks ; ++it) {
    ibin = (log(track[it].E/GeV)-TLEmin)*(TNBins/(TLEmax-TLEmin));
    if(ibin<0)      ibin = 0;
    if(ibin>TNBins) ibin = TNBins -1;
    track[it].tv = fXsecTable[ibin];
  }
}

FQUALIFIER 
void GXeBremsstrahlung::GetLambda_soa(G4int ntracks,
				      GXSoATrack track)
{
  // loop that may not be auto-vectorized
  // track[].tv is a place-holder for any transit variable which can be utilized
  // for vector instructions if necessary

  size_t ibin = 0;
  for(G4int it = 0 ; it < ntracks ; ++it) {
    ibin = (log((track.E)[it]/GeV)-TLEmin)*(TNBins/(TLEmax-TLEmin));
    if(ibin<0)      ibin = 0;
    if(ibin>TNBins) ibin = TNBins -1;
    (track.tv)[it] = fXsecTable[ibin]; //this line prevents vectorization
  }
}
