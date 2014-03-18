#include "GPeBremsstrahlung.h"
#include "GPRandom.h"

#include "GPSeltzerBergerRelModel.h"

#include "stdio.h"

#ifndef __CUDA_ARCH__
#ifdef GPUPLOTS
  #include "GPHistoManager.hh"
#endif
#ifdef GPUDEBUG
  #include <iostream>
#endif
#endif // __CUDA_ARCH__


FQUALIFIER
GPeBremsstrahlung::GPeBremsstrahlung(curandState* devStates,
				     int threadId,
				     GPPhysicsTable* lamdaTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  // G4eBremsstrahlung
  isInitialised = false;

  // G4VProcess
  theNumberOfInteractionLengthLeft = -1;
  currentInteractionLength = -1;
  theInitialNumberOfInteractionLength = -1.0;

  //G4VEnergyLossProcess
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
  fRange        = DBL_MAX;
  massRatio = 1.0;
  fFactor = 1.0;
  reduceFactor = 1.0;

  lowestKinEnergy  = 1.*eV;

  // default dRoverRange and finalRange
  //  SetStepFunction(0.2, 1.0*mm);
  dRoverRange = 0.2;
  finalRange = 1.0*mm;
  minKinEnergy = 0.1*keV;

  //model specific
  theLambdaTable = lamdaTable;
  theRangeTable = 0;
  theDEDXTable = 0;
  theInverseRangeTable = 0;

  currentModel = 0;
  isIonisation = false;
  aGPILSelection = CandidateForSelection;

  SetIonisation(isIonisation);

}

FQUALIFIER GPeBremsstrahlung::~GPeBremsstrahlung()
{
}

//G4eBremsstrahlung::InitialiseEnergyLossProcess(const G4ParticleDefinition*,
//                                               const G4ParticleDefinition*)
FQUALIFIER void 
GPeBremsstrahlung::InitialiseProcess(GPSeltzerBergerRelModel* model)
{
  if(!isInitialised) {

    currentModel = model;

    // energyLimit between SeltzerBerger eBremsstrahlungRel = 1*GeV;
    currentModel->SetEnergyLimitModes(1.0*GeV);

    //minKinEnergy and maxKinEnergy of G4LossTableManager
    currentModel->SetLowEnergyLimit(0.1*keV);
    currentModel->SetHighEnergyLimit(10.0*TeV);

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

//---------------------------------------------------------------------------
//
// G4VProcess
//
//---------------------------------------------------------------------------

FQUALIFIER
void GPeBremsstrahlung::ResetNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft = -log(GPUniformRand(fDevStates,fThreadId));
  theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 
}

FQUALIFIER
void GPeBremsstrahlung::EndTracking()
{
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = -1.0;
  theInitialNumberOfInteractionLength=-1.0;
}

//---------------------------------------------------------------------------
//
// G4VEnergyLossProcess:
//
//---------------------------------------------------------------------------

FQUALIFIER
void GPeBremsstrahlung::StartTracking()
{
  // reset parameters for the new track
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 
}

FQUALIFIER
void GPeBremsstrahlung::SetIonisation(G4bool val)
{
  isIonisation = val;
  if(val) { aGPILSelection = CandidateForSelection; }
  else    { aGPILSelection = NotCandidateForSelection; }
}

FQUALIFIER
G4bool GPeBremsstrahlung::IsIonisationProcess()
{
  return isIonisation;
}

FQUALIFIER
G4double GPeBremsstrahlung::AlongStepGetPhysicalInteractionLength(
					        GPGPILSelection* selection)
{
  G4double x = DBL_MAX;
  *selection = aGPILSelection;

  if(isIonisation) {
    fRange = GetScaledRangeForScaledEnergy(preStepScaledEnergy)*reduceFactor;

    x = fRange;
    G4double y = x*dRoverRange;
    G4double finR = finalRange;
    //@@@ rndmStepFlag=false
    //if(rndmStepFlag) { 
    //finR = min(finR,currentCouple->GetProductionCuts()->GetProductionCut(1)); 
    //}
    if(x > finR) { x = y + finR*(1.0 - dRoverRange)*(2.0 - finR/fRange); }
  }

  return x;
}

FQUALIFIER
G4double GPeBremsstrahlung::GetScaledRangeForScaledEnergy(G4double e)
{
  //  G4double x = ((*theRangeTableForLoss)[basedCoupleIndex])->Value(e);
  G4double x =  theRangeTable->physicsVectors[1].Value(e);
  if(e < minKinEnergy) { x *= sqrt(e/minKinEnergy); }
  return x;
}

FQUALIFIER
G4double GPeBremsstrahlung::PostStepGetPhysicalInteractionLength(GXTrack* track,
 					          G4double previousStepSize,
					          GPForceCondition* condition)
{
  // condition is set to "Not Forced" 
  *condition = NotForced;
  G4double x = DBL_MAX;

  DefineMaterial();
  preStepKinEnergy    = track->E ; //track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;

  //@@@ we select SeltzerBergerModel and eBremsstrahlungRelModel
  //  SelectModel(preStepScaledEnergy);
  if(!currentModel->IsActive(preStepScaledEnergy)) { return x; }

#ifdef GPUPLOTS
#ifndef __CUDA_ARCH__
  GPHistoManager& hmgr = GPHistoManager::getInstance();
#endif
#endif

  // compute mean free path
  if(preStepScaledEnergy < mfpKinEnergy) {
    //@@@ we use lambda table
    //    if (integral) { ComputeLambdaForScaledEnergy(preStepScaledEnergy); }
    //    else  { 
    preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy); 
#ifdef GPUPLOTS
#ifndef __CUDA_ARCH__
    hmgr.getHisto("hbremPreStepLambda").Fill(preStepLambda);
#endif
#endif
    //}

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
      ResetNumberOfInteractionLengthLeft();

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

#ifdef  GPUPLOTS
#ifndef __CUDA_ARCH__
    hmgr.getHisto("hbremNbOfIntLengthLeft").Fill(theNumberOfInteractionLengthLeft);
    hmgr.getHisto("hbremPreStepScaledEnergy").Fill(preStepScaledEnergy);
#endif
#endif

  return x;
}

FQUALIFIER
void GPeBremsstrahlung::DefineMaterial()
{
  //fFactor = chargeSqRatio*biasFactor*(*theDensityFactor)[currentCoupleIndex];
  fFactor = 1.0;
  reduceFactor = 1.0/(fFactor*massRatio);
  mfpKinEnergy = DBL_MAX;
}

FQUALIFIER
G4double GPeBremsstrahlung::GetLambdaForScaledEnergy(G4double e)
{
  //  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e);
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPeBremsstrahlung.cu: "<< e <<' '<< fFactor <<' '<< theLambdaTable->physicsVectors[1].Value(e)
	   << std::endl;
#endif // __CUDA_ARCH__
#endif // GPUDEBUG
  return fFactor* theLambdaTable->physicsVectors[1].Value(e);
}

//G4VParticleChange* G4VEnergyLossProcess::AlongStepDoIt(const G4Track& track,
//                                                       const G4Step& step)
FQUALIFIER
GPVParticleChange& GPeBremsstrahlung::AlongStepDoIt(GXTrack* track, 
						    GPMaterial* material, 
						    G4double stepLength)
{
  fParticleChange.ParticleChangeForLoss_InitializeForAlongStep(track);
  // The process has range table - calculate energy loss
  if(!isIonisation || !currentModel->IsActive(preStepScaledEnergy)) {
    return fParticleChange;
  }
  //@@@ isIonisation = false
  return fParticleChange;
}

FQUALIFIER
G4double GPeBremsstrahlung::GetDEDXForScaledEnergy(G4double e)
{
  //  G4double x = fFactor*(*theDEDXTable)[basedCoupleIndex]->Value(e);
  G4double x = fFactor* theDEDXTable->physicsVectors[1].Value(e);
  if(e < minKinEnergy) { x *= sqrt(e/minKinEnergy); }
  return x;
}

FQUALIFIER
G4double GPeBremsstrahlung::ScaledKinEnergyForLoss(G4double r)
{
  //  G4PhysicsVector* v = (*theInverseRangeTable)[basedCoupleIndex];
  //  G4double rmin = v->Energy(0);
  G4double rmin = theInverseRangeTable->physicsVectors[1].Energy(0);
  G4double e = 0.0; 
  if(r >= rmin) { 
    //    e = v->Value(r); 
    e = theInverseRangeTable->physicsVectors[1].InvRangeValue(r); 
  }
  else if(r > 0.0) {
    G4double x = r/rmin;
    e = minKinEnergy*x*x;
  }
  return e;
}

//G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
//                                                      const G4Step& step)
FQUALIFIER
GPVParticleChange& GPeBremsstrahlung::PostStepDoIt(GXTrack* track, 
						   GPMaterial* material)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 

  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);

  G4double finalT = track->E ; //track.GetKineticEnergy();
  if(finalT <= lowestKinEnergy) { return fParticleChange; }

  G4double postStepScaledEnergy = finalT*massRatio;

  if(!currentModel->IsActive(postStepScaledEnergy)) { 
    return fParticleChange; 
  }

  //@@@ we use neither the bias nor the integral method
  /*
  // Integral approach
  if (integral) {
    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy);
    if(lx <= 0.0) {
      return &fParticleChange;
    } else if(preStepLambda*G4UniformRand() > lx) {
      return &fParticleChange;
    }
  }
  */

  //@@@ we select G4eBremsstrahlungRelModel and G4SeltzerBergerModel
  //  SelectModel(postStepScaledEnergy);

  // define new weight for primary and secondaries
  //  G4double weight = fParticleChange.GetParentWeight();
  //  if(weightFlag) {
  //    weight /= biasFactor;
  //    fParticleChange.ProposeWeight(weight);
  //  }

  //@@@ the particle is electron
  //  const G4DynamicParticle* dynParticle = track.GetDynamicParticle(); 
  //@@@ use the default cut = 1mm
  //  G4double tcut = (*theCuts)[currentCoupleIndex]; 
  G4double tcut = 0.0893347; //MeV

  // sample secondaries
  //  secParticles.clear();

  //@@@ use SampleSecondaries of SeltzerBergerModel or eBremsstrahlungRelModel
  //    currentModel->SampleSecondaries(&secParticles, currentCouple, 
  //				    dynParticle, tcut);

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
    fParticleChange.AddSecondary(t);
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
