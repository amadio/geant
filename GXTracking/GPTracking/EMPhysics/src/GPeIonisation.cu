#include "GPeIonisation.h"
#include "GPRandom.h"

#include "GPMollerBhabhaModel.h"
#include "GPUniversalFluctuation.h"

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
GPeIonisation::GPeIonisation(curandState* devStates,
			     int threadId,
			     GPPhysicsTable* lamdaTable, 
			     GPPhysicsTable* rangeTable, 
			     GPPhysicsTable* dedxTable, 
			     GPPhysicsTable* invrTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

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
  theRangeTable = rangeTable;
  theDEDXTable = dedxTable;
  theInverseRangeTable = invrTable;

  currentModel = 0;
  isIonisation = true;
  aGPILSelection = CandidateForSelection;

  SetIonisation(isIonisation);

  //G4eIonisation
  isInitialised = false;

}

FQUALIFIER 
GPeIonisation::~GPeIonisation()
{
  ;
}

FQUALIFIER 
void GPeIonisation::InitialiseProcess(GPMollerBhabhaModel* model)
{
  if(!isInitialised) {

    //SetEmModel
    currentModel = model;

    //minKinEnergy and maxKinEnergy of G4LossTableManager
    currentModel->SetLowEnergyLimit(0.1*keV); 
    currentModel->SetHighEnergyLimit(10.0*TeV);
    
    //set ParticleChange
    currentModel->SetParticleChange(&fParticleChange);

    //if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }
    //we set the fluctuation model inside AlongStepDoIt where it is used

    isInitialised = true;
  }
}

//---------------------------------------------------------------------------
//
// G4VProcess
//
//---------------------------------------------------------------------------

FQUALIFIER
void GPeIonisation::ResetNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft = -log(GPUniformRand(fDevStates,fThreadId));
  theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 
}

void GPeIonisation::EndTracking()
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
void GPeIonisation::StartTracking()
{
  // reset parameters for the new track
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 
}

FQUALIFIER
void GPeIonisation::SetIonisation(G4bool val)
{
  isIonisation = val;
  if(val) { aGPILSelection = CandidateForSelection; }
  else    { aGPILSelection = NotCandidateForSelection; }
}

FQUALIFIER
G4bool GPeIonisation::IsIonisationProcess()
{
  return isIonisation;
}

FQUALIFIER
G4double GPeIonisation::AlongStepGetPhysicalInteractionLength(
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
G4double GPeIonisation::GetScaledRangeForScaledEnergy(G4double e)
{
  //  G4double x = ((*theRangeTableForLoss)[basedCoupleIndex])->Value(e);
  G4double x =  theRangeTable->physicsVectors[1].Value(e);
  if(e < minKinEnergy) { x *= sqrt(e/minKinEnergy); }
  return x;
}

FQUALIFIER
G4double GPeIonisation::PostStepGetPhysicalInteractionLength(GXTrack* track,
					          G4double previousStepSize,
					          GPForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double x = DBL_MAX;

  DefineMaterial();
  preStepKinEnergy    = track->E ; // track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;

  //@@@ we select G4MollerBhabhaModel for G4eIonisation model
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

#ifdef GPUPLOTS
#ifndef __CUDA_ARCH__
    hmgr.getHisto("hioniNbOfIntLengthLeft").Fill(theNumberOfInteractionLengthLeft);
    hmgr.getHisto("hioniPreStepScaledEnergy").Fill(preStepScaledEnergy);
    hmgr.getHisto("hioniPreStepLambda").Fill(preStepLambda);
    hmgr.getHisto("hioniMassRatio").Fill(massRatio);
#endif
#endif

  return x;
}

FQUALIFIER
void GPeIonisation::DefineMaterial()
{
  //fFactor = chargeSqRatio*biasFactor*(*theDensityFactor)[currentCoupleIndex];
  fFactor = 1.0;
  reduceFactor = 1.0/(fFactor*massRatio);
  mfpKinEnergy = DBL_MAX;
}

FQUALIFIER
G4double GPeIonisation::GetLambdaForScaledEnergy(G4double e)
{
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPeIonisation.cu: Ene="<< e <<" fFactor="<< fFactor
	   <<" lambda(E)="<< theLambdaTable->physicsVectors[1].Value(e)
	   << std::endl;
#endif // __CUDA_ARCH__
#endif // GPUDEBUG
  //  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e);
  return fFactor* theLambdaTable->physicsVectors[1].Value(e);
}

//G4VParticleChange* G4VEnergyLossProcess::AlongStepDoIt(const G4Track& track,
//                                                       const G4Step& step)
FQUALIFIER
GPVParticleChange& GPeIonisation::AlongStepDoIt(GXTrack* track, 
						GPMaterial* material, 
						G4double stepLength)
{
  //@@@G4FWP - redefine Material due to the split GPIL and DoIt in 
  //           the kernel level: remove this by using GXTrackLiason
  DefineMaterial();
  preStepKinEnergy    = track->E ; // track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  //@@@G4FWP

  fParticleChange.ParticleChangeForLoss_InitializeForAlongStep(track);
  // The process has range table - calculate energy loss
  if(!isIonisation || !currentModel->IsActive(preStepScaledEnergy)) {
    return fParticleChange;
  }

  // Get the actual (true) Step length
  //  G4double length = step.GetStepLength();
  //  G4double length = stepLength;
  G4double length = track->s;
  if(length <= 0.0) { return fParticleChange; }
  G4double eloss  = 0.0;

  //  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  // define new weight for primary and secondaries
  //  G4double weight = fParticleChange.GetParentWeight();
  //  if(weightFlag) {
  //    weight /= biasFactor;
  //    fParticleChange.ProposeWeight(weight);
  //  }

  // stopping
  if (length >= fRange) {
    eloss = preStepKinEnergy;
    //@@@check this out
    //    if (useDeexcitation) {
    //      atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
    //                                              eloss, currentCoupleIndex);
    //      if(scTracks.size() > 0) { FillSecondariesAlongStep(eloss, weight); }
    //      if(eloss < 0.0) { eloss = 0.0; }
    //    }
    //    fParticleChange.SetProposedKineticEnergy(0.0);
    //    fParticleChange.ProposeLocalEnergyDeposit(eloss);
    //    return &fParticleChange;
  }
  // Short step
  eloss = GetDEDXForScaledEnergy(preStepScaledEnergy)*length;
#ifdef GPUPLOTS
#ifndef __CUDA_ARCH__
  GPHistoManager& hmgr = GPHistoManager::getInstance();
  hmgr.getHisto("hioniDedxForScaledEnergyTimesLength").Fill(eloss);
#endif
#endif

  // Long step
  G4double linLossLimit = 0.01; //used once at here
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPeIonisation.cu: Eloss="<< eloss
	   <<" preStepKinE * linLossLimit="<< preStepKinEnergy*linLossLimit
	   <<" - Eloss required to be larger for ElossFromKinE-ScaledE calculation!"
	   << std::endl;
#endif
#endif
  if(eloss > preStepKinEnergy*linLossLimit) {
    //G4double x = GetScaledRangeForScaledEnergy(preStepScaledEnergy) 
    //  - length/reduceFactor;
    G4double x = (fRange - length)/reduceFactor;
    eloss = preStepKinEnergy - ScaledKinEnergyForLoss(x)/massRatio;
#ifndef __CUDA_ARCH__
#ifdef GPUDEBUG
  std::cout<<" fRange="<< fRange <<" length="<< length <<" reduceFactor="<< reduceFactor
	   <<" x="<< x <<" commentedX="<< (fRange - length/reduceFactor)
	   <<" Eloss="<< eloss << std::endl;
#endif
#ifdef GPUPLOTS
    hmgr.getHisto("hioniElossFromKinEnergyMinusScaledEnergyForLoss").Fill(eloss);
#endif
#endif // __CUDA_ARCH__
  }

  //  G4double cut  = (*theCuts)[currentCoupleIndex];
  G4double esec = 0.0;

  //useSubCutoff = false by default
  /*
  if(useSubCutoff) {
    if(idxSCoffRegions[currentCoupleIndex]) {

      G4bool yes = false;
      G4StepPoint* prePoint = step.GetPreStepPoint();

      // Check boundary
      if(prePoint->GetStepStatus() == fGeomBoundary) { yes = true; }

      // Check PrePoint
      else {
        G4double preSafety  = prePoint->GetSafety();
        G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);

        // recompute presafety
        if(preSafety < rcut) {
          preSafety = safetyHelper->ComputeSafety(prePoint->GetPosition());
        }

        if(preSafety < rcut) { yes = true; }

        // Check PostPoint
        else {
          G4double postSafety = preSafety - length; 
          if(postSafety < rcut) {
            postSafety = 
              safetyHelper->ComputeSafety(step.GetPostStepPoint()->GetPosition()
);
            if(postSafety < rcut) { yes = true; }
          }
        }
      }
  
      // Decided to start subcut sampling
      if(yes) {

        cut = (*theSubCuts)[currentCoupleIndex];
        eloss -= GetSubDEDXForScaledEnergy(preStepScaledEnergy)*length;
        esec = SampleSubCutSecondaries(scTracks, step, 
                                       currentModel,currentCoupleIndex);
      }   
    }
  }
  */

  // Sample fluctuations
  // @@@ lossFluctuationFlag = true
  // currentModel->GetModelOfFluctuations() -> G4UniversalFluctuation()
  // see G4eIonisation::InitialiseEnergyLossProcess

  //  if (lossFluctuationFlag) {
  //    G4VEmFluctuationModel* fluc = currentModel->GetModelOfFluctuations();
  //    if(fluc && 
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPeIonisation.cu: Eloss="<< eloss
	   <<" esec="<< esec <<" and lowestKinE="<< lowestKinEnergy
	   <<" --> sum="<< (eloss + esec + lowestKinEnergy)
	   <<" vs. preStepKinE="<< preStepKinEnergy << std::endl;
  std::cout<<"GPeIoni.cu: If above sum < preStepKinE - then sample fluctuations will be calculated\n";
#endif
#endif
  if(  (eloss + esec + lowestKinEnergy) < preStepKinEnergy) {

    G4double tmax = 0.5*preStepKinEnergy;
    //        min(currentModel->MaxSecondaryKinEnergy(dynParticle),cut);
    G4double emean = eloss;
    GPUniversalFluctuation fluc(fDevStates,fThreadId);
    eloss = fluc.SampleFluctuations(material,preStepKinEnergy,
					  tmax,length,emean);
#ifndef __CUDA_ARCH__
#ifdef GPUDEBUG
    std::cout<<"GPeIoni.cu: tmax="<< tmax
	     <<" emean="<< emean <<" Eloss=" << eloss << std::endl;
#endif
#ifdef GPUPLOTS
    hmgr.getHisto("hioniElossFromSampleFluctuations").Fill(eloss);
#endif
#endif // __CUDA_ARCH__
  }

  //  }

  // deexcitation
  // useDeexcitation = false unless atomDeexcitation->IsPIXEActive() = true
  /*
  if (useDeexcitation) {
    G4double esecfluo = preStepKinEnergy - esec;
    G4double de = esecfluo;
    atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
                                            de, currentCoupleIndex);

    // sum of de-excitation energies
    esecfluo -= de;

    // subtracted from energy loss
    if(eloss >= esecfluo) {
      esec  += esecfluo;
      eloss -= esecfluo;
    } else {
      esec += esecfluo;
      eloss = 0.0; 
    } 
  }
  if(scTracks.size() > 0) { FillSecondariesAlongStep(eloss, weight); }
  */

  // Energy balance
  G4double finalT = preStepKinEnergy - eloss - esec;
  if (finalT <= lowestKinEnergy) {
    eloss += finalT;
    finalT = 0.0;
  } 

  if(eloss < 0.0) { eloss = 0.0; }
  fParticleChange.SetProposedKineticEnergy(finalT);
  fParticleChange.ProposeLocalEnergyDeposit(eloss);
#ifndef __CUDA_ARCH__
#ifdef GPUDEBUG
  std::cout<<"GPeIoni.cu: finalT="<< finalT <<" Eloss="<< eloss << std::endl;
#endif
#ifdef GPUPLOTS
  hmgr.getHisto("hioniEloss").Fill(eloss);
#endif
#endif // __CUDA_ARCH__
  return fParticleChange;
}

FQUALIFIER
G4double GPeIonisation::GetDEDXForScaledEnergy(G4double e)
{
  //  G4double x = fFactor*(*theDEDXTable)[basedCoupleIndex]->Value(e);
  G4double x = fFactor* theDEDXTable->physicsVectors[1].Value(e);
  if(e < minKinEnergy) { x *= sqrt(e/minKinEnergy); }
  return x;
}

FQUALIFIER
G4double GPeIonisation::ScaledKinEnergyForLoss(G4double r)
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
GPVParticleChange& GPeIonisation::PostStepDoIt(GXTrack* track, 
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
  // secParticles.clear();

  //  currentModel->SampleSecondaries(&secParticles, currentCouple, 
  //				  dynParticle, tcut);
  currentModel->SampleSecondaries(track, material, tcut,100*GeV);  

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
