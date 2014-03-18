#include "GPPhotonProcess.h"

#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"
#include "GPPhysicsTable.h"
#include "GPVParticleChange.h"
#include "GPRandom.h"

FQUALIFIER
GPPhotonProcess::GPPhotonProcess() {}

FQUALIFIER
GPPhotonProcess::GPPhotonProcess(curandState* devStates,
				 int threadId,
				 GPPhysicsTable* lamdaTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  //Process and model
  isInitialised = false;
  theProcessType = kNullType;
  currentModel = 0;

  theLambdaTable = lamdaTable;
  theDensityFactor = 0;
  theDensityIdx = 0;

  // Size of tables assuming spline
  minKinEnergy = 0.1*keV;
  maxKinEnergy = 10.0*TeV;
  nLambdaBins  = 77;
  lambdaFactor  = 0.8;
  polarAngleLimit = 0.0;

  currentMaterial = 0;

  mfpKinEnergy  = DBL_MAX;
  preStepKinEnergy = 0.0;
  preStepLambda = 0.0;
  fFactor = 1.0;

}

FQUALIFIER
GPPhotonProcess::~GPPhotonProcess()
{
  ;
}

FQUALIFIER
GPEmProcessType GPPhotonProcess::GetProcessType()
{
  return theProcessType;
}

FQUALIFIER
void GPPhotonProcess::Clear()
{
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
}

FQUALIFIER 
void GPPhotonProcess::InitialiseProcess(GPEmProcessType kType, 
					GPPhotonModel *aModel)
{
  if(!isInitialised) {
    //@@@ Process sepecific implementation
    //SetEmModel
    theProcessType = kType;
    currentModel = aModel;

    currentModel->SetProcessType(theProcessType);

    if(theProcessType == kConversion) {
      const G4double limit = 80*GeV;
      G4double emin = fmax(MinKinEnergy(), 2*electron_mass_c2);
      G4double emax = MaxKinEnergy();
      SetMinKinEnergy(emin);

      currentModel->SetLowEnergyLimit(emin);
      G4double ehigh = fmin(emax,limit);
      ehigh = fmin(ehigh,currentModel->HighEnergyLimit());

      currentModel->SetHighEnergyLimit(ehigh);

      if(emax > ehigh) {
	//currentModel should be G4PairProductionRelModel - not implemented yet
	//currentModel->SetLowEnergyLimit(ehigh);
	//currentModel->SetHighEnergyLimit(emax);
      }
    }
    else { // theProcessType == kCompton || theProcessType == kPhotoElectric
      currentModel->SetLowEnergyLimit(MinKinEnergy());
      currentModel->SetHighEnergyLimit(MaxKinEnergy());
    }

    currentModel->SetParticleChange(&fParticleChange);

    isInitialised = true;
  }
}

FQUALIFIER
void GPPhotonProcess::StartTracking(GXTrack* track)
{
  // reset parameters for the new track
  //  currentParticle = track->GetParticleDefinition();
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 
}

FQUALIFIER
void GPPhotonProcess::StartStepping()
{
  //@@@G4FWP 
  //For multiple steppings, this should be invoked at the beginning of 
  //each step to clean up secondaries produced in the previous step.
  fParticleChange.Clear();
}

FQUALIFIER
void GPPhotonProcess::SetLambdaTable(GPPhysicsTable* table)
{
  theLambdaTable = table;
}

FQUALIFIER
G4double GPPhotonProcess::PostStepGetPhysicalInteractionLength(
                             GXTrack* track,
                             GPMaterial* material,
                             G4double   previousStepSize,
                             GPForceCondition* condition)
{
  *condition = NotForced;
  G4double x = DBL_MAX;

  preStepKinEnergy = track->E; //track.GetKineticEnergy();
  DefineMaterial(material);
  //@@@  SelectModel(preStepKinEnergy, currentCoupleIndex);

  if(!currentModel->IsActive(preStepKinEnergy)) { return x; }

   // compute mean free path
  if(preStepKinEnergy < mfpKinEnergy) {
    //    if (integral) { ComputeIntegralLambda(preStepKinEnergy); }
    //    else { 
    preStepLambda = GetCurrentLambda(preStepKinEnergy); 
    //    }

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

      // subtract NumberOfInteractionLengthLeft using previous step
      theNumberOfInteractionLengthLeft -= previousStepSize/currentInteractionLength;
      //SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft < 0.) {
	theNumberOfInteractionLengthLeft = 0.0;
      }
    }

    // new mean free path and step limit for the next step
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;
 }
  return x;
}

FQUALIFIER
GPVParticleChange& GPPhotonProcess::PostStepDoIt(GXTrack* track,
						 GPMaterial* material)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 

  //@@@implement for particle change for a specific model
  //  fParticleChange.Clear();
  fParticleChange.ParticleChangeForGamma_InitializeForPostStep(track);

  // Do not make anything if particle is stopped, the annihilation then
  // should be performed by the AtRestDoIt!

  //@@@  if (track.GetTrackStatus() == fStopButAlive) { return &fParticleChange; }

  G4double finalT = track->E;//track.GetKineticEnergy();

  // Integral approach
  /*
  if (integral) {
    G4double lx = GetLambda(finalT, currentCouple);
    if(preStepLambda<lx && 1 < verboseLevel) {
      //warning
    }

    if(preStepLambda*G4UniformRand() > lx) {
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChange;
    }
  }
  */

  //@@@  SelectModel(finalT, currentCoupleIndex);
  if(!currentModel->IsActive(finalT)) { return fParticleChange; }

  // sample secondaries
  //  secParticles.clear();
  //  currentModel->SampleSecondaries(&secParticles, currentCouple, 
  //				  track.GetDynamicParticle(),
  //				  (*theCuts)[currentCoupleIndex]);
  currentModel->SampleSecondaries(track, material,1.0*mm);

  // save secondaries
  //  G4int num = secParticles.size();
  G4int num =1; //@@@ = 2 for GammaConversion
  if(num > 0) {

    //    fParticleChange.SetNumberOfSecondaries(num);
    G4double edep = fParticleChange.GetLocalEnergyDeposit();
    /*
    for (G4int i=0; i<num; ++i) {
      if (secParticles[i]) {
        G4DynamicParticle* dp = secParticles[i];
        const G4ParticleDefinition* p = dp->GetParticleDefinition();
        G4double e = dp->GetKineticEnergy();
        G4bool good = true;
        if(applyCuts) {
	  if (p == theGamma) {
	    if (e < (*theCutsGamma)[currentCoupleIndex]) { good = false; }
	    
	  } else if (p == theElectron) {
	    if (e < (*theCutsElectron)[currentCoupleIndex]) { good = false; }
	    
	  } else if (p == thePositron) {
	    if (electron_mass_c2 < (*theCutsGamma)[currentCoupleIndex] &&
		e < (*theCutsPositron)[currentCoupleIndex]) {
	      good = false;
	      e += 2.0*electron_mass_c2;
	    }
	  }
	  // added secondary if it is good
        }
        if (good) { 
          G4Track* t = new G4Track(dp, track.GetGlobalTime(),track.GetPosition());
          t->SetTouchableHandle(track.GetTouchableHandle());
          t->SetWeight(weight);
          pParticleChange->AddSecondary(t); 
        } else {
	  delete dp;
	  edep += e;
	}
      } 
    }
    */
    //@@@ a secondary electron produced by a model and 
    //    additional secondary positron by conversion
    GXTrack electron = currentModel->GetSecondaryElectron();
    fParticleChange.AddSecondary(electron);

    if(theProcessType == kConversion) {
      GXTrack positron = currentModel->GetSecondaryPositron();
      fParticleChange.AddSecondary(positron);
    }
      
    fParticleChange.ProposeLocalEnergyDeposit(edep);
  }
  
  if(0.0 == fParticleChange.GetProposedKineticEnergy() &&
     fAlive == fParticleChange.GetTrackStatus()) {
    //  if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
    //       { fParticleChange.ProposeTrackStatus(fStopButAlive); }
    //  else { 
    fParticleChange.ProposeTrackStatus(fStopAndKill); 
    //}
  }

  return fParticleChange;
}

FQUALIFIER
G4double 
GPPhotonProcess::CrossSectionPerVolume(G4double kineticEnergy,
				    GPMaterial* material)
{
  // Cross section per atom is calculated
  DefineMaterial(material);
  G4double cross = 0.0;
  if(theLambdaTable) {
    //    cross = (*theDensityFactor)[currentCoupleIndex]*
    //      (((*theLambdaTable)[basedCoupleIndex])->Value(kineticEnergy));
    cross = 1.0*
      theLambdaTable->physicsVectors[1].Value(kineticEnergy);

  } else {
    //    SelectModel(kineticEnergy, currentCoupleIndex);
    cross = currentModel->CrossSectionPerVolume(currentMaterial,
						kineticEnergy);
  }

  if(cross < 0.0) { cross = 0.0; }
  return cross;
}

FQUALIFIER
G4double GPPhotonProcess::GetMeanFreePath(GPMaterial* material,
				       G4double kineticEnergy,
				       GPForceCondition* condition)
{
  *condition = NotForced;
  return GPPhotonProcess::MeanFreePath(material,kineticEnergy);
}

FQUALIFIER
G4double GPPhotonProcess::MeanFreePath(GPMaterial* material,
				    G4double kineticEnergy)
{
  //  DefineMaterial(track.GetMaterialCutsCouple());
  DefineMaterial(material);
  preStepLambda = GetCurrentLambda(kineticEnergy);
  G4double x = DBL_MAX;
  if(0.0 < preStepLambda) { x = 1.0/preStepLambda; }
  return x;
}

FQUALIFIER
G4double 
GPPhotonProcess::ComputeCrossSectionPerAtom(G4double kineticEnergy, 
					 G4double Z, G4double A, G4double cut)
{
  // SelectModel(kineticEnergy, currentCoupleIndex);
  G4double x = 0.0;
  if(currentModel) {
    x = currentModel->ComputeCrossSectionPerAtom(kineticEnergy,
						 Z,A,cut);
  }
  return x;
}

FQUALIFIER
GPElement* GPPhotonProcess::GetCurrentElement()
{
  GPElement* elm = 0;
  if(currentModel) {elm = currentModel->GetCurrentElement(); }
  return elm;
}


FQUALIFIER
void GPPhotonProcess::SetMinKinEnergy(G4double e)
{
  nLambdaBins = GPlrint(nLambdaBins*log(maxKinEnergy/e)
			/log(maxKinEnergy/minKinEnergy));
  minKinEnergy = e;
}

FQUALIFIER
void GPPhotonProcess::SetMaxKinEnergy(G4double e)
{
  nLambdaBins = GPlrint(nLambdaBins*log(e/minKinEnergy)
			/log(maxKinEnergy/minKinEnergy));
  maxKinEnergy = e;
}

//move to util
FQUALIFIER 
int GPPhotonProcess::GPlrint(double ad)
{
  return (ad>0) ? (G4int)(ad+.5) : (G4int)(ad-.5);
}

// inline

FQUALIFIER 
void GPPhotonProcess::DefineMaterial(GPMaterial* material)
{
  currentMaterial = material;
  fFactor = 1.0;
  mfpKinEnergy = DBL_MAX;
}

FQUALIFIER G4double GPPhotonProcess::GetLambdaFromTable(G4double e)
{
  //  return ((*theLambdaTable)[basedCoupleIndex])->Value(e);
  return theLambdaTable->physicsVectors[1].Value(e);
}

FQUALIFIER G4double GPPhotonProcess::ComputeCurrentLambda(G4double e)
{
  // return currentModel->CrossSectionPerVolume(baseMaterial,currentParticle,
  //                                           e,(*theCuts)[currentCoupleIndex]);
  //@@@ define theCuts value 
  return currentModel->CrossSectionPerVolume(currentMaterial,e,1.0*mm);
}

FQUALIFIER 
G4double GPPhotonProcess::GetCurrentLambda(G4double e)
{
  G4double x;
  if(theLambdaTable)   { x = GetLambdaFromTable(e); }
  else                 { x = ComputeCurrentLambda(e); }
  return fFactor*x;
}

FQUALIFIER
G4double GPPhotonProcess::GetLambda(G4double kinEnergy, 
				 GPMaterial* material)
{
  DefineMaterial(material);
  return GetCurrentLambda(kinEnergy);
}

FQUALIFIER 
G4double GPPhotonProcess::MinKinEnergy()
{
  return minKinEnergy;
}

FQUALIFIER 
G4double GPPhotonProcess::MaxKinEnergy()
{
  return maxKinEnergy;
}

FQUALIFIER
void GPPhotonProcess::SetPolarAngleLimit(G4double val)
{
  if(val < 0.0)     { polarAngleLimit = 0.0; }
  else if(val > pi) { polarAngleLimit = pi;  }
  else              { polarAngleLimit = val; }
}

FQUALIFIER 
G4double GPPhotonProcess::PolarAngleLimit()
{
  return polarAngleLimit;
}

FQUALIFIER 
void GPPhotonProcess::SetLambdaFactor(G4double val)
{
  if(val > 0.0 && val <= 1.0) { lambdaFactor = val; }
}

//---------------------------------------------------------------------------
//
// G4VProcess
//
//---------------------------------------------------------------------------

FQUALIFIER
void GPPhotonProcess::ResetNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft = -log(GPUniformRand(fDevStates,fThreadId));
  theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 
}

FQUALIFIER
void GPPhotonProcess::EndTracking()
{
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = -1.0;
  theInitialNumberOfInteractionLength=-1.0;
}
