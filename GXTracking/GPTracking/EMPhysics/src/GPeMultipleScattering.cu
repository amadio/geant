#include "GPeMultipleScattering.h"
#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"

//#include "G4LossTableManager.hh"
//#include "G4MaterialCutsCouple.hh"
//#include "G4Step.hh"
//#include "G4ParticleDefinition.hh"
//#include "G4VEmFluctuationModel.hh"
//#include "G4UnitsTable.hh"
//#include "G4ProductionCutsTable.hh"
//#include "G4Electron.hh"
//#include "G4GenericIon.hh"
//#include "G4TransportationManager.hh"
//#include "G4SafetyHelper.hh"

#include "stdio.h"

FQUALIFIER
GPeMultipleScattering::GPeMultipleScattering(curandState* devStates,
					     int threadId,
					     GPPhysicsTable* lambdaTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  geomMin = 1.e-6*mm;
  lowestKinEnergy = 1*eV;

  physStepLimit = gPathLength = tPathLength = 0.0;

  fPositionChanged = false;
  isActive = false;

  // G4VProcess
  theNumberOfInteractionLengthLeft = -1;
  currentInteractionLength = -1;
  theInitialNumberOfInteractionLength = -1.0;

  //model specific
  currentModel = 0;
  theLambdaTable = lambdaTable;

  //G4eMultipleScattering
  isInitialized = false;

}

FQUALIFIER
GPeMultipleScattering::~GPeMultipleScattering()
{
}

FQUALIFIER
void GPeMultipleScattering::InitialiseProcess(GPUrbanMscModel95* model)
{
  if(isInitialized) { return; }

  //SetEmModel
  currentModel = model;

  //set ParticleChange
  currentModel->SetParticleChange(&fParticleChange);

  currentModel->SetLambdaTable(theLambdaTable);
  //  currentModel->SetHighEnergyLimit(100*MeV);

  isInitialized = true;
}

//---------------------------------------------------------------------------
//
// G4VMultipleScattering
//
//---------------------------------------------------------------------------

FQUALIFIER
void GPeMultipleScattering::StartTracking()
{
  currentModel->StartTracking();
}

FQUALIFIER
void GPeMultipleScattering::EndTracking()
{
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = -1.0;
  theInitialNumberOfInteractionLength=-1.0;
}


FQUALIFIER
G4double GPeMultipleScattering::AlongStepGetPhysicalInteractionLength(
                                                   GPMaterial* material,
						   G4double kineticEnergy,
						   G4double currentMinimalStep,
						   GPGPILSelection* selection)
{
  // get Step limit proposed by the process
  *selection = NotCandidateForSelection;
  physStepLimit = gPathLength = tPathLength = currentMinimalStep;

  G4double ekin = kineticEnergy; //track.GetKineticEnergy();

  // select new model
  //  if(1 < numberOfModels) {
  //    currentModel = static_cast<G4VMscModel*>(
  //      SelectModel(ekin,track.GetMaterialCutsCouple()->GetIndex()));
  //  }

  // step limit
  //  if(currentModel->IsActive(ekin) && gPathLength >= geomMin 
  if(gPathLength >= geomMin && ekin >= lowestKinEnergy) {
    isActive = true;

    //tPathLength=currentModel->ComputeTruePathLengthLimit(track, gPathLength);
    G4double xlambda = GetTransportMeanFreePath(ekin);

    tPathLength = currentModel->ComputeTruePathLengthLimit(material,
					   ekin,xlambda,gPathLength);

    if (tPathLength < physStepLimit) { 
      *selection = CandidateForSelection; 
    }
  } 
  else { isActive = false; }

  return gPathLength;
}

FQUALIFIER 
G4double 
GPeMultipleScattering::PostStepGetPhysicalInteractionLength(GXTrack* track, 
					       GPForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

FQUALIFIER
GPVParticleChange&  GPeMultipleScattering::AlongStepDoIt(GPMaterial* material,
							 GXTrack* track)
{
  //  fParticleChange.ProposeMomentumDirection(
  //                  step.GetPostStepPoint()->GetMomentumDirection());
  //  fNewPosition = step.GetPostStepPoint()->GetPosition();
  GPThreeVector direction = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));
  fNewPosition = GPThreeVector_create(track->x,track->y,track->z);

  fParticleChange.SetProposedMomentumDirection(direction);
  fParticleChange.SetProposedPosition(fNewPosition);

  fPositionChanged = false;

  //  G4double geomLength = step.GetStepLength();
  G4double geomLength = track->s;

  // very small step - no msc
  if(!isActive) {
    tPathLength = geomLength;

    // sample msc
  } else {
    G4double range = currentModel->GetRange(material, track->E );
    //    currentModel->GetRange(currParticle,track.GetKineticEnergy(),
    //			   track.GetMaterialCutsCouple());
    G4double trueLength = currentModel->ComputeTrueStepLength(geomLength);
    
    // protection against wrong t->g->t conversion
    //    if(trueLength > tPathLength) 

    if (trueLength <= physStepLimit) {
      tPathLength = trueLength; 
    } else {
      tPathLength = physStepLimit - 0.5*geomMin; 
    }

    // do not sample scattering at the last or at a small step
    if(tPathLength + geomMin < range && tPathLength > geomMin) {

      //      G4double preSafety = step.GetPreStepPoint()->GetSafety();
      //      G4double postSafety= preSafety - geomLength; 
      //      G4bool safetyRecomputed = false;
      //      if( postSafety < geomMin ) {
      //	safetyRecomputed = true;
      //@@@	postSafety = currentModel->ComputeSafety(fNewPosition,0.0); 
      //      } 

      GPThreeVector displacement = 
      	currentModel->SampleScattering(material,direction,0.1*mm);
	// currentModel->SampleScattering(step.GetPostStepPoint()->GetMomentumDirection(),postSafety);

      G4double r2 = GPThreeVector_mag2(displacement);

      // make correction for displacement
      if(r2 > 0.0) {

	fPositionChanged = true;
        G4double fac = 1.0;

	// displaced point is definitely within the volume

	//@@@ check this part
	//	if(r2 > postSafety*postSafety) {
	//          if(!safetyRecomputed) {
	//	    postSafety = currentModel->ComputeSafety(fNewPosition, 0.0);
	//	  } 
	  // add a factor which ensure numerical stability
	//	  if(r2 > postSafety*postSafety) { 
	//	    fac = 0.99*postSafety/sqrt(r2); }
	//	}


	// compute new endpoint of the Step
	//	fNewPosition += fac*displacement;
	fNewPosition = GPThreeVector_add(fNewPosition,
					 GPThreeVector_mult(displacement,fac));
	//safetyHelper->ReLocateWithinVolume(fNewPosition);
      }

    }
  }

  fParticleChange.ProposeTrueStepLength(tPathLength);
  fParticleChange.SetProposedPosition(fNewPosition);
  return fParticleChange;

}

FQUALIFIER
GPVParticleChange& GPeMultipleScattering::PostStepDoIt(GXTrack* track)
{
  fParticleChange.ParticleChangeForMSC_Initialize(track);  
  
  if(fPositionChanged) { 
    //    safetyHelper->ReLocateWithinVolume(fNewPosition);
    fParticleChange.SetProposedPosition(fNewPosition); 
  }
  
  return fParticleChange;
}

FQUALIFIER
G4double GPeMultipleScattering::GetTransportMeanFreePath(G4double ekin)
{
  G4double x = theLambdaTable->physicsVectors[1].Value(ekin)*1.0/(ekin*ekin);
  if(0.0 >= x) { x = DBL_MAX; }
  else { x = 1.0/x; }
  return x;
}
