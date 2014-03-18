#include "GXeSteppingManager.h"

#include "GXeProcessManager.h"
#include "GXeBremsstrahlung.h"
#include "GXeIonisation.h"
#include "GXeMultipleScattering.h"

FQUALIFIER GXeSteppingManager::GXeSteppingManager(GPGeomManager *geomManager,
						  GXFieldMap *magMap) 
{
  fStepStatus = fUndefined;

  fTrack = 0;
  fStep = 0;
  fMaterial = 0;

  fBrem = 0;
  fIoni = 0;
  fMsc = 0;

  physIntLength = DBL_MAX; 

  ElectronMAXofPostStepLoops = 0;
  ElectronMAXofAlongStepLoops = 0;

  //  fAlongStepDoItProcTriggered = 0;
  fPostStepDoItProcTriggered = 0;

  //interface to the secondary stack on the global(GPG) or heap(CPU) memory
  theSecondaryStack =0;
  theStackSize = 0;

  //Navigator and Transportation process
  fNavigator = 0;
  fTransportationProcess = 0;

}

FQUALIFIER
void GXeSteppingManager::SetTransportation(GXTransportation* aTransportation,
					   GPNavigator* aNavigator)
{
  fNavigator = aNavigator;
  fTransportationProcess = aTransportation;
}

GPStepStatus GXeSteppingManager::SteppingGPIL()
{
  //--------
  // Prelude
  //--------
  // Store last PostStepPoint to PreStepPoint, and swap current and nex
  // volume information of G4Track. Reset total energy deposit in one Step. 
  fStep->CopyPostToPreStepPoint();
  fStep->ResetTotalEnergyDeposit();
  
  // Switch next touchable in track to current one
  //  fTrack->SetTouchableHandle(fTrack->GetNextTouchableHandle());
  
  //JA Set the volume before it is used (in DefineStepLength() for User Limit) 
  //  fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
  fCurrentVolume = GPNavigator_LocateGlobalPointAndSetup(fNavigator,
                   GPThreeVector_create(fTrack->x,fTrack->y,fTrack->z),
                                        NULL,false,true);
  fMaterial = GPLogicalVolume_GetMaterial(
              GPVPhysicalVolume_GetLogicalVolume(fCurrentVolume));

  //  fLiason->material = fMaterial;
  //---------------------------------
  // AlongStep and PostStep Processes
  //---------------------------------
  
  // Find minimum Step length demanded by active disc./cont. processes
  DefinePhysicalStepLength2();

  //  fTrack->SetStepLength( PhysicalStep );
  fTrack->s =  PhysicalStep;

  //temporary tweak to pass the step status to the SteppingDoIt step
  fTrack->status =  fStepStatus;

  // Stepping process finish. Return the value of the StepStatus.
  return fStepStatus;
  
}

FQUALIFIER
GPStepStatus GXeSteppingManager::SteppingDoIt()
{
  //--------
  // Prelude
  //--------
  // Store last PostStepPoint to PreStepPoint, and swap current and nex
  // volume information of G4Track. Reset total energy deposit in one Step. 
  //  fStep->CopyPostToPreStepPoint();
  //  fStep->ResetTotalEnergyDeposit();
  
  // Switch next touchable in track to current one
  //  fTrack->SetTouchableHandle(fTrack->GetNextTouchableHandle());
  
  //JA Set the volume before it is used (in DefineStepLength() for User Limit) 
  //  fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();

  //continute ater SteppingGPIL
  PhysicalStep = fTrack->s;
  fStepStatus = static_cast<GPStepStatus>(fTrack->status);

  fStep->SetStepLength( PhysicalStep );
  G4double GeomStepLength = PhysicalStep;

  // Store StepStatus to PostStepPoint
  (fStep->GetPostStepPoint()).SetStepStatus( fStepStatus );

  // Invoke AlongStepDoIt 
  InvokeAlongStepDoItProcs2();

  // Update track by taking into account all changes by AlongStepDoIt
  //@@@G4FWP is this redundant?
  fStep->UpdateTrack();

  // Update safety after invocation of all AlongStepDoIts
  //  endpointSafOrigin= fPostStepPoint->GetPosition();
  endpointSafOrigin= (fStep->GetPostStepPoint()).GetPosition();
  endpointSafety=  max( proposedSafety - GeomStepLength, kCarTolerance);
  
  (fStep->GetPostStepPoint()).SetSafety( endpointSafety );
  
  // Invoke PostStepDoIt
  InvokePostStepDoItProcs2(fTrack->proc);

  //-------
  // Finale
  //-------
  
  // Update 'TrackLength' and remeber the Step length of the current Step
  //  fTrack->AddTrackLength(fStep->GetStepLength());
  fPreviousStepSize = fStep->GetStepLength();
  fStep->SetTrack(fTrack);

  // Stepping process finish. Return the value of the StepStatus.
  return fStepStatus;  
}

FQUALIFIER
void GXeSteppingManager::GetProcessNumber(GXeProcessManager* aProcessManager)
{
  ElectronMAXofPostStepLoops = aProcessManager->GetNumberOfElectronProcess(); 
  if(ElectronMAXofPostStepLoops>0) {
    fBrem = aProcessManager->GetBremsstrahlungProcess();
    fIoni = aProcessManager->GetIonisationProcess();
    fMsc = aProcessManager->GetMultipleScatteringProcess();
    ElectronMAXofAlongStepLoops = ElectronMAXofPostStepLoops;
  }
}

FQUALIFIER
void GXeSteppingManager::DefinePhysicalStepLength2()
{
  PhysicalStep  = DBL_MAX;
  physIntLength = DBL_MAX;

  fPostStepDoItProcTriggered = 0;

  //Do Electron PostStepGPIL

  fTrack->proc = -1;

  if( ElectronMAXofPostStepLoops > 0) {

    if(fBrem) {
      physIntLength = 
        fBrem->PostStepGetPhysicalInteractionLength(fTrack,
                                                    fPreviousStepSize,
                                                    &fCondition);
      if(physIntLength < PhysicalStep ){
        PhysicalStep = physIntLength;
        fStepStatus = fPostStepDoItProc;
        fPostStepDoItProcTriggered = kBremsstrahlung;

        //temp word
        fTrack->proc = 0;
      }
    }    

    if(fIoni) {
      physIntLength = 
        fIoni->PostStepGetPhysicalInteractionLength(fTrack,
                                                    fPreviousStepSize,
                                                    &fCondition);
      if(physIntLength < PhysicalStep ){
        PhysicalStep = physIntLength;
        fStepStatus = fPostStepDoItProc;
        fPostStepDoItProcTriggered = kIonisation;

        //temp word
        fTrack->proc = 1;
      }
    }

    //multiple scattering return DBL_MAX, so never be the winner

  }

  // AlongStepGPIL  
  proposedSafety = DBL_MAX;
  G4double safetyProposedToAndByProcess = proposedSafety;

  // AlongStepGPIL for electron processes

  // do nothing for the brem process since it returns DBL_MAX

  if(fIoni) {
    physIntLength = 
      fIoni->AlongStepGetPhysicalInteractionLength(&fGPILSelection);
    
    if(physIntLength < PhysicalStep){
      PhysicalStep = physIntLength;
      fTrack->proc = 1;
    }
  }

  if(fMsc) {
    //move it to DoIt
  }

  // Transportation is assumed to be the last process in the vector

  if(ElectronMAXofPostStepLoops>0) {
    physIntLength = GXTransportation_AlongStepGPIL_Electron(
                                     fTransportationProcess,
                                     fTrack,
                                     fPreviousStepSize,
                                     PhysicalStep,
                                     safetyProposedToAndByProcess,
                                     &fGPILSelection);
  }

  if(physIntLength < PhysicalStep){
    PhysicalStep = physIntLength;
    fStepStatus = fGeomBoundary;
  }

  // Make sure to check the safety, even if Step is not limited 
  // by this process. J. Apostolakis, June 20, 1998
  // 

  if (safetyProposedToAndByProcess < proposedSafety) {
    // proposedSafety keeps the smallest value:
    proposedSafety               = safetyProposedToAndByProcess;
  }
  else {
    // safetyProposedToAndByProcess always proposes a valid safety:
    safetyProposedToAndByProcess = proposedSafety;
  } 

}

FQUALIFIER
void GXeSteppingManager::SetInitialStep(GXTrack* valueTrack)
{
  // Set up several local variables.
  fPreviousStepSize = 0.;
  fStepStatus = fUndefined;

  fTrack = valueTrack;

  PhysicalStep = 0.;
  GeometricalStep = 0.;
  CorrectedStep = 0.;

  if(fTrack->E <= 0.0){
    fTrack->status =  fStopButAlive;
  }

  //@@@G4FWP Set Touchable to track here if necessary

  // Initial set up for attribues of 'Step'
  fStep->InitializeStep(valueTrack);

  fTrack->s = 0.;

}

FQUALIFIER
void GXeSteppingManager::SetGPILStep(GXTrack* valueTrack)
{
  // Set up several local variables.
  fPreviousStepSize = 0.;
  fStepStatus = fUndefined;

  fTrack = valueTrack;
  //  fLiason = valueLiason;

  PhysicalStep = 0.;
  GeometricalStep = 0.;
  CorrectedStep = 0.;

  if(fTrack->E <= 0.0){
    fTrack->status =  fStopButAlive;
  }

  //@@@G4FWP Set Touchable to track here if necessary

  // Initial set up for attribues of 'Step'
  fStep->InitializeStep(valueTrack);

  fTrack->s = 0.;

}

FQUALIFIER
void GXeSteppingManager::SetDoItStep(GXTrack* valueTrack)
{
  // Set up several local variables.
  fPreviousStepSize = 0.;

  fTrack = valueTrack;
  //  fLiason = valueLiason;

  //passed from SteppingGPIL
  fStepStatus = static_cast<GPStepStatus>(fTrack->status);
  PhysicalStep = fTrack->s;

  GeometricalStep = 0.;
  CorrectedStep = 0.;

  //do this one more time because fTrack->status was set to fStepStatus
  //for the tweak used in SteppingGPIL
  if(fTrack->E <= 0.0){
    fTrack->status =  fStopButAlive;
  }

  //@@@G4FWP Set Touchable to track here if necessary

  // Initial set up again for attribues of 'Step'
  fStep->InitializeStep(valueTrack);

}

void GXeSteppingManager::InvokeAlongStepDoItProcs2()
{
  // If the current Step is defined by a 'ExclusivelyForced' 
  // PostStepDoIt, then don't invoke any AlongStepDoIt
  if(fStepStatus == fExclusivelyForcedProc){
    return; 
  }
  
  // no AlongStepDoIt for photon processes
  // none of AlongStepDoIt of electron processes has secondaries,
  // so, the secondary handling is skipped

  if(fBrem) {
    fParticleChange = fBrem->AlongStepDoIt(fTrack,fMaterial);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  if(fIoni) {
    fParticleChange = fIoni->AlongStepDoIt(fTrack,fMaterial);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  if(fMsc) {

    physIntLength = 
      fMsc->AlongStepGetPhysicalInteractionLength(
                                                  &fGPILSelection);
    if(physIntLength < PhysicalStep){
      PhysicalStep = physIntLength;

      fTrack->s =  PhysicalStep;
      // Check if the process wants to be the GPIL winner. For example,
      // multi-scattering proposes Step limit, but won't be the winner.
      if(fGPILSelection==CandidateForSelection){
        fStepStatus = fAlongStepDoItProc;
        fTrack->proc = 2;
      }
    }

    fParticleChange = fMsc->AlongStepDoIt(fTrack,fMaterial);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  
  //@@@G4FWP insert AlongstepDoIt of the transporation process here

  fStep->UpdateTrack();

  if ( fTrack->status == fAlive && fTrack->E <= DBL_MIN ) {
    fTrack->status = fStopAndKill;
  }
}

FQUALIFIER
void GXeSteppingManager::InvokePostStepDoItProcs2(G4int proc)
{
  // Invoke the specified discrete processes
  if(ElectronMAXofPostStepLoops > 0) {
    switch (proc) {
    case 0 : 
      InvokePSDIPeBremsstrahlung();
      break;
    case 1 :
      InvokePSDIPeIonisation();
      break;
    case 2 :
      InvokePSDIPeMultipleScattering();
      break;
    default :
      break;
    }
  }
}

FQUALIFIER
void GXeSteppingManager::InvokePSDIPeBremsstrahlung()
{
  fParticleChange = fBrem->PostStepDoIt( fTrack);
  
  // Update PostStepPoint of Step according to ParticleChange
  fParticleChange.UpdateStepForPostStep(fStep);
  
  // Update G4Track according to ParticleChange after each PostStepDoIt
  fStep->UpdateTrack();
  
  // Update safety after each invocation of PostStepDoIts
  (fStep->GetPostStepPoint()).SetSafety( CalculateSafety() );
  
  // Now Store the secondaries from ParticleChange to SecondaryList
  GXTrack* tempSecondaryTrack;
  G4int    num2ndaries;
  
  num2ndaries = fParticleChange.GetNumberOfSecondaries();
  
  for(G4int DSecLoop=0 ; DSecLoop < num2ndaries; DSecLoop++){

    tempSecondaryTrack = fParticleChange.GetSecondary(DSecLoop);

    if(tempSecondaryTrack->E <= DBL_MIN){
      //do not store this secondary
      ;
    } else {
      StoreSecondary(tempSecondaryTrack);
      fN2ndariesPostStepDoIt++;
    }
  } //end of loop on secondary 
  
  // Set the track status according to what the process defined
  fTrack->status = fParticleChange.GetTrackStatus();
  
  // clear ParticleChange
  fParticleChange.Clear();

}

FQUALIFIER
void GXeSteppingManager::InvokePSDIPeIonisation() 
{
  fParticleChange = fIoni->PostStepDoIt( fTrack);
  
  // Update PostStepPoint of Step according to ParticleChange
  fParticleChange.UpdateStepForPostStep(fStep);
  
  // Update G4Track according to ParticleChange after each PostStepDoIt
  fStep->UpdateTrack();
  
  // Update safety after each invocation of PostStepDoIts
  (fStep->GetPostStepPoint()).SetSafety( CalculateSafety() );
  
  // Now Store the secondaries from ParticleChange to SecondaryList
  GXTrack* tempSecondaryTrack;
  G4int    num2ndaries;
  
  num2ndaries = fParticleChange.GetNumberOfSecondaries();
  
  for(G4int DSecLoop=0 ; DSecLoop < num2ndaries; DSecLoop++){

    tempSecondaryTrack = fParticleChange.GetSecondary(DSecLoop);

    if(tempSecondaryTrack->E <= DBL_MIN){
      //do not store this secondary
      ;
    } else {
      StoreSecondary(tempSecondaryTrack);
      fN2ndariesPostStepDoIt++;
    }
  } //end of loop on secondary 
  
  // Set the track status according to what the process defined
  fTrack->status = fParticleChange.GetTrackStatus();
  
  // clear ParticleChange
  fParticleChange.Clear();

}

FQUALIFIER
void GXeSteppingManager::InvokePSDIPeMultipleScattering() 
{
  fParticleChange = fMsc->PostStepDoIt( fTrack);
  
  // Update PostStepPoint of Step according to ParticleChange
  fParticleChange.UpdateStepForPostStep(fStep);
  
  // Update G4Track according to ParticleChange after each PostStepDoIt
  fStep->UpdateTrack();
  
  // Update safety after each invocation of PostStepDoIts
  (fStep->GetPostStepPoint()).SetSafety( CalculateSafety() );

  //no secondaries produced by the GPeMultipleScattering
  
  // Set the track status according to what the process defined
  fTrack->status = fParticleChange.GetTrackStatus();
  
  // clear ParticleChange
  fParticleChange.Clear();

}

FQUALIFIER
void GXeSteppingManager::SetSecondaryStack(GXTrack *secTracks,
					   G4int *stackSize) 
{
  theSecondaryStack = secTracks;
  theStackSize = stackSize;
}

FQUALIFIER
void GXeSteppingManager::StoreSecondary(GXTrack* aTrack) 
{
#ifndef __CUDA_ARCH__
  theSecondaryStack[*theStackSize] = *aTrack;
  ++(*theStackSize);
#else
  int theOffset = atomicAdd(theStackSize,1);
  theSecondaryStack[theOffset] = *aTrack;
#endif
}

FQUALIFIER
void GXeSteppingManager::SetStep(GPStep* aStep)
{
  fStep = aStep;
}

FQUALIFIER
GPStep* GXeSteppingManager::GetStep()
{
  return fStep;
}

FQUALIFIER
void GXeSteppingManager::SetMaterial(GPMaterial* aMaterial)
{
  fMaterial = aMaterial;
}

FQUALIFIER
GPMaterial* GXeSteppingManager::GetMaterial()
{
  return fMaterial;
}

FQUALIFIER
G4double GXeSteppingManager::CalculateSafety()
{
  return max( endpointSafety -
	      GPThreeVector_mag(
	      GPThreeVector_sub(endpointSafOrigin,
				(fStep->GetPostStepPoint()).GetPosition())),
	      kCarTolerance );
}
