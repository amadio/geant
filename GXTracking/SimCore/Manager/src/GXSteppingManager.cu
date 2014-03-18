#include "GXSteppingManager.h"

FQUALIFIER GXSteppingManager::GXSteppingManager(GPGeomManager *geomManager,
                                                GXFieldMap *magMap) 
{
  fStepStatus = fUndefined;
  fTrack = 0;
  fStep = 0;
  fMaterial = 0;

  fCurrentProcess = 0;
  fPhotonProcessVector = 0;
  fBrem = 0;
  fIoni = 0;
  fMsc = 0;

  physIntLength = DBL_MAX; 

  PhotonMAXofPostStepLoops = 0;
  PhotonMAXofAlongStepLoops = 0;
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
void GXSteppingManager::SetTransportation(GXTransportation* aTransportation,
					  GPNavigator* aNavigator)
{
  fNavigator = aNavigator;
  fTransportationProcess = aTransportation;
}

FQUALIFIER
GPStepStatus GXSteppingManager::Stepping()
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
  
  // Reset the secondary particles
  fN2ndariesPostStepDoIt = 0;
  
  //JA Set the volume before it is used (in DefineStepLength() for User Limit) 
  //  fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
  fCurrentVolume = GPNavigator_LocateGlobalPointAndSetup(fNavigator,
                   GPThreeVector_create(fTrack->x,fTrack->y,fTrack->z),
                                             NULL,false,true);
  fMaterial = GPLogicalVolume_GetMaterial(
              GPVPhysicalVolume_GetLogicalVolume(fCurrentVolume));

  //-----------------
  // AtRest Processes
  //-----------------

  //no AtRest Processes for Electrons and Photons

  //---------------------------------
  // AlongStep and PostStep Processes
  //---------------------------------
  
  // Find minimum Step length demanded by active disc./cont. processes
  DefinePhysicalStepLength();

  // Store the Step length (geometrical length) to G4Step and G4Track

  fStep->SetStepLength( PhysicalStep );
  //  fTrack->SetStepLength( PhysicalStep );
  fTrack->s =  PhysicalStep;
  G4double GeomStepLength = PhysicalStep;

  // Store StepStatus to PostStepPoint
  (fStep->GetPostStepPoint()).SetStepStatus( fStepStatus );

  // Invoke AlongStepDoIt 
  InvokeAlongStepDoItProcs();

  // Update track by taking into account all changes by AlongStepDoIt
  //@@@G4FWP is this redundant?
  fStep->UpdateTrack();

  // Update safety after invocation of all AlongStepDoIts
  //  endpointSafOrigin= fPostStepPoint->GetPosition();
  endpointSafOrigin= (fStep->GetPostStepPoint()).GetPosition();
  endpointSafety=  fmax( proposedSafety - GeomStepLength, kCarTolerance);
  
  (fStep->GetPostStepPoint()).SetSafety( endpointSafety );
  
  // Invoke PostStepDoIt
  InvokePostStepDoItProcs();

  //-------
  // Finale
  //-------
  
  // Update 'TrackLength' and remeber the Step length of the current Step
  //  fTrack->AddTrackLength(fStep->GetStepLength());
  fPreviousStepSize = fStep->GetStepLength();
  fStep->SetTrack(fTrack);

  //@@@G4FWP - skipping hit making for now
  /*
  // Send G4Step information to Hit/Dig if the volume is sensitive
  fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
  StepControlFlag =  fStep->GetControlFlag();
  if( fCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation) {
    fSensitive = fStep->GetPreStepPoint()->
      GetSensitiveDetector();
    if( fSensitive != 0 ) {
      fSensitive->Hit(fStep);
    }
  }
  */
  
  // Stepping process finish. Return the value of the StepStatus.
  return fStepStatus;
  
}

GPStepStatus GXSteppingManager::SteppingGPIL()
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
GPStepStatus GXSteppingManager::SteppingDoIt()
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
  endpointSafety=  fmax( proposedSafety - GeomStepLength, kCarTolerance);
  
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
void GXSteppingManager::GetProcessNumber(GXProcessManager* aProcessManager)
{
  PhotonMAXofPostStepLoops = aProcessManager->GetNumberOfPhotonProcess(); 
  if(PhotonMAXofPostStepLoops > 0) {
    fPhotonProcessVector = aProcessManager->GetPhotonProcessVector();
    //photon has no AlongStepGPIL
    PhotonMAXofAlongStepLoops = 0;
  }

  ElectronMAXofPostStepLoops = aProcessManager->GetNumberOfElectronProcess(); 
  if(ElectronMAXofPostStepLoops>0) {
    fBrem = aProcessManager->GetBremsstrahlungProcess();
    fIoni = aProcessManager->GetIonisationProcess();
    fMsc = aProcessManager->GetMultipleScatteringProcess();
    ElectronMAXofAlongStepLoops = ElectronMAXofPostStepLoops;
  }
}

FQUALIFIER
void GXSteppingManager::DefinePhysicalStepLength()
{
  PhysicalStep  = DBL_MAX;
  physIntLength = DBL_MAX;

  //Do Photon PostStepGPIL

  fPostStepDoItProcTriggered = 0;
  for(G4int np = 0 ; np <  PhotonMAXofPostStepLoops ; ++np) {
    fCurrentProcess = &(fPhotonProcessVector[np]);

    physIntLength = 
      fCurrentProcess->PostStepGetPhysicalInteractionLength(fTrack,
							   fMaterial,
							   fPreviousStepSize,
							   &fCondition );
    if(physIntLength < PhysicalStep ){
      PhysicalStep = physIntLength;
      fStepStatus = fPostStepDoItProc;
      fPostStepDoItProcTriggered = G4int(np);
    }
  }

  //Do Electron PostStepGPIL

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
      }
    }

    //multiple scattering return DBL_MAX, so never be the winner

  }

  // AlongStepGPIL  

  proposedSafety = DBL_MAX;
  G4double safetyProposedToAndByProcess = proposedSafety;

  for(size_t kp=0; kp < PhotonMAXofAlongStepLoops; kp++){
    //All photon processes considered are discrete processes,
    //so, do nothing for AlongStepGPIL for photon processes
    ;
  }

  // AlongStepGPIL for electron processes

  // do nothing for the brem process since it returns DBL_MAX
  if(fIoni) {
    physIntLength = 
      fIoni->AlongStepGetPhysicalInteractionLength(&fGPILSelection);
    
    if(physIntLength < PhysicalStep){
      PhysicalStep = physIntLength;
    }
  }

  if(fMsc) {
    physIntLength = 
      fMsc->AlongStepGetPhysicalInteractionLength(fMaterial,fTrack->E,
						  PhysicalStep,
						  &fGPILSelection);
    if(physIntLength < PhysicalStep){
      PhysicalStep = physIntLength;

      // Check if the process wants to be the GPIL winner. For example,
      // multi-scattering proposes Step limit, but won't be the winner.
      if(fGPILSelection==CandidateForSelection){
	fStepStatus = fAlongStepDoItProc;
      }
    }
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
  else {
    physIntLength = GXTransportation_AlongStepGPIL_Photon(
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
void GXSteppingManager::DefinePhysicalStepLength2()
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
void GXSteppingManager::SetInitialStep(GXTrack* valueTrack)
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
void GXSteppingManager::SetGPILStep(GXTrack* valueTrack)
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
void GXSteppingManager::SetDoItStep(GXTrack* valueTrack)
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

void GXSteppingManager::InvokeAlongStepDoItProcs()
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
    fParticleChange = fBrem->AlongStepDoIt(fTrack,fMaterial,PhysicalStep);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  if(fIoni) {
    fParticleChange = fIoni->AlongStepDoIt(fTrack,fMaterial,PhysicalStep);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  if(fMsc) {
    fParticleChange = fMsc->AlongStepDoIt(fMaterial,fTrack);
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

void GXSteppingManager::InvokeAlongStepDoItProcs2()
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
    fParticleChange = fBrem->AlongStepDoIt(fTrack,fMaterial,PhysicalStep);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  if(fIoni) {
    fParticleChange = fIoni->AlongStepDoIt(fTrack,fMaterial,PhysicalStep);
    fParticleChange.UpdateStepForAlongStep(fStep);
    fTrack->status = fParticleChange.GetTrackStatus();
    fParticleChange.Clear();
  }
  if(fMsc) {

    physIntLength = 
      fMsc->AlongStepGetPhysicalInteractionLength(fMaterial,fTrack->E,
                                                  PhysicalStep,
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

    fParticleChange = fMsc->AlongStepDoIt(fMaterial,fTrack);
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
void GXSteppingManager::InvokePostStepDoItProcs()
{
  // Invoke the specified discrete processes
  
  if(PhotonMAXofPostStepLoops > 0) {
    InvokePSDIP(fPostStepDoItProcTriggered);
  }  
  if(ElectronMAXofPostStepLoops > 0) {
    switch (fPostStepDoItProcTriggered) {
    case kBremsstrahlung : 
      InvokePSDIPeBremsstrahlung();
      break;
    case kIonisation :
      InvokePSDIPeIonisation();
      break;
    case kMultipleScattering :
      InvokePSDIPeMultipleScattering();
      break;
    default :
      break;
    }
  }
}

FQUALIFIER
void GXSteppingManager::InvokePostStepDoItProcs2(G4int proc)
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
void GXSteppingManager::InvokePSDIP(size_t np)
{
  //  fCurrentProcess = (*fPostStepDoItVector)[np];
  fCurrentProcess = &(fPhotonProcessVector[np]);

  fParticleChange 
    = fCurrentProcess->PostStepDoIt( fTrack, fMaterial);
  //    = fCurrentProcess->PostStepDoIt( *fTrack, *fStep);
  
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

    /*
    if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
      { ApplyProductionCut(tempSecondaryTrack); }
    
    // Set parentID 
    tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );
    
    // Set the process pointer which created this track 
    tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );
    */

    // If this 2ndry particle has 'zero' kinetic energy, make sure
    // it invokes a rest process at the beginning of the tracking
    if(tempSecondaryTrack->E <= DBL_MIN){
      ;
      //delete tempSecondaryTrack
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
void GXSteppingManager::InvokePSDIPeBremsstrahlung()
{
  fParticleChange = fBrem->PostStepDoIt( fTrack, fMaterial);
  
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
void GXSteppingManager::InvokePSDIPeIonisation() 
{
  fParticleChange = fIoni->PostStepDoIt( fTrack, fMaterial);
  
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
void GXSteppingManager::InvokePSDIPeMultipleScattering() 
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
void GXSteppingManager::SetSecondaryStack(GXTrack *secTracks,
					  G4int *stackSize) 
{
  theSecondaryStack = secTracks;
  theStackSize = stackSize;
}

FQUALIFIER
void GXSteppingManager::StoreSecondary(GXTrack* aTrack) 
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
void GXSteppingManager::SetStep(GPStep* aStep)
{
  fStep = aStep;
}

FQUALIFIER
GPStep* GXSteppingManager::GetStep()
{
  return fStep;
}

FQUALIFIER
void GXSteppingManager::SetMaterial(GPMaterial* aMaterial)
{
  fMaterial = aMaterial;
}

FQUALIFIER
GPMaterial* GXSteppingManager::GetMaterial()
{
  return fMaterial;
}

FQUALIFIER
G4double GXSteppingManager::CalculateSafety()
{
  return fmax( endpointSafety -
	      GPThreeVector_mag(
	      GPThreeVector_sub(endpointSafOrigin,
				(fStep->GetPostStepPoint()).GetPosition())),
	      kCarTolerance );
}
