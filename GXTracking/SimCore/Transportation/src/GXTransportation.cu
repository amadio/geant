#include "GXTransportation.h"
#include "GXFieldTrack.h"

#include "GXTrack.h"
#include "GPUtils.h"
#include "GXFieldTrack.h"

// Constructor

FQUALIFIER
void GXTransportation_Constructor( GXTransportation *This, 
                                   GXPropagatorInField *propagator,
				   G4int verboseLevel )
{
  This->fParticleIsLooping = false ;
  This->fPreviousSftOrigin = GPThreeVector_create(0.,0.,0.);
  This->fPreviousSafety = 0.0 ;
  This->fThreshold_Warning_Energy = 100 * MeV ;  
  This->fThreshold_Important_Energy = 250 * MeV ; 
  This->fThresholdTrials = 10 ; 
  This->fUnimportant_Energy = 1 * MeV ;  // Not used
  This->fNoLooperTrials = 0;
  This->fSumEnergyKilled = 0.0 ;
  This->fMaxEnergyKilled = 0.0 ;
  This->fShortStepOptimisation = false; // Old default: true (=fast short steps)
  This->fUseMagneticMoment = false;  

  This->fFieldPropagator = propagator;
  This->fLinearNavigator = 
    GXPropagatorInField_GetNavigatorForPropagating(propagator) ;

}

// 1. General

FQUALIFIER
G4double  GXTransportation_AlongStepGetPhysicalInteractionLength( 
            GXTransportation *This,
			   GXTrack   *track,
				G4double  previousStepSize,
				G4double  currentMinimumStep,
				G4double  currentSafety,
				GPGPILSelection* selection )
{
  G4double geometryStepLength= -1.0, newSafety= -1.0; 
  This->fParticleIsLooping = false ;
  
  // Initial actions moved to  StartTrack()   
  // --------------------------------------
  // Note: in case another process changes touchable handle
  //    it will be necessary to add here (for all steps)   
  // fCurrentTouchableHandle = aTrack->GetTouchableHandle();
  
  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  //
  *selection = CandidateForSelection ;

  // Get initial Energy/Momentum of the track
  //
  //  const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
  //  const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;

  //G4FWP On device, GXTrack is used instead of G4Track

  //  GPThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  //  GPThreeVector startPosition          = track.GetPosition() ;

  GPThreeVector startMomentum = GPThreeVector_create(track->px,
						     track->py,
						     track->pz);

  GPThreeVector startPosition    = GPThreeVector_create(track->x,
                                                        track->y,
                                                        track->z);

  // The Step Point safety can be limited by other geometries and/or the 
  // assumptions of any process - it's not always the geometrical safety.
  // We calculate the starting point's isotropic safety here.
  //
  GPThreeVector OriginShift = GPThreeVector_sub(startPosition, 
						This->fPreviousSftOrigin) ;
  G4double      MagSqShift  = GPThreeVector_mag2(OriginShift) ;
  if( MagSqShift >= GPsqr(This->fPreviousSafety) ) {
    currentSafety = 0.0 ;
  }
  else {
    currentSafety = This->fPreviousSafety - sqrt(MagSqShift) ;
  }
  
  // Is the particle charged or has it a magnetic moment?
  G4double particleCharge = track->q ; 
  G4double magneticMoment = 0.0 ;//@@@G4FWP temporary value

  G4double       restMass = 0.0 ;//@@@G4FWP temporary value
  G4double       E = 0.0 ;

  if(fabs(particleCharge)<0.5) { //photon
    E = GPThreeVector_mag(startMomentum); 
  }
  else { //e+/e-
    restMass = 0.510998910;
    E = sqrt(GPThreeVector_mag2(startMomentum) + restMass*restMass); 
  }
  
  This->fGeometryLimitedStep = false ;
  
  // Check if the particle has a force, EM or gravitational, exerted on it
  //
  GXFieldManager* fieldMgr = 0;
  G4bool          fieldExertsForce = false ;
  
  G4bool gravityOn = false;
  G4bool fieldExists= false;  // Field is not 0 (null pointer)
  
  //********************************************************************
  //@@@G4FWP use the magnetic field map on the device 
  //  fieldMgr = fFieldPropagator->FindAndSetFieldManager( track.GetVolume() );
  fieldMgr = GXPropagatorInField_FindAndSetFieldManager(This->fFieldPropagator);
  //@@@G4FWP************************************************************
  
  if( fieldMgr != 0 ) {
    // Message the field Manager, to configure it for this track
    GXFieldManager_ConfigureForTrack(fieldMgr
		                     //	&track 
				     );
    GXMagneticField* ptrField = GXFieldManager_GetDetectorField(fieldMgr);
    fieldExists = (ptrField!=0) ;
    
    if( fieldExists ) {
      gravityOn= false;
      
      if(  (particleCharge != 0.0) 
	   || (This->fUseMagneticMoment && (magneticMoment != 0.0) )
	   || (gravityOn          && (restMass != 0.0) )) {
	fieldExertsForce = fieldExists; 
      }
    }
  }
  
  if( !fieldExertsForce ) {
    G4double linearStepLength ;
    if(This->fShortStepOptimisation && (currentMinimumStep <= currentSafety) ) {
      // The Step is guaranteed to be taken
      //
      geometryStepLength   = currentMinimumStep ;
      This->fGeometryLimitedStep = false ;
    }
    else {
      //  Find whether the straight path intersects a volume
      //
      linearStepLength = GPNavigator_ComputeStep(This->fLinearNavigator, 
						 startPosition, 
						 startMomentum,
						 currentMinimumStep, 
						 &newSafety); 
      // Remember last safety origin & value.
      //
      This->fPreviousSftOrigin = startPosition ;
      This->fPreviousSafety    = newSafety ; 
      
      currentSafety = newSafety ;
      
      This->fGeometryLimitedStep= (linearStepLength <= currentMinimumStep); 
      if( This->fGeometryLimitedStep ) {
	// The geometry limits the Step size (an intersection was found.)
	geometryStepLength   = linearStepLength ;
      } 
      else {
	// The full Step is taken.
	geometryStepLength   = currentMinimumStep ;
      }
    }
    
    This->endpointDistance = geometryStepLength ;
    
    // Calculate final position
    //
    This->fTransportEndPosition = GPThreeVector_add(startPosition,
						    GPThreeVector_mult(GPThreeVector_unit(startMomentum),geometryStepLength)) ;
    
    // Momentum direction, energy and polarisation are unchanged by transport
    //
    //     This->fTransportEndMomentumDir   = startMomentumDir ; 
    //     This->fTransportEndKineticEnergy = track.GetKineticEnergy() ;
    //     This->fTransportEndSpin          = track.GetPolarization();
    This->fTransportEndKineticEnergy = E ;
    This->fParticleIsLooping         = false ;
    This->fMomentumChanged           = false ; 
  }
  else {  //  A field exerts force 
    //     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
    G4double       momentumMagnitude = E ; //G4FWP temporary
    //    GPThreeVector  EndUnitMomentum ;
    G4double       lengthAlongCurve ;
    
    GXPropagatorInField_SetChargeMomentumMass( This->fFieldPropagator,
    					       particleCharge,    // in e+ units
    					       momentumMagnitude, // in Mev/c 
    					       restMass           ) ;  
    
    GXFieldTrack  aFieldTrack;
    GXFieldTrack_Constructor(&aFieldTrack,
			     startPosition, 
			     startMomentum, //track.GetMomentumDirection(),
			     E, //track.GetKineticEnergy(),
			     restMass,
			     particleCharge,0.0);
    if( currentMinimumStep > 0 ) {
      // Do the Transport in the field (non recti-linear)
      //
      // @@@G4FWP 
      // get the physical volume of the current position and substitute 
      // track.GetVolume() in the original code
      
      GPVPhysicalVolume* cur_phy =
	GPNavigator_LocateGlobalPointAndSetup(This->fLinearNavigator,
					      startPosition,NULL,false,true);
      
      lengthAlongCurve = GXPropagatorInField_ComputeStep(
							 This->fFieldPropagator, 
							 aFieldTrack,
							 currentMinimumStep, 
							 currentSafety, 
							 cur_phy // track.GetVolume() //G4FWP 
							 ) ;
      This->fGeometryLimitedStep= lengthAlongCurve < currentMinimumStep; 
      if( This->fGeometryLimitedStep ) {
	geometryStepLength   = lengthAlongCurve ;
      } else {
	geometryStepLength   = currentMinimumStep ;
      }
      
      // Remember last safety origin & value.
      //
      This->fPreviousSftOrigin = startPosition ;
      This->fPreviousSafety    = currentSafety ;         
    }
    else {
      geometryStepLength   = lengthAlongCurve= 0.0 ;
      This->fGeometryLimitedStep = false ;
    }
    
    // Get the End-Position and End-Momentum (Dir-ection)
    //
    This->fTransportEndPosition = GXFieldTrack_GetPosition(&aFieldTrack) ;
    
    // Momentum:  Magnitude and direction can be changed too now ...
    //
    This->fMomentumChanged         = true ; 
    //     This->fTransportEndMomentumDir = GXFieldTrack_GetMomentumDir(&aFieldTrack);
    
    This->fTransportEndKineticEnergy = 
      GXFieldTrack_GetKineticEnergy(&aFieldTrack); 
    
    //    G4double  startEnergy= E; //G4FWP track.GetKineticEnergy();
    //    G4double  endEnergy= This->fTransportEndKineticEnergy; 
    
    This->fTransportEndKineticEnergy= E; // track.GetKineticEnergy(); 
  
    This->fParticleIsLooping = 
      GXPropagatorInField_IsParticleLooping(This->fFieldPropagator) ;
    This->endpointDistance = 
      GPThreeVector_mag(GPThreeVector_sub(This->fTransportEndPosition,
					  startPosition)) ;
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if( currentMinimumStep == 0.0 ) 
  {
      if( currentSafety == 0.0 )  This->fGeometryLimitedStep = true ;
  }

  // Update the safety starting from the end-point,
  // if it will become negative at the end-point.
  //
  if( currentSafety < This->endpointDistance ) 
  {
    if( particleCharge != 0.0 ) 
    {
      G4double endSafety =
            	GPNavigator_ComputeSafety(This->fLinearNavigator,
            				  This->fTransportEndPosition,
					  DBL_MAX, false) ;

      currentSafety      = endSafety ;
      This->fPreviousSftOrigin = This->fTransportEndPosition ;
      This->fPreviousSafety    = currentSafety ; 

      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength
      //
      currentSafety     += This->endpointDistance ;
      
    }
  }            

  //  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

//-----------------------------------------------------------------------------
// 2. Photon Only
//-----------------------------------------------------------------------------

FQUALIFIER
G4double  GXTransportation_AlongStepGPIL_Photon( 
                                GXTransportation *This,
			        GXTrack *track,
				G4double  previousStepSize,
				G4double  currentMinimumStep,
				G4double  currentSafety,
				GPGPILSelection* selection )
{
  G4double geometryStepLength= -1.0, newSafety= -1.0; 
  This->fParticleIsLooping = false ;
  
  *selection = CandidateForSelection ;

  GPThreeVector startMomentum = GPThreeVector_create(track->px,
						     track->py,
						     track->pz);

  GPThreeVector startPosition    = GPThreeVector_create(track->x,
                                                        track->y,
                                                        track->z);

  GPThreeVector OriginShift = GPThreeVector_sub(startPosition, 
						This->fPreviousSftOrigin) ;
  G4double      MagSqShift  = GPThreeVector_mag2(OriginShift) ;

  G4double  E = GPThreeVector_mag(startMomentum); 

  if( MagSqShift >= GPsqr(This->fPreviousSafety) ) {
    currentSafety = 0.0 ;
  }
  else {
    currentSafety = This->fPreviousSafety - sqrt(MagSqShift) ;
  }
  
  This->fGeometryLimitedStep = false ;
  
  G4double linearStepLength ;
  if(This->fShortStepOptimisation && (currentMinimumStep <= currentSafety) ) {
    // The Step is guaranteed to be taken
    //
    geometryStepLength   = currentMinimumStep ;
    This->fGeometryLimitedStep = false ;
  }
  else {
    //  Find whether the straight path intersects a volume
    //
    linearStepLength = GPNavigator_ComputeStep(This->fLinearNavigator, 
					       startPosition, 
					       startMomentum,
					       currentMinimumStep, 
					       &newSafety); 
    // Remember last safety origin & value.
    //
    This->fPreviousSftOrigin = startPosition ;
    This->fPreviousSafety    = newSafety ; 
    
    currentSafety = newSafety ;
    
    This->fGeometryLimitedStep= (linearStepLength <= currentMinimumStep); 
    if( This->fGeometryLimitedStep ) {
      // The geometry limits the Step size (an intersection was found.)
      geometryStepLength   = linearStepLength ;
    } 
    else {
      // The full Step is taken.
      geometryStepLength   = currentMinimumStep ;
    }
  }
  
  This->endpointDistance = geometryStepLength ;
  
  // Calculate final position
  //
  This->fTransportEndPosition = GPThreeVector_add(startPosition,
						  GPThreeVector_mult(GPThreeVector_unit(startMomentum),geometryStepLength)) ;
  
  This->fTransportEndKineticEnergy = E ;
  This->fParticleIsLooping         = false ;
  This->fMomentumChanged           = false ; 
  
  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if( currentMinimumStep == 0.0 ) 
    {
      if( currentSafety == 0.0 )  This->fGeometryLimitedStep = true ;
  }

  //  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

//-----------------------------------------------------------------------------
// 3. Electron Only
//-----------------------------------------------------------------------------

FQUALIFIER
G4double  GXTransportation_AlongStepGPIL_Electron( 
                                GXTransportation *This,
			        GXTrack *track,
				G4double  previousStepSize,
				G4double  currentMinimumStep,
				G4double  currentSafety,
				GPGPILSelection* selection )
{
  G4double geometryStepLength= -1.0;
  //G4double newSafety= -1.0; 
  This->fParticleIsLooping = false ;
  
  // Initial actions moved to  StartTrack()   
  // --------------------------------------
  // Note: in case another process changes touchable handle
  //    it will be necessary to add here (for all steps)   
  // fCurrentTouchableHandle = aTrack->GetTouchableHandle();
  
  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  //
  *selection = CandidateForSelection ;

  // Get initial Energy/Momentum of the track
  //
  //  const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
  //  const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;

  //G4FWP On device, GXTrack is used instead of G4Track

  //  GPThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  //  GPThreeVector startPosition          = track.GetPosition() ;

  GPThreeVector startMomentum = GPThreeVector_create(track->px,
						     track->py,
						     track->pz);

  GPThreeVector startPosition    = GPThreeVector_create(track->x,
                                                        track->y,
                                                        track->z);

  // The Step Point safety can be limited by other geometries and/or the 
  // assumptions of any process - it's not always the geometrical safety.
  // We calculate the starting point's isotropic safety here.
  //
  GPThreeVector OriginShift = GPThreeVector_sub(startPosition, 
						This->fPreviousSftOrigin) ;

  G4double      MagSqShift  = GPThreeVector_mag2(OriginShift) ;
  if( MagSqShift >= GPsqr(This->fPreviousSafety) ) {
    currentSafety = 0.0 ;
  }
  else {
    currentSafety = This->fPreviousSafety - sqrt(MagSqShift) ;
  }
  
  // Is the particle charged or has it a magnetic moment?
  G4double particleCharge = -1.0 ;  //Electron

  //  G4double magneticMoment = 0.0 ;//@@@G4FWP temporary value
  G4double       restMass =   0.51099891;//@@@G4FWP Electron
  G4double       E = sqrt(GPThreeVector_mag2(startMomentum) + restMass*restMass) ;
  
  This->fGeometryLimitedStep = false ;
  
  // Check if the particle has a force, EM or gravitational, exerted on it
  //
  GXFieldManager* fieldMgr = 0;
  //  G4bool          fieldExertsForce = false;
  
  //  G4bool gravityOn = false;
  //  G4bool fieldExists= false;  // Field is not 0 (null pointer)
  
  //********************************************************************
  //@@@G4FWP use the magnetic field map on the device 
  //  fieldMgr = fFieldPropagator->FindAndSetFieldManager( track.GetVolume() );
  fieldMgr = GXPropagatorInField_FindAndSetFieldManager(This->fFieldPropagator);
  //@@@G4FWP************************************************************
  
  if( fieldMgr != 0 ) {
    // Message the field Manager, to configure it for this track
    GXFieldManager_ConfigureForTrack(fieldMgr
		                     //	&track 
				     );
    /*
    GXMagneticField* ptrField = GXFieldManager_GetDetectorField(fieldMgr);
    fieldExists = (ptrField!=0) ;
    
    if( fieldExists ) {
    gravityOn= false;
      
      if(  (particleCharge != 0.0) 
	   || (This->fUseMagneticMoment && (magneticMoment != 0.0) )
	   || (gravityOn          && (restMass != 0.0) )) {
	fieldExertsForce = fieldExists; 
      }
    }
    */
  }
  
  //  A field exerts force 
  //     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
  G4double       momentumMagnitude = E ; //G4FWP temporary
  //    GPThreeVector  EndUnitMomentum ;
  G4double       lengthAlongCurve ;
  
  GXPropagatorInField_SetChargeMomentumMass( This->fFieldPropagator,
					     particleCharge,    // in e+ units
					     momentumMagnitude, // in Mev/c 
					     restMass           ) ;  
  
  GXFieldTrack  aFieldTrack;
  GXFieldTrack_Constructor(&aFieldTrack,
			   startPosition, 
			   startMomentum, //track.GetMomentumDirection(),
			   E, //track.GetKineticEnergy(),
			   restMass,
			   particleCharge,0.0);
  if( currentMinimumStep > 0 ) {
    // Do the Transport in the field (non recti-linear)
    //
    // @@@G4FWP 
    // get the physical volume of the current position and substitute 
    // track.GetVolume() in the original code
    
    GPVPhysicalVolume* cur_phy =
      GPNavigator_LocateGlobalPointAndSetup(This->fLinearNavigator,
					    startPosition,NULL,false,true);
    
    lengthAlongCurve = GXPropagatorInField_ComputeStep(
						       This->fFieldPropagator, 
						       aFieldTrack,
						       currentMinimumStep, 
						       currentSafety, 
						       cur_phy // track.GetVolume() //G4FWP 
						       ) ;
    This->fGeometryLimitedStep= lengthAlongCurve < currentMinimumStep; 
    if( This->fGeometryLimitedStep ) {
      geometryStepLength   = lengthAlongCurve ;
    } else {
      geometryStepLength   = currentMinimumStep ;
    }
    
    // Remember last safety origin & value.
    //
    This->fPreviousSftOrigin = startPosition ;
    This->fPreviousSafety    = currentSafety ;         
  }
  else {
    geometryStepLength   = lengthAlongCurve= 0.0 ;
    This->fGeometryLimitedStep = false ;
  }
  
  // Get the End-Position and End-Momentum (Dir-ection)
  //
  This->fTransportEndPosition = GXFieldTrack_GetPosition(&aFieldTrack) ;
  
  // Momentum:  Magnitude and direction can be changed too now ...
  //
  This->fMomentumChanged         = true ; 
  //     This->fTransportEndMomentumDir = GXFieldTrack_GetMomentumDir(&aFieldTrack);
  
  This->fTransportEndKineticEnergy = 
    GXFieldTrack_GetKineticEnergy(&aFieldTrack); 
  
  //  G4double  startEnergy= E; //G4FWP track.GetKineticEnergy();
  //  G4double  endEnergy= This->fTransportEndKineticEnergy; 
  
  This->fTransportEndKineticEnergy= E; // track.GetKineticEnergy(); 
  
  This->fParticleIsLooping = 
    GXPropagatorInField_IsParticleLooping(This->fFieldPropagator) ;
  This->endpointDistance = 
    GPThreeVector_mag(GPThreeVector_sub(This->fTransportEndPosition,
					startPosition)) ;

  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if( currentMinimumStep == 0.0 ) 
  {
      if( currentSafety == 0.0 )  This->fGeometryLimitedStep = true ;
  }

  // Update the safety starting from the end-point,
  // if it will become negative at the end-point.
  //
  if( currentSafety < This->endpointDistance ) 
  {
    if( particleCharge != 0.0 ) 
    {
      G4double endSafety =
            	GPNavigator_ComputeSafety(This->fLinearNavigator,
            				  This->fTransportEndPosition,
					  DBL_MAX, false) ;

      currentSafety      = endSafety ;
      This->fPreviousSftOrigin = This->fTransportEndPosition ;
      This->fPreviousSafety    = currentSafety ; 

      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength
      //
      currentSafety     += This->endpointDistance ;
      
    }
  }            

  //  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}
