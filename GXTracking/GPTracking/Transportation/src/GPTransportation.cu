//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4Transportation.cc 2011/06/10 16:19:46 japost Exp japost $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// ------------------------------------------------------------
//  GEANT 4  include file implementation
//
// ------------------------------------------------------------
//
// This class is a process responsible for the transportation of 
// a particle, ie the geometrical propagation that encounters the 
// geometrical sub-volumes of the detectors.
//
// It is also tasked with the key role of proposing the "isotropic safety",
//   which will be used to update the post-step point's safety.
//
// =======================================================================
// Modified:   
//   28 Oct  2011, P.Gumpl./J.Ap: Detect gravity field, use magnetic moment 
//   20 Nov  2008, J.Apostolakis: Push safety to helper - after ComputeSafety
//    9 Nov  2007, J.Apostolakis: Flag for short steps, push safety to helper
//   19 Jan  2006, P.MoraDeFreitas: Fix for suspended tracks (StartTracking)
//   11 Aug  2004, M.Asai: Add G4VSensitiveDetector* for updating stepPoint.
//   21 June 2003, J.Apostolakis: Calling field manager with 
//                     track, to enable it to configure its accuracy
//   13 May  2003, J.Apostolakis: Zero field areas now taken into
//                     account correclty in all cases (thanks to W Pokorski).
//   29 June 2001, J.Apostolakis, D.Cote-Ahern, P.Gumplinger: 
//                     correction for spin tracking   
//   20 Febr 2001, J.Apostolakis:  update for new FieldTrack
//   22 Sept 2000, V.Grichine:     update of Kinetic Energy
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "GPTransportation.h"
//#include "G4TransportationProcessType.hh"

//#include "G4ProductionCutsTable.hh"
//#include "G4ParticleTable.hh"
#include "GPChordFinder.h"
#include "GPFieldTrack.h"
//#include "G4SafetyHelper.hh"
//#include "G4FieldManagerStore.hh"

//class G4VSensitiveDetector;
#include "GPTrack.h"
#include "GPUtils.h"

//#include "stdio.h"

//////////////////////////////////////////////////////////////////////////
//
// Constructor

FQUALIFIER
void GPTransportation_Constructor( GPTransportation *This, 
                                   GPPropagatorInField *propagator,
				   G4int verboseLevel )
{
  //  G4VProcess( G4String("Transportation"), fTransportation ),

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
  This->fVerboseLevel = verboseLevel;

  // set Process Sub Type
  //  SetProcessSubType(static_cast<int>(TRANSPORTATION));

  //@@@G4FWP: G4TransportationManager is passed as an input
  //  G4TransportationManager* transportMgr ; 
  //  transportMgr = G4TransportationManager::GetTransportationManager() ; 

  //  fLinearNavigator = transportMgr->GetNavigatorForTracking() ; 

  This->fFieldPropagator = propagator;
  This->fLinearNavigator = 
    GPPropagatorInField_GetNavigatorForPropagating(propagator) ;

  //  printf("address of the navigator %ld\n",This->fLinearNavigator);

  //  fpSafetyHelper =   transportMgr->GetSafetyHelper();  // New 

  // Cannot determine whether a field exists here, as it would 
  //  depend on the relative order of creating the detector's 
  //  field and this process. That order is not guaranted.
  // Instead later the method DoesGlobalFieldExist() is called

  //  static G4TouchableHandle nullTouchableHandle;  // Points to (G4VTouchable*) 0
  //  fCurrentTouchableHandle = nullTouchableHandle; 

  This->fEndGlobalTimeComputed  = false;
  This->fCandidateEndGlobalTime = 0;

}

FQUALIFIER
void GPTransportation_Constructor2( GPTransportation *This, 
				    GPTransportationManager *transportMgr,
				    G4int verboseLevel )
{
  //  G4VProcess( G4String("Transportation"), fTransportation ),

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
  This->fVerboseLevel = verboseLevel;

  // set Process Sub Type
  //  SetProcessSubType(static_cast<int>(TRANSPORTATION));

  //@@@G4FWP: G4TransportationManager is passed as an input
  //  G4TransportationManager* transportMgr ; 
  //  transportMgr = G4TransportationManager::GetTransportationManager() ; 

  //  fLinearNavigator = transportMgr->GetNavigatorForTracking() ; 

  This->fFieldPropagator = 
    GPTransportationManager_GetPropagatorInField(transportMgr) ;

  //  fpSafetyHelper =   transportMgr->GetSafetyHelper();  // New 

  // Cannot determine whether a field exists here, as it would 
  //  depend on the relative order of creating the detector's 
  //  field and this process. That order is not guaranted.
  // Instead later the method DoesGlobalFieldExist() is called

  //  static G4TouchableHandle nullTouchableHandle;  // Points to (G4VTouchable*) 0
  //  fCurrentTouchableHandle = nullTouchableHandle; 

  This->fEndGlobalTimeComputed  = false;
  This->fCandidateEndGlobalTime = 0;

}

//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.

FQUALIFIER
G4double  GPTransportation_AlongStepGetPhysicalInteractionLength( 
                                GPTransportation *This,
			        GPTrack *track,
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

  //G4FWP On device, GPTrack is used instead of G4Track

  //  GPThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  //  GPThreeVector startPosition          = track.GetPosition() ;

  GPThreeVector startMomentumDir = GPThreeVector_unit(
				   GPThreeVector_create(track->px,
							track->py,
							track->pz));
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
  if( MagSqShift >= GPsqr(This->fPreviousSafety) )
  {
     currentSafety = 0.0 ;
  }
  else
  {
     currentSafety = This->fPreviousSafety - sqrt(MagSqShift) ;
  }

  // Is the particle charged or has it a magnetic moment?
  //
  //  G4double particleCharge = pParticle->GetCharge() ; 
  //  G4double magneticMoment = pParticle->GetMagneticMoment() ;
  //  G4double       restMass = pParticleDef->GetPDGMass() ;
  G4double particleCharge = track->q ; 
  G4double magneticMoment = 0.0 ;//@@@G4FWP temporary value
  G4double       restMass = 1.0 ;//@@@G4FWP temporary value

  This->fGeometryLimitedStep = false ;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"

  // Check if the particle has a force, EM or gravitational, exerted on it
  //
  GPFieldManager* fieldMgr = 0;
  G4bool          fieldExertsForce = false ;

  G4bool gravityOn = false;
  G4bool fieldExists= false;  // Field is not 0 (null pointer)

  //********************************************************************
  //@@@G4FWP use the magnetic field map on the device 
  //  fieldMgr = fFieldPropagator->FindAndSetFieldManager( track.GetVolume() );
  fieldMgr = GPPropagatorInField_FindAndSetFieldManager(This->fFieldPropagator);
  //@@@G4FWP************************************************************

  if( fieldMgr != 0 )
  {
    // Message the field Manager, to configure it for this track
    GPFieldManager_ConfigureForTrack(fieldMgr
		                     //	&track 
				     );
    // Is here to allow a transition from no-field pointer 
    //   to finite field (non-zero pointer).
    
    // If the field manager has no field ptr, the field is zero 
    //     by definition ( = there is no field ! )
    
    GPMagneticField* ptrField 
      = GPFieldManager_GetDetectorField(fieldMgr);
    fieldExists = (ptrField!=0) ;

    if( fieldExists ) 
    {
      gravityOn= GPMagneticField_IsGravityActive(ptrField);

      if(  (particleCharge != 0.0) 
	   || (This->fUseMagneticMoment && (magneticMoment != 0.0) )
	   || (gravityOn          && (restMass != 0.0) )
	   )
      {
	fieldExertsForce = fieldExists; 
      }
    }
  }
  
  if( !fieldExertsForce ) 
  {
     G4double linearStepLength ;
     if(This->fShortStepOptimisation && (currentMinimumStep <= currentSafety) )
     {
       // The Step is guaranteed to be taken
       //
       geometryStepLength   = currentMinimumStep ;
       This->fGeometryLimitedStep = false ;
     }
     else
     {
       //  Find whether the straight path intersects a volume
       //
       linearStepLength = GPNavigator_ComputeStep(This->fLinearNavigator, 
       						  startPosition, 
       						  startMomentumDir,
       						  currentMinimumStep, 
						  &newSafety); 
       // Remember last safety origin & value.
       //
       This->fPreviousSftOrigin = startPosition ;
       This->fPreviousSafety    = newSafety ; 
       //       fpSafetyHelper->SetCurrentSafety( newSafety, startPosition);

       currentSafety = newSafety ;
          
       This->fGeometryLimitedStep= (linearStepLength <= currentMinimumStep); 
       if( This->fGeometryLimitedStep )
       {
         // The geometry limits the Step size (an intersection was found.)
         geometryStepLength   = linearStepLength ;
       } 
       else
       {
         // The full Step is taken.
         geometryStepLength   = currentMinimumStep ;
       }
     }

     This->endpointDistance = geometryStepLength ;

     // Calculate final position
     //
     This->fTransportEndPosition = GPThreeVector_add(startPosition,
		   GPThreeVector_mult(startMomentumDir,geometryStepLength)) ;

     // Momentum direction, energy and polarisation are unchanged by transport
     //
     This->fTransportEndMomentumDir   = startMomentumDir ; 
     //     This->fTransportEndKineticEnergy = track.GetKineticEnergy() ;
     //     This->fTransportEndSpin          = track.GetPolarization();
     This->fTransportEndKineticEnergy = track->E ;
     This->fTransportEndSpin = GPThreeVector_create(0.,0.,0.);//G4FWP temporary
     This->fParticleIsLooping         = false ;
     This->fMomentumChanged           = false ; 
     This->fEndGlobalTimeComputed     = false ;
  }
  else   //  A field exerts force
  {
    //     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
    G4double       momentumMagnitude = track->E ; //G4FWP temporary
    //    GPThreeVector  EndUnitMomentum ;
    G4double       lengthAlongCurve ;
 
    GPPropagatorInField_SetChargeMomentumMass( This->fFieldPropagator,
    					       particleCharge,    // in e+ units
    					       momentumMagnitude, // in Mev/c 
    					       restMass           ) ;  

    //     GPThreeVector spin        = track.GetPolarization() ;
    GPThreeVector spin = GPThreeVector_create(0.,0.,0.) ;//G4FWP temporary

    GPFieldTrack  aFieldTrack;
    GPFieldTrack_Constructor2(&aFieldTrack,
			      startPosition, 
			      startMomentumDir, //track.GetMomentumDirection(),
			      0.0, 
			      track->E, //track.GetKineticEnergy(),
			      restMass,
			      c_light, //track.GetVelocity(),
			      0.0,     //track.GetGlobalTime(), // Lab.
			      0.0,     //track.GetProperTime(), // Part.
			      &spin                  ) ;
     if( currentMinimumStep > 0 ) 
     {
       // Do the Transport in the field (non recti-linear)
       //
       // @@@G4FWP 
       // get the physical volume of the current position and substitute 
       // track.GetVolume() in the original code

       GPVPhysicalVolume* cur_phy =
       	 GPNavigator_LocateGlobalPointAndSetup(This->fLinearNavigator,
					       startPosition,NULL,false,true);
       
       lengthAlongCurve = GPPropagatorInField_ComputeStep(
                                              This->fFieldPropagator, 
					      aFieldTrack,
					      currentMinimumStep, 
					      currentSafety, 
					      cur_phy //0 or track.GetVolume() //G4FWP 
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
       //	fpSafetyHelper->SetCurrentSafety( currentSafety, startPosition);
     }
     else
     {
       geometryStepLength   = lengthAlongCurve= 0.0 ;
       This->fGeometryLimitedStep = false ;
     }
      
     // Get the End-Position and End-Momentum (Dir-ection)
     //
     This->fTransportEndPosition = GPFieldTrack_GetPosition(&aFieldTrack) ;

     // Momentum:  Magnitude and direction can be changed too now ...
     //
     This->fMomentumChanged         = true ; 
     This->fTransportEndMomentumDir = GPFieldTrack_GetMomentumDir(&aFieldTrack);

     This->fTransportEndKineticEnergy = 
       GPFieldTrack_GetKineticEnergy(&aFieldTrack); 

     G4bool fieldChangeE = GPFieldManager_DoesFieldChangeEnergy(
       GPPropagatorInField_GetCurrentFieldManager(This->fFieldPropagator));
     if( fieldChangeE )
     {
       // If the field can change energy, then the time must be integrated
       //    - so this should have been updated
       //
       This->fCandidateEndGlobalTime  =
	 GPFieldTrack_GetLabTimeOfFlight(&aFieldTrack);
       This->fEndGlobalTimeComputed    = true;
       
       // was ( fCandidateEndGlobalTime != track.GetGlobalTime() );
       // a cleaner way is to have FieldTrack knowing whether time is updated.
     }
     else
     {
       // The energy should be unchanged by field transport,
       //    - so the time changed will be calculated elsewhere
       //
       This->fEndGlobalTimeComputed = false;

       // Check that the integration preserved the energy 
       //     -  and if not correct this!
       G4double  startEnergy= track->E; //G4FWP track.GetKineticEnergy();
       G4double  endEnergy= This->fTransportEndKineticEnergy; 

       //       static 
       G4int no_inexact_steps=0;
       G4int no_large_ediff = 0;
       G4double absEdiff = fabs(startEnergy- endEnergy);
       if( absEdiff > perMillion * endEnergy )
       {
	 no_inexact_steps++;
	 // Possible statistics keeping here ...
       }
       if( This->fVerboseLevel > 1 )
       {
	 if( fabs(startEnergy- endEnergy) > perThousand * endEnergy )
         {
	   //	   static 
	   G4int no_warnings= 0, warnModulo=1,  moduloFactor= 10; 
	   no_large_ediff++;
	   if( (no_large_ediff% warnModulo) == 0 )
           {
	     no_warnings++;
	     if( (This->fVerboseLevel > 2 ) || (no_warnings<4) || 
		 (no_large_ediff == warnModulo * moduloFactor) ){
	     }
	     if( no_large_ediff == warnModulo * moduloFactor )
             {
	       warnModulo *= moduloFactor;
	     }
	   }
	 }
       }  // end of if (fVerboseLevel)

       // Correct the energy for fields that conserve it
       //  This - hides the integration error
       //       - but gives a better physical answer
       This->fTransportEndKineticEnergy= track->E; // track.GetKineticEnergy(); 
     }

     This->fTransportEndSpin = GPFieldTrack_GetSpin(&aFieldTrack);
     This->fParticleIsLooping = 
       GPPropagatorInField_IsParticleLooping(This->fFieldPropagator) ;
     This->endpointDistance = GPThreeVector_mag(
	                GPThreeVector_sub(This->fTransportEndPosition,
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
      //G4FWP temporary
      // fpSafetyHelper->SetCurrentSafety(currentSafety,fTransportEndPosition);

      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength
      //
      currentSafety     += This->endpointDistance ;
      
    }
  }            

  //  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

