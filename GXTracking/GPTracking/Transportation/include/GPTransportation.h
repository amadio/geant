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
// $Id: G4Transportation.hh,v 1.17 2007-11-09 15:39:20 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//        GEANT 4  include file implementation
// ------------------------------------------------------------
//
// Class description:
//
// G4Transportation is a process responsible for the transportation of 
// a particle, i.e. the geometrical propagation encountering the 
// geometrical sub-volumes of the detectors.
// It is also tasked with part of updating the "safety".

// =======================================================================
// Created:  19 March 1997, J. Apostolakis
// =======================================================================
#ifndef GPTransportation_hh
#define GPTransportation_hh 1

//#include "G4VProcess.hh"
#include "GPFieldManager.h"

#include "GPNavigator.h"
#include "GPTransportationManager.h"
#include "GPTransportation.h"
#include "GPPropagatorInField.h"
#include "GPTrack.h"
//#include "G4Step.hh"
//#include "G4ParticleChangeForTransport.hh"
//class G4SafetyHelper; 

#include "GPTypeDef.h"
#include "GPThreeVector.h"

#include "GPGPILSelection.h"

struct GPTransportation
{
  GPNavigator*         fLinearNavigator;
  GPPropagatorInField* fFieldPropagator;
  // The Propagators used to transport the particle
  
  GPThreeVector        fTransportEndPosition;
  GPThreeVector        fTransportEndMomentumDir;
  G4double             fTransportEndKineticEnergy;
  GPThreeVector        fTransportEndSpin;
  G4bool               fMomentumChanged;
  G4bool               fEnergyChanged;
  G4bool               fEndGlobalTimeComputed; 
  G4double             fCandidateEndGlobalTime;
  // The particle's state after this Step, Store for DoIt
  
  G4bool               fParticleIsLooping;
  
  //  G4TouchableHandle    fCurrentTouchableHandle;
  
  G4bool fGeometryLimitedStep;
  // Flag to determine whether a boundary was reached.
  
  GPThreeVector  fPreviousSftOrigin;
  G4double       fPreviousSafety; 
  // Remember last safety origin & value.
  
  //  G4ParticleChangeForTransport fParticleChange;
  // New ParticleChange
  
  G4double endpointDistance;
  
  // Thresholds for looping particles: 
  // 
  G4double fThreshold_Warning_Energy;     //  Warn above this energy
  G4double fThreshold_Important_Energy;   //  Hesitate above this
  G4int    fThresholdTrials;              //    for this no of trials
  // Above 'important' energy a 'looping' particle in field will 
  //   *NOT* be abandoned, except after fThresholdTrials attempts.
  G4double fUnimportant_Energy;
  //  Below this energy, no verbosity for looping particles is issued
  
  // Counter for steps in which particle reports 'looping',
  //   if it is above 'Important' Energy 
  G4int    fNoLooperTrials; 
  // Statistics for tracks abandoned
  G4double fSumEnergyKilled;
  G4double fMaxEnergyKilled;
  
  // Whether to avoid calling G4Navigator for short step ( < safety)
  //   If using it, the safety estimate for endpoint will likely be smaller.
  G4bool   fShortStepOptimisation; 
  
  // Whether to track state change from magnetic moment in a B-field
  G4bool   fUseMagneticMoment; 
  
  //  G4SafetyHelper* fpSafetyHelper;  // To pass it the safety value obtained
  
  // Verbosity 
  G4int    fVerboseLevel;
  // Verbosity level for warnings
  // eg about energy non-conservation in magnetic field.
};

extern "C" {

FQUALIFIER
void GPTransportation_Constructor( GPTransportation *This, 
                                   GPPropagatorInField *propagator,
                                   G4int verboseLevel );

FQUALIFIER
void GPTransportation_Constructor2( GPTransportation *This, 
                                    GPTransportationManager *transportMgr,
                                    G4int verboseLevel );

FQUALIFIER
G4double  GPTransportation_AlongStepGetPhysicalInteractionLength( 
                                GPTransportation *This,
                                GPTrack *track,
                                G4double  previousStepSize,
                                G4double  currentMinimumStep,
                                G4double  currentSafety,
                                GPGPILSelection* selection );

}

//#include "GPTransportation.c"

#endif  
