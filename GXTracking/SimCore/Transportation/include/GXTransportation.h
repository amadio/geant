#ifndef GXTransportation_hh
#define GXTransportation_hh 1

#include "GXTransportation.h"
#include "GXPropagatorInField.h"
#include "GXFieldManager.h"

#include "GPNavigator.h"

#include "GXTrack.h"
#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GPGPILSelection.h"

struct GXTransportation
{
  GPNavigator*         fLinearNavigator;
  GXPropagatorInField* fFieldPropagator;
  // The Propagators used to transport the particle
  
  GPThreeVector        fTransportEndPosition;
  GPThreeVector        fTransportEndMomentumDir;
  G4double             fTransportEndKineticEnergy;
  G4bool               fMomentumChanged;
  //  G4bool               fEnergyChanged;
  //  G4bool               fEndGlobalTimeComputed; 
  //  G4double             fCandidateEndGlobalTime;
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
  
};

extern "C" {

FQUALIFIER
void GXTransportation_Constructor( GXTransportation *This, 
                                   GXPropagatorInField *propagator,
                                   G4int verboseLevel );

FQUALIFIER
G4double  GXTransportation_AlongStepGetPhysicalInteractionLength( 
                                GXTransportation *This,
                                GXTrack *track,
                                G4double  previousStepSize,
                                G4double  currentMinimumStep,
                                G4double  currentSafety,
                                GPGPILSelection* selection );

FQUALIFIER
G4double  GXTransportation_AlongStepGPIL_Photon( 
                                GXTransportation *This,
                                GXTrack *track,
                                G4double  previousStepSize,
                                G4double  currentMinimumStep,
                                G4double  currentSafety,
                                GPGPILSelection* selection );

FQUALIFIER
G4double  GXTransportation_AlongStepGPIL_Electron( 
                                GXTransportation *This,
                                GXTrack *track,
                                G4double  previousStepSize,
                                G4double  currentMinimumStep,
                                G4double  currentSafety,
                                GPGPILSelection* selection );
}

#endif  
