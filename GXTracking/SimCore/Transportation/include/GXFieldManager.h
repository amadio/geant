#ifndef GXFIELDMANAGER_HH
#define GXFIELDMANAGER_HH 1

//#include "globals.hh"
#include "GPTypeDef.h"

#include "GXMagneticField.h"
#include "GXChordFinder.h"
#include "GXFieldManager.h"
#include "GPChargeState.h"

// Forward reference for parameter configuration
//class G4Field;
//class G4MagneticField;
//class G4ChordFinder;
//class G4Track;

struct GXFieldManager
{
     // Dependent objects -- with state that depends on tracking
     GXMagneticField*        fDetectorField;
     GXChordFinder*  fChordFinder;

     G4bool          fAllocatedChordFinder; // Did we used "new" to
					    // create fChordFinder ?
     // INVARIANTS of tracking  ---------------------------------------
     // 
     //  1. CONSTANTS 
     G4double  fEpsilonMinDefault;   // Can be 1.0e-5 to 1.0e-10 ...
     G4double  fEpsilonMaxDefault;   // Can be 1.0e-3 to 1.0e-8 ...

     //  2. CHARACTERISTIC of field
     G4bool          fFieldChangesEnergy;

     //  3. PARAMETERS 
     // 
     //     Values for the required accuracies
     G4double  fDelta_One_Step_Value;      //  for one tracking/physics step
     G4double  fDelta_Intersection_Val;    //  for boundary intersection

     G4double  fDefault_Delta_One_Step_Value;   // = 0.25 * mm;
     G4double  fDefault_Delta_Intersection_Val; // = 0.1 * mm;

     //     Values for the small possible relative accuracy of a step
     //     (corresponding to the greatest possible integration accuracy)
     G4double  fEpsilonMin; 
     G4double  fEpsilonMax;

};

// Our current design and implementation expect that a particular
// geometrical region has a Field manager.
//  By default a Field Manager is created for the world volume, and
//  will be utilised for all volumes unless it is overridden by a 'local'
//  field manager.

// Note also that a region with both electric E and magnetic B field will 
//  have these treated as one field.
// Similarly it could be extended to treat other fields as additional components
//  of a single field type.

// Implementation of inline functions

extern "C" {

FQUALIFIER
void GXFieldManager_Constructor( GXFieldManager *This,
				 GXMagneticField       *detectorField, 
				 GXChordFinder *pChordFinder, 
				 G4bool        fieldChangesEnergy);

FQUALIFIER
void GXFieldManager_Constructor2(GXFieldManager *This,
				 GXMagneticField *detectorField);

FQUALIFIER
void GXFieldManager_ConfigureForTrack( GXFieldManager *This);

FQUALIFIER
void GXFieldManager_CreateChordFinder( GXFieldManager *This,
				       GXMagneticField *detectorMagField);

FQUALIFIER
G4bool GXFieldManager_SetDetectorField(GXFieldManager *This,
				       GXMagneticField *pDetectorField);

FQUALIFIER
GXMagneticField* GXFieldManager_GetDetectorField(GXFieldManager *This);

FQUALIFIER
G4bool GXFieldManager_DoesFieldExist(GXFieldManager *This);

FQUALIFIER  
void GXFieldManager_SetChordFinder(GXFieldManager *This, 
				   GXChordFinder *aChordFinder);

FQUALIFIER  
GXChordFinder*  GXFieldManager_GetChordFinder(GXFieldManager *This);

FQUALIFIER
G4double GXFieldManager_GetDeltaIntersection(GXFieldManager *This);

FQUALIFIER
G4double GXFieldManager_GetDeltaOneStep(GXFieldManager *This);

FQUALIFIER
void GXFieldManager_SetDeltaOneStep(GXFieldManager *This,
				    G4double valDeltaOneStep);

FQUALIFIER
void GXFieldManager_SetDeltaIntersection(GXFieldManager *This,
					 G4double valDeltaIntersection);

FQUALIFIER
void GXFieldManager_SetAccuraciesWithDeltaOneStep(GXFieldManager *This,
						  G4double valDeltaOneStep);

FQUALIFIER 
G4bool GXFieldManager_DoesFieldChangeEnergy(GXFieldManager *This);

FQUALIFIER 
void GXFieldManager_SetFieldChangesEnergy(GXFieldManager *This,
					  G4bool value);

FQUALIFIER 
G4double  GXFieldManager_GetMinimumEpsilonStep(GXFieldManager *This);

FQUALIFIER 
void GXFieldManager_SetMinimumEpsilonStep(GXFieldManager *This,
					  G4double newEpsMin );
FQUALIFIER 
G4double  GXFieldManager_GetMaximumEpsilonStep(GXFieldManager *This);

FQUALIFIER 
void GXFieldManager_SetMaximumEpsilonStep(GXFieldManager *This,
					  G4double newEpsMax );

}

#endif 
