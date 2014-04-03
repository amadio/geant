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
// $Id: G4FieldManager.hh,v 1.16 2006-06-29 18:22:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
// class G4FieldManager
//
// Class description:
//
// A class to manage (Store) a pointer to the Field subclass that
// describes the field of a detector (magnetic, electric or other).
// Also stores a reference to the chord finder.
//
// The G4FieldManager class exists to allow the user program to specify 
// the electric, magnetic and/or other field(s) of the detector.
// 
// A field manager can be set to a logical volume (or to more than one), 
// in order to vary its field from that of the world.  In this manner
// a zero or constant field can override a global field,  a more or 
// less exact version can override the external approximation, lower
// or higher precision for tracking can be specified, a different 
// stepper can be chosen for different volumes, ...
//
// It also stores a pointer to the ChordFinder object that can do the
// propagation in this field. All geometrical track "advancement" 
// in the field is handled by this ChordFinder object.
//
// G4FieldManager allows the other classes/object (of the MagneticField 
// & other class categories) to find out whether a detector field object 
// exists and what that object is.
//
// The Chord Finder must be created either by calling CreateChordFinder
// for a Magnetic Field or by the user creating a  a Chord Finder object
// "manually" and setting this pointer.
//
// A default FieldManager is created by the singleton class
// G4NavigatorForTracking and exists before main is called.
// However a new one can be created and given to G4NavigatorForTracking.
//
// Our current design envisions that one Field manager is 
// valid for each region detector.

// History:
// - 05.11.03 John Apostolakis, Added Min/MaximumEpsilonStep
// - 20.06.03 John Apostolakis, Abstract & ability to ConfigureForTrack
// - 10.03.97 John Apostolakis, design and implementation.
// -------------------------------------------------------------------

#ifndef GPFIELDMANAGER_HH
#define GPFIELDMANAGER_HH 1

//#include "globals.hh"
#include "GPTypeDef.h"
#include "GPMagneticField.h"
#include "GPChordFinder.h"
#include "GPFieldManager.h"
#include "GPChargeState.h"

// Forward reference for parameter configuration
//class G4Field;
//class G4MagneticField;
//class G4ChordFinder;
//class G4Track;

struct GPFieldManager
{
     // Dependent objects -- with state that depends on tracking
     GPMagneticField*        fDetectorField;
     GPChordFinder*  fChordFinder;

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
void GPFieldManager_Constructor( GPFieldManager *This,
				 GPMagneticField       *detectorField, 
				 GPChordFinder *pChordFinder, 
				 G4bool        fieldChangesEnergy);

FQUALIFIER
void GPFieldManager_Constructor2(GPFieldManager *This,
				 GPMagneticField *detectorField);

FQUALIFIER
void GPFieldManager_ConfigureForTrack( GPFieldManager *This);

FQUALIFIER
void GPFieldManager_CreateChordFinder( GPFieldManager *This,
				       GPMagneticField *detectorMagField);

FQUALIFIER
G4bool GPFieldManager_SetDetectorField(GPFieldManager *This,
				       GPMagneticField *pDetectorField);

FQUALIFIER
GPMagneticField* GPFieldManager_GetDetectorField(GPFieldManager *This);

FQUALIFIER
G4bool GPFieldManager_DoesFieldExist(GPFieldManager *This);

FQUALIFIER  
void GPFieldManager_SetChordFinder(GPFieldManager *This, 
				   GPChordFinder *aChordFinder);

FQUALIFIER  
GPChordFinder*  GPFieldManager_GetChordFinder(GPFieldManager *This);

FQUALIFIER
G4double GPFieldManager_GetDeltaIntersection(GPFieldManager *This);

FQUALIFIER
G4double GPFieldManager_GetDeltaOneStep(GPFieldManager *This);

FQUALIFIER
void GPFieldManager_SetDeltaOneStep(GPFieldManager *This,
				    G4double valDeltaOneStep);

FQUALIFIER
void GPFieldManager_SetDeltaIntersection(GPFieldManager *This,
					 G4double valDeltaIntersection);

FQUALIFIER
void GPFieldManager_SetAccuraciesWithDeltaOneStep(GPFieldManager *This,
						  G4double valDeltaOneStep);

FQUALIFIER 
G4bool GPFieldManager_DoesFieldChangeEnergy(GPFieldManager *This);

FQUALIFIER 
void GPFieldManager_SetFieldChangesEnergy(GPFieldManager *This,
					  G4bool value);

FQUALIFIER 
G4double  GPFieldManager_GetMinimumEpsilonStep(GPFieldManager *This);

FQUALIFIER 
void GPFieldManager_SetMinimumEpsilonStep(GPFieldManager *This,
					  G4double newEpsMin );
FQUALIFIER 
G4double  GPFieldManager_GetMaximumEpsilonStep(GPFieldManager *This);

FQUALIFIER 
void GPFieldManager_SetMaximumEpsilonStep(GPFieldManager *This,
					  G4double newEpsMax );

}

#endif 
