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
// $Id: G4FieldManager.cc,v 1.15 2007-12-07 15:34:10 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//#include "G4Field.hh"
//#include "GPFieldManagerStore.h"

#include "GPFieldManager.h"
#include "GPMagneticField.h"
#include "GPChordFinder.h"

FQUALIFIER
void GPFieldManager_Constructor( GPFieldManager *This,
				 GPMagneticField       *detectorField, 
				 GPChordFinder *pChordFinder, 
				 G4bool        fieldChangesEnergy
				 )
{ 
  This->fDetectorField = detectorField; 
  This->fChordFinder = pChordFinder; 
  This->fAllocatedChordFinder = false;
  This->fEpsilonMinDefault = 5.0e-5; 
  This->fEpsilonMaxDefault = 0.001;
  This->fDefault_Delta_One_Step_Value = 0.01*millimeter; 
  This->fDefault_Delta_Intersection_Val = 0.001*millimeter;
  This->fEpsilonMin = This->fEpsilonMinDefault ;
  This->fEpsilonMax = This->fEpsilonMaxDefault;

  This->fDelta_One_Step_Value = This->fDefault_Delta_One_Step_Value;
  This->fDelta_Intersection_Val= This->fDefault_Delta_Intersection_Val;
  if ( detectorField )
    This->fFieldChangesEnergy= 
      GPMagneticField_DoesFieldChangeEnergy(detectorField);
  else
    This->fFieldChangesEnergy= fieldChangesEnergy;

  // Add to store
  //  GPFieldManagerStore_Register(This);
}

FQUALIFIER
void GPFieldManager_Constructor2(GPFieldManager *This,
				 GPMagneticField *detectorField)
{
  This->fDetectorField = detectorField;

  This->fFieldChangesEnergy = false; 
  This->fAllocatedChordFinder = false;
  This->fEpsilonMinDefault = 5.0e-5; 
  This->fEpsilonMaxDefault = 0.001;
  This->fDefault_Delta_One_Step_Value = 0.01*millimeter; 
  This->fDefault_Delta_Intersection_Val = 0.001*millimeter;
  This->fEpsilonMin = This->fEpsilonMinDefault ;
  This->fEpsilonMax = This->fEpsilonMaxDefault;


  //   fChordFinder= new G4ChordFinder( detectorField );
  GPChordFinder fChordFinder;
  GPChordFinder_Constructor2(&fChordFinder, detectorField,0,0 );

  This->fDelta_One_Step_Value= This->fDefault_Delta_One_Step_Value;
  This->fDelta_Intersection_Val= This->fDefault_Delta_Intersection_Val;
  // Add to store
  //  GPFieldManagerStore_Register(This);
}

FQUALIFIER
void GPFieldManager_ConfigureForTrack( GPFieldManager *This
				       //const G4Track * 
				       ) 
{
   // Default is to do nothing!
   ;
}

//GPFieldManager::~GPFieldManager()
//{
//   if( fAllocatedChordFinder ){
//      delete fChordFinder;
//   }
//   GPFieldManagerStore::DeRegister(this);
//}

FQUALIFIER
void GPFieldManager_CreateChordFinder( GPFieldManager *This,
				       GPMagneticField *detectorMagField)
{
  //   if ( This->fAllocatedChordFinder ) delete fChordFinder;
  //   fChordFinder= new G4ChordFinder( detectorMagField );
  GPChordFinder newChordFinder;
  GPChordFinder_Constructor2(&newChordFinder, detectorMagField, 0,0 );
  This->fChordFinder = &newChordFinder;
  This->fAllocatedChordFinder= true;
}

FQUALIFIER
G4bool GPFieldManager_SetDetectorField(GPFieldManager *This,
				       GPMagneticField *pDetectorField)
{
   This->fDetectorField= pDetectorField;

   if ( pDetectorField )
     This->fFieldChangesEnergy = 
       GPMagneticField_DoesFieldChangeEnergy(pDetectorField);
   else
     This->fFieldChangesEnergy= false;   //  No field 

   return false;
}

//
// G4FieldManager inline implementation 
//
// -------------------------------------------------------------------


FQUALIFIER
GPMagneticField* GPFieldManager_GetDetectorField(GPFieldManager *This)
{ 
  // If pointer is null, should this raise an exception ??
  return This->fDetectorField;
} 

FQUALIFIER
G4bool GPFieldManager_DoesFieldExist(GPFieldManager *This)
{ 
  return (This->fDetectorField != 0);
} 

FQUALIFIER  
void GPFieldManager_SetChordFinder(GPFieldManager *This, 
				   GPChordFinder *aChordFinder)
{
  This->fChordFinder= aChordFinder;
}

FQUALIFIER  
GPChordFinder*  GPFieldManager_GetChordFinder(GPFieldManager *This)
{  
  return This->fChordFinder;
}

FQUALIFIER
G4double GPFieldManager_GetDeltaIntersection(GPFieldManager *This)
{
  return This->fDelta_Intersection_Val;
}

FQUALIFIER
G4double GPFieldManager_GetDeltaOneStep(GPFieldManager *This)
{
  return This->fDelta_One_Step_Value;
}

FQUALIFIER
void GPFieldManager_SetDeltaOneStep(GPFieldManager *This,
				    G4double valDeltaOneStep)
{ 
  This->fDelta_One_Step_Value= valDeltaOneStep;  
}

FQUALIFIER
void GPFieldManager_SetDeltaIntersection(GPFieldManager *This,
					 G4double valDeltaIntersection)
{
  This->fDelta_Intersection_Val = valDeltaIntersection;
}

FQUALIFIER
void GPFieldManager_SetAccuraciesWithDeltaOneStep(GPFieldManager *This,
						  G4double valDeltaOneStep)
{ 
  This->fDelta_One_Step_Value= valDeltaOneStep;  
  This->fDelta_Intersection_Val = 0.4 * This->fDelta_One_Step_Value;
}

FQUALIFIER 
G4bool GPFieldManager_DoesFieldChangeEnergy(GPFieldManager *This)
{ 
  return This->fFieldChangesEnergy;
}

FQUALIFIER 
void GPFieldManager_SetFieldChangesEnergy(GPFieldManager *This,
					  G4bool value)
{ 
  This->fFieldChangesEnergy = value; 
}

// Minimum for Relative accuracy of any Step 
FQUALIFIER 
G4double  GPFieldManager_GetMinimumEpsilonStep(GPFieldManager *This)
{
  return This->fEpsilonMin; 
}

FQUALIFIER 
void GPFieldManager_SetMinimumEpsilonStep(GPFieldManager *This,
					  G4double newEpsMin )
{
  if( (newEpsMin > 0.0) && (fabs(1.0+newEpsMin) > 1.0) )
  {
    This->fEpsilonMin = newEpsMin;
  }
}

// Maximum for Relative accuracy of any Step 
FQUALIFIER 
G4double  GPFieldManager_GetMaximumEpsilonStep(GPFieldManager *This)
{
  return This->fEpsilonMax; 
}

FQUALIFIER 
void GPFieldManager_SetMaximumEpsilonStep(GPFieldManager *This,
					  G4double newEpsMax )
{
  if(    (newEpsMax > 0.0) 
      && (newEpsMax >= This->fEpsilonMin ) 
      && (fabs(1.0+newEpsMax)>1.0) )
  {
    This->fEpsilonMax = newEpsMax;
  }
}
