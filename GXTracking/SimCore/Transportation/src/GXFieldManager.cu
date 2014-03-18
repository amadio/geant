#include "GXFieldManager.h"
#include "GXMagneticField.h"
#include "GXChordFinder.h"

FQUALIFIER
void GXFieldManager_Constructor( GXFieldManager *This,
				 GXMagneticField       *detectorField, 
				 GXChordFinder *pChordFinder, 
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

    This->fFieldChangesEnergy=  false;
  //      GXMagneticField_DoesFieldChangeEnergy(detectorField);
  else
    This->fFieldChangesEnergy= fieldChangesEnergy;

  // Add to store
  //  GXFieldManagerStore_Register(This);
}

FQUALIFIER
void GXFieldManager_Constructor2(GXFieldManager *This,
				 GXMagneticField *detectorField)
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
  GXChordFinder fChordFinder;
  GXChordFinder_Constructor2(&fChordFinder, detectorField,0,0 );

  This->fDelta_One_Step_Value= This->fDefault_Delta_One_Step_Value;
  This->fDelta_Intersection_Val= This->fDefault_Delta_Intersection_Val;
  // Add to store
  //  GXFieldManagerStore_Register(This);
}

FQUALIFIER
void GXFieldManager_ConfigureForTrack( GXFieldManager *This
				       //const G4Track * 
				       ) 
{
   // Default is to do nothing!
   ;
}

FQUALIFIER
void GXFieldManager_CreateChordFinder( GXFieldManager *This,
				       GXMagneticField *detectorMagField)
{
  //   if ( This->fAllocatedChordFinder ) delete fChordFinder;
  //   fChordFinder= new G4ChordFinder( detectorMagField );
  GXChordFinder newChordFinder;
  GXChordFinder_Constructor2(&newChordFinder, detectorMagField, 0,0 );
  This->fChordFinder = &newChordFinder;
  This->fAllocatedChordFinder= true;
}

FQUALIFIER
G4bool GXFieldManager_SetDetectorField(GXFieldManager *This,
				       GXMagneticField *pDetectorField)
{
   This->fDetectorField= pDetectorField;

   if ( pDetectorField )
     This->fFieldChangesEnergy = false;
   //       GXMagneticField_DoesFieldChangeEnergy(pDetectorField);
   else
     This->fFieldChangesEnergy= false;   //  No field 

   return false;
}

//
// G4FieldManager inline implementation 
//
// -------------------------------------------------------------------


FQUALIFIER
GXMagneticField* GXFieldManager_GetDetectorField(GXFieldManager *This)
{ 
  // If pointer is null, should this raise an exception ??
  return This->fDetectorField;
} 

FQUALIFIER
G4bool GXFieldManager_DoesFieldExist(GXFieldManager *This)
{ 
  return (This->fDetectorField != 0);
} 

FQUALIFIER  
void GXFieldManager_SetChordFinder(GXFieldManager *This, 
				   GXChordFinder *aChordFinder)
{
  This->fChordFinder= aChordFinder;
}

FQUALIFIER  
GXChordFinder*  GXFieldManager_GetChordFinder(GXFieldManager *This)
{  
  return This->fChordFinder;
}

FQUALIFIER
G4double GXFieldManager_GetDeltaIntersection(GXFieldManager *This)
{
  return This->fDelta_Intersection_Val;
}

FQUALIFIER
G4double GXFieldManager_GetDeltaOneStep(GXFieldManager *This)
{
  return This->fDelta_One_Step_Value;
}

FQUALIFIER
void GXFieldManager_SetDeltaOneStep(GXFieldManager *This,
				    G4double valDeltaOneStep)
{ 
  This->fDelta_One_Step_Value= valDeltaOneStep;  
}

FQUALIFIER
void GXFieldManager_SetDeltaIntersection(GXFieldManager *This,
					 G4double valDeltaIntersection)
{
  This->fDelta_Intersection_Val = valDeltaIntersection;
}

FQUALIFIER
void GXFieldManager_SetAccuraciesWithDeltaOneStep(GXFieldManager *This,
						  G4double valDeltaOneStep)
{ 
  This->fDelta_One_Step_Value= valDeltaOneStep;  
  This->fDelta_Intersection_Val = 0.4 * This->fDelta_One_Step_Value;
}

FQUALIFIER 
G4bool GXFieldManager_DoesFieldChangeEnergy(GXFieldManager *This)
{ 
  return This->fFieldChangesEnergy;
}

FQUALIFIER 
void GXFieldManager_SetFieldChangesEnergy(GXFieldManager *This,
					  G4bool value)
{ 
  This->fFieldChangesEnergy = value; 
}

// Minimum for Relative accuracy of any Step 
FQUALIFIER 
G4double  GXFieldManager_GetMinimumEpsilonStep(GXFieldManager *This)
{
  return This->fEpsilonMin; 
}

FQUALIFIER 
void GXFieldManager_SetMinimumEpsilonStep(GXFieldManager *This,
					  G4double newEpsMin )
{
  if( (newEpsMin > 0.0) && (fabs(1.0+newEpsMin) > 1.0) )
  {
    This->fEpsilonMin = newEpsMin;
  }
}

// Maximum for Relative accuracy of any Step 
FQUALIFIER 
G4double  GXFieldManager_GetMaximumEpsilonStep(GXFieldManager *This)
{
  return This->fEpsilonMax; 
}

FQUALIFIER 
void GXFieldManager_SetMaximumEpsilonStep(GXFieldManager *This,
					  G4double newEpsMax )
{
  if(    (newEpsMax > 0.0) 
      && (newEpsMax >= This->fEpsilonMin ) 
      && (fabs(1.0+newEpsMax)>1.0) )
  {
    This->fEpsilonMax = newEpsMax;
  }
}
