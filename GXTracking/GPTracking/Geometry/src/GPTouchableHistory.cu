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
// $Id: G4TouchableHistory.cc,v 1.15 2009-11-06 11:10:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4TouchableHistory Implementation
//
// ----------------------------------------------------------------------

#include "GPVPhysicalVolume.h"
#include "GPTouchableHistory.h"

//G4Allocator<GPTouchableHistory> aTouchableHistoryAllocator;

FQUALIFIER
void GPTouchableHistory_Constructor(GPTouchableHistory *This)
{
  This->frot = GPRotationMatrix_create(1,0,0,
				       0,1,0,
				       0,0,1);
  This->ftlate = GPThreeVector_create(0,0,0);

  GPNavigationHistory_Constructor(&(This->fhistory));

  GPVPhysicalVolume* pPhysVol = 0;
  GPNavigationHistory_SetFirstEntry(&(This->fhistory),pPhysVol);
}

/*
FQUALIFIER
void GPTouchableHistory_Constructor2(GPTouchableHistory *This,
				     GPNavigator *navigator )
//  : fhistory(history)
{
  This->fhistory = history;

  //  G4AffineTransform tf(fhistory.GetTopTransform().Inverse());
  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(&(This->fhistory));
  GPAffineTransform tf = GPAffineTransform_Inverse(&aT); 
  This->ftlate = GPAffineTransform_NetTranslation(&tf);
  This->frot = GPAffineTransform_NetRotation(&tf);
}
*/

/*
GPTouchableHistory_~GPTouchableHistory()
{
}
*/

FQUALIFIER
GPThreeVector
GPTouchableHistory_GetTranslation(GPTouchableHistory *This,
				  G4int depth)
{
  // The value returned will change at the next call
  // Copy it if you want to use it!
  //
  //  static GPThreeVector currTranslation;

  if(depth==0.0)
  {
    return This->ftlate;
  }
  else
  {
    GPAffineTransform aT= GPNavigationHistory_GetTransform(&(This->fhistory), 
			  GPTouchableHistory_CalculateHistoryIndex(This,depth));
    GPThreeVector currTranslation = GPAffineTransform_NetTranslation(&aT);
    return currTranslation;
  }
}

FQUALIFIER
GPRotationMatrix
GPTouchableHistory_GetRotation(GPTouchableHistory *This,
			       G4int depth)
{
  // The value returned will change at the next call
  // Copy it if you want to use it!
  //
  //  static GPRotationMatrix rotM;

  if(depth==0.0)
  {
    return This->frot;
  }
  else
  {
    GPAffineTransform aT= GPNavigationHistory_GetTransform(&(This->fhistory),
			  GPTouchableHistory_CalculateHistoryIndex(This,depth));
    GPRotationMatrix rotM = GPAffineTransform_NetRotation(&aT);
    return rotM;
  }
}
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
// $Id: G4TouchableHistory.icc,v 1.15 2010-10-27 07:34:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4TouchableHistory inline implementation
// ----------------------------------------------------------------------

//#if defined G4GEOM_ALLOC_EXPORT
//extern G4DLLEXPORT G4Allocator<G4TouchableHistory> aTouchableHistoryAllocator;
//#else
//extern G4DLLIMPORT G4Allocator<G4TouchableHistory> aTouchableHistoryAllocator;
//#endif

FQUALIFIER
void  GPTouchableHistory_UpdateYourself( GPTouchableHistory *This,
					 GPVPhysicalVolume*  pPhysVol,
					 GPNavigationHistory* pHistory ) 
{
  This->fhistory = *pHistory;  
  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(&(This->fhistory));
  GPAffineTransform tf = GPAffineTransform_Inverse(&aT); 

  if( pPhysVol == 0 )
  {
    // This means that the track has left the World Volume.
    // Since the Navigation History does not already reflect this,
    // we must correct this problem here.
    //
    GPNavigationHistory_SetFirstEntry(&(This->fhistory),pPhysVol);
  }
  This->ftlate = GPAffineTransform_NetTranslation(&tf);
  This->frot = GPAffineTransform_NetRotation(&tf);
}

FQUALIFIER
G4int GPTouchableHistory_CalculateHistoryIndex( GPTouchableHistory *This,
						G4int stackDepth )
{
  return (GPNavigationHistory_GetDepth(&(This->fhistory))-stackDepth); // was -1
}

FQUALIFIER
GPVPhysicalVolume* GPTouchableHistory_GetVolume( GPTouchableHistory *This,
						 G4int depth )
{
  return GPNavigationHistory_GetVolume(&(This->fhistory),
	 GPTouchableHistory_CalculateHistoryIndex(This,depth));
}

FQUALIFIER
GPVSolid* GPTouchableHistory_GetSolid( GPTouchableHistory *This,
				       G4int depth )
{
  return GPLogicalVolume_GetSolid(
           GPVPhysicalVolume_GetLogicalVolume(
             GPNavigationHistory_GetVolume(&(This->fhistory),
	     GPTouchableHistory_CalculateHistoryIndex(This,depth))
					      ));
}

FQUALIFIER
G4int GPTouchableHistory_GetReplicaNumber( GPTouchableHistory *This,
					   G4int depth )
{
  return GPNavigationHistory_GetReplicaNo(&(This->fhistory),
	 GPTouchableHistory_CalculateHistoryIndex(This,depth));
}

FQUALIFIER
G4int GPTouchableHistory_GetHistoryDepth( GPTouchableHistory *This)
{
  return  GPNavigationHistory_GetDepth(&(This->fhistory));
}

FQUALIFIER
G4int GPTouchableHistory_MoveUpHistory( GPTouchableHistory *This,
					G4int num_levels )
{
  G4int maxLevelsMove = GPNavigationHistory_GetDepth(&(This->fhistory));
  G4int minLevelsMove = 0;              // Cannot redescend today!
                                        // Soon it will be possible
                                        // by adding a data member here
                                        //     fCurrentDepth;
  if( num_levels > maxLevelsMove )
  {
    num_levels = maxLevelsMove;
  }
  else if( num_levels < minLevelsMove )
  {
    num_levels = minLevelsMove;
  }
  GPNavigationHistory_BackLevel2(&(This->fhistory), num_levels ); 

  return num_levels;
}

FQUALIFIER
GPNavigationHistory* GPTouchableHistory_GetHistory( GPTouchableHistory *This )
{
  return &(This->fhistory);
}

// There is no provision in case this class is subclassed.
// If it is subclassed, this will fail and may not give errors!
//
/*
FQUALIFIER
void* GPTouchableHistory_operator new(size_t)
{
  return (void *) aTouchableHistoryAllocator.MallocSingle();
}

FQUALIFIER
void GPTouchableHistory_operator delete(void *aTH)
{
  aTouchableHistoryAllocator.FreeSingle((G4TouchableHistory *) aTH);
}
*/

//---------------------------------------------------------------------
//from G4VTouchable.icc
//---------------------------------------------------------------------

//G4int G4VTouchable::GetCopyNumber(G4int depth) const
FQUALIFIER
G4int GPTouchableHistory_GetCopyNumber(GPTouchableHistory *This, 
				       G4int depth)
{
  return GPTouchableHistory_GetReplicaNumber(This,depth);
}
