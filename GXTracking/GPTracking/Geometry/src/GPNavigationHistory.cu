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
// $Id: G4NavigationHistory.cc,v 1.17 2010-12-15 17:05:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4NavigationHistory Implementation
//
// Author: P.Kent, August 96
//
// ----------------------------------------------------------------------

#include "GPNavigationHistory.h"
//#include "G4ios.hh"

//#ifndef WIN32
// Initialise static data for the specialized memory pool
// for the internal STL vector of histories  ...
//
//  G4ChunkIndexType* G4AllocStats::allocStat = 0;
//  G4int             G4AllocStats::totSpace = 0;
//  G4int             G4AllocStats::numCat = 0;
//#endif

FQUALIFIER
void GPNavigationHistory_Constructor(GPNavigationHistory *This)
{
  //  fNavHistory(kHistoryMax)
  This->fStackDepth = 0 ; 
  GPNavigationHistory_Clear(This);
}

//G4NavigationHistory::G4NavigationHistory(const G4NavigationHistory &h)
//  : fNavHistory(h.fNavHistory), fStackDepth(h.fStackDepth)
//{
//}

//G4NavigationHistory::~G4NavigationHistory()
//{
//}

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
// $Id: G4NavigationHistory.icc,v 1.14 2009-08-04 08:27:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4NavigationHistory Inline implementation
//
// ----------------------------------------------------------------------

FQUALIFIER
void GPNavigationHistory_Reset(GPNavigationHistory *This)
{
  This->fStackDepth=0;
}

FQUALIFIER
void GPNavigationHistory_Clear(GPNavigationHistory *This)
{
  GPAffineTransform origin;
  GPAffineTransform_Vector(&origin, GPThreeVector_create(0.,0.,0.));

  GPNavigationLevel tmpNavLevel;
  GPNavigationLevel_Constructor(&tmpNavLevel, NULL, origin, kNormal, -1) ;

  GPNavigationHistory_Reset(This);

  //  for (register G4int ilev=fNavHistory.size()-1; ilev>=0; ilev--)
  for (G4int ilev=kHistoryMax-1; ilev>=0; ilev--)
  {
     This->fNavHistory[ilev] = tmpNavLevel;
  }
}

FQUALIFIER
void GPNavigationHistory_SetFirstEntry(GPNavigationHistory *This,
				       GEOMETRYLOC GPVPhysicalVolume* pVol)
{
  GPThreeVector translation = GPThreeVector_create(0.,0.,0.);
  G4int copyNo = -1;

  // Protection needed in case pVol=null 
  // so that a touchable-history can signal OutOfWorld 
  //

  if( pVol != GEOMETRYNULL )
  {
    translation = GPVPhysicalVolume_GetTranslation(pVol);
    //    copyNo = GPVPhysicalVolume_GetCopyNo(pVol);
    copyNo = 0;
  }

  GPNavigationLevel navLevel;

  GPAffineTransform aT;
  GPAffineTransform_Vector(&aT,translation);

  GPNavigationLevel_Constructor(&navLevel,
				pVol, aT, kNormal, copyNo );

  This->fNavHistory[0] = navLevel;

}

FQUALIFIER
const GPAffineTransform* 
GPNavigationHistory_GetPtrTopTransform(GPNavigationHistory *This)
{
  return GPNavigationLevel_GetPtrTransform(
	 &(This->fNavHistory[This->fStackDepth]));
}

FQUALIFIER
GPAffineTransform 
GPNavigationHistory_GetTopTransform(GPNavigationHistory *This)
{
  return  GPNavigationLevel_GetTransform(
	 &(This->fNavHistory[This->fStackDepth]));
}

FQUALIFIER
G4int GPNavigationHistory_GetTopReplicaNo(GPNavigationHistory *This)
{
  return  GPNavigationLevel_GetReplicaNo(
	 &(This->fNavHistory[This->fStackDepth]));
}

FQUALIFIER
EVolume GPNavigationHistory_GetTopVolumeType(GPNavigationHistory *This)
{
  return  GPNavigationLevel_GetVolumeType(
	 &(This->fNavHistory[This->fStackDepth]));
}

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigationHistory_GetTopVolume(GPNavigationHistory *This)
{
  return  GPNavigationLevel_GetPhysicalVolume(
	 &(This->fNavHistory[This->fStackDepth]));
}

FQUALIFIER
G4int GPNavigationHistory_GetDepth(GPNavigationHistory *This)
{
  return This->fStackDepth;
}

FQUALIFIER
GPAffineTransform GPNavigationHistory_GetTransform(GPNavigationHistory *This,
						   G4int n)
{
  return GPNavigationLevel_GetTransform(&(This->fNavHistory[n]));
}

FQUALIFIER
G4int GPNavigationHistory_GetReplicaNo(GPNavigationHistory *This,
				       G4int n)
{
  return GPNavigationLevel_GetReplicaNo(&(This->fNavHistory[n]));
}

FQUALIFIER
EVolume GPNavigationHistory_GetVolumeType(GPNavigationHistory *This,
					  G4int n)
{
  return GPNavigationLevel_GetVolumeType(&(This->fNavHistory[n]));
}

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigationHistory_GetVolume(GPNavigationHistory *This,
						 G4int n)
{
  return GPNavigationLevel_GetPhysicalVolume(&(This->fNavHistory[n]));
}

FQUALIFIER
G4int GPNavigationHistory_GetMaxDepth(GPNavigationHistory *This)
{
  //  return fNavHistory.size();
  return kHistoryMax;
}

FQUALIFIER
void GPNavigationHistory_BackLevel(GPNavigationHistory *This)
{
  //  assert( fStackDepth>0 );

  // Tell  the  level  that I am forgetting it
  // delete fNavHistory(fStackDepth);
  //
  (This->fStackDepth)--;
}

FQUALIFIER
void GPNavigationHistory_BackLevel2(GPNavigationHistory *This,
				   G4int n)
{
  //  assert( n<=fStackDepth );
  This->fStackDepth-=n;
}

/*
///@@@G4FWP fixed History Size
void GPNavigationHistory_EnlargeHistory(GPNavigationHistory *This)
{
  //  G4int len = fNavHistory.size();
  G4int len = kHistoryMax;
  if ( len==This->fStackDepth )
  {
    // Note: Resize operation clears additional entries
    //
    G4int nlen = len+kHistoryStride;
    fNavHistory.resize(nlen);
  }  
}
*/

FQUALIFIER
void GPNavigationHistory_NewLevel( GPNavigationHistory *This,
				   GEOMETRYLOC GPVPhysicalVolume *pNewMother,
				   EVolume vType,
				   G4int nReplica )
{
  This->fStackDepth++;
  //  EnlargeHistory();  // Enlarge if required

  GPNavigationLevel navLevel;

  GPAffineTransform aT;
  GPAffineTransform_Constructor3(&aT,
				 GPVPhysicalVolume_GetRotation(pNewMother),
                                 GPVPhysicalVolume_GetTranslation(pNewMother));

  GPAffineTransform bT = GPNavigationLevel_GetTransform(&(This->fNavHistory[This->fStackDepth-1]));

  GPNavigationLevel_Constructor2( &navLevel,
		     pNewMother, 
		     bT,
		     aT,
		     vType,
		     nReplica ); 

  This->fNavHistory[This->fStackDepth] = navLevel;

  // The constructor computes the new global->local transform
}
