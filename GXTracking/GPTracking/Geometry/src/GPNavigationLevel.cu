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
// $Id: GPNavigationLevel.cc,v 1.5 2009-08-03 16:13:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 30.09.97 J.Apostolakis Initial version. 
//                    
// ----------------------------------------------------------------------

#include "GPNavigationLevel.h"

//G4Allocator<GPNavigationLevel> aNavigationLevelAllocator;

FQUALIFIER
void GPNavigationLevel_Constructor( GPNavigationLevel *This,
				    GEOMETRYLOC GPVPhysicalVolume* pPhysVol,
				    GPAffineTransform  afTransform,
				    EVolume            volTp,
				    G4int              repNo )
{
  //  GPNavigationLevelRep_Constructor(This->fLevelRep, 
  //				   pPhysVol, afTransform, volTp, repNo );
  This->sTransform = afTransform;
  This->sPhysicalVolumePtr = pPhysVol;

  This->sReplicaNo = repNo;
  This->sVolumeType = volTp;
  This->fCountRef = 1; 

}

FQUALIFIER
void GPNavigationLevel_Constructor2( GPNavigationLevel *This,
				     GEOMETRYLOC GPVPhysicalVolume *pPhysVol,
				     GPAffineTransform levelAbove,
				     GPAffineTransform relativeCurrent,
				     EVolume            volTp,
				     G4int              repNo )
{

  //  GPNavigationLevelRep_Constructor2(This->fLevelRep, 
  //				    pPhysVol, 
  //				    levelAbove, 
  //				    relativeCurrent, 
  //				    volTp, 
  //				    repNo );

  This->sPhysicalVolumePtr = pPhysVol;
  This->sReplicaNo = repNo;
  This->sVolumeType = volTp;
  This->fCountRef = 1; 
  GPAffineTransform_InverseProduct(&(This->sTransform),&levelAbove,
				   &relativeCurrent );

}

FQUALIFIER
void GPNavigationLevel_GPNavigationLevel0(GPNavigationLevel *This)
{
  //  GPNavigationLevelRep_Constructor0(This->fLevelRep);
  This->sPhysicalVolumePtr = NULL;
  This->sReplicaNo = -1;
  This->sVolumeType = kReplica;
  This->fCountRef = 1; 

}

/*
GPNavigationLevel_GPNavigationLevel(const GPNavigationLevel& right)
  : fLevelRep( right.fLevelRep )
{
  fLevelRep->AddAReference(); 
}

GPNavigationLevel_~GPNavigationLevel()
{
  if( fLevelRep->RemoveAReference() )  { delete fLevelRep; }
}

GPNavigationLevel& GPNavigationLevel_operator=(const GPNavigationLevel &right)
{ 
  if ( &right != this )
  {
    right.fLevelRep->AddAReference(); 
    if( fLevelRep->RemoveAReference() )  { delete fLevelRep; }
    fLevelRep = right.fLevelRep;
  }
  return *this;
} 
*/

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
// $Id: GPNavigationLevel.icc,v 1.24 2010-10-27 07:34:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 30.09.97 J.Apostolakis Initial version. 
//                        
// ----------------------------------------------------------------------

//#if defined G4GEOM_ALLOC_EXPORT
//extern G4DLLEXPORT G4Allocator<GPNavigationLevel> aNavigationLevelAllocator;
//#else
//extern G4DLLIMPORT G4Allocator<GPNavigationLevel> aNavigationLevelAllocator;
//#endif

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigationLevel_GetPhysicalVolume(GPNavigationLevel *This) 
{ 
  //  return GPNavigationLevelRep_GetPhysicalVolume(This->fLevelRep); 
  return This->sPhysicalVolumePtr;
}

FQUALIFIER
GPAffineTransform GPNavigationLevel_GetTransform(GPNavigationLevel *This) 
{ 
  //  return GPNavigationLevelRep_GetTransform(This->fLevelRep) ; 
  return This->sTransform;
} 

FQUALIFIER
GPAffineTransform* GPNavigationLevel_GetPtrTransform(GPNavigationLevel *This) 
{ 
  //  return GPNavigationLevelRep_GetTransformPtr(This->fLevelRep) ; 
  return &(This->sTransform);
} 

FQUALIFIER
EVolume GPNavigationLevel_GetVolumeType(GPNavigationLevel *This) 
{ 
  //  return GPNavigationLevelRep_GetVolumeType(This->fLevelRep) ; 
  return This->sVolumeType;
}

FQUALIFIER
G4int GPNavigationLevel_GetReplicaNo(GPNavigationLevel *This) 
{ 
  //  return GPNavigationLevelRep_GetReplicaNo(This->fLevelRep) ; 
  return This->sReplicaNo;
}

FQUALIFIER
void GPNavigationLevel_AddAReference(GPNavigationLevel *This) 
{
  (This->fCountRef)++; 
}

FQUALIFIER
G4bool GPNavigationLevel_RemoveAReference(GPNavigationLevel *This) 
{
  return( --(This->fCountRef) <= 0 ); 
}
