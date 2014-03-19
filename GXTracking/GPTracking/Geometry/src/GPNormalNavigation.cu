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
// $Id: GPNormalNavigation.cc,v 1.11 2010-11-04 08:57:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class GPNormalNavigation Implementation
//
// Author: P.Kent, 1996
//
// --------------------------------------------------------------------

#include "GPNormalNavigation.h"
#include "GPAffineTransform.h"
#include "GPAuxiliaryNavServices.h"

// ********************************************************************
// Constructor
// ********************************************************************
//
FQUALIFIER
void GPNormalNavigation_GPNormalNavigation(GPNormalNavigation *This)
{
  This->fCheck = false;
  //  fLogger = new G4NavigationLogger("GPNormalNavigation");
}

// ********************************************************************
// Destructor
// ********************************************************************
//
//GPNormalNavigation_~GPNormalNavigation()
//{
//  delete fLogger;
//}

// ********************************************************************
// ComputeStep
// ********************************************************************
//
FQUALIFIER
G4double
GPNormalNavigation_ComputeStep( GPNormalNavigation *This,
				GPThreeVector localPoint,
                                GPThreeVector localDirection,
                                G4double currentProposedStepLength,
				G4double *newSafety,
				GPNavigationHistory *history,
				G4bool *validExitNormal,
				GPThreeVector *exitNormal,
				G4bool *exiting,
				G4bool *entering,
				GEOMETRYLOC GPVPhysicalVolume *(*pBlockedPhysical),
				G4int *blockedReplicaNo)
{
  GEOMETRYLOC GPVPhysicalVolume *motherPhysical, *samplePhysical, *blockedExitedVol=0;
  GEOMETRYLOC GPLogicalVolume *motherLogical;
  GEOMETRYLOC GPVSolid *motherSolid;
  GPThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;

  motherPhysical = GPNavigationHistory_GetTopVolume(history);
  motherLogical  = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid    = GPLogicalVolume_GetSolid(motherLogical);

  // Compute mother safety
  //
  motherSafety = GPVSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety; // Working isotropic safety
  
  //
  // Compute daughter safeties & intersections
  //

  // Exiting normal optimisation
  //
  if ( *exiting && *validExitNormal )
  {
    if ( GPThreeVector_dot(localDirection,*exitNormal)>=kMinExitingNormalCosine )
    {
      // Block exited daughter volume
      //
      blockedExitedVol =* pBlockedPhysical;
      ourSafety = 0;
    }
  }
  *exiting  = false;
  *entering = false;

  localNoDaughters = GPLogicalVolume_GetNoDaughters(motherLogical);

  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo--)
  {
    samplePhysical = GPLogicalVolume_GetDaughter(motherLogical,sampleNo);
    if ( samplePhysical != blockedExitedVol )
    {

      GPAffineTransform sampleTf;
      GPAffineTransform_Constructor3(&sampleTf,
		        GPVPhysicalVolume_GetRotation(samplePhysical),
		        GPVPhysicalVolume_GetTranslation(samplePhysical));

      GPAffineTransform_Invert(&sampleTf);

      GPThreeVector samplePoint =
	GPAffineTransform_TransformPoint(&sampleTf,localPoint);
      GEOMETRYLOC GPVSolid *sampleSolid =
	GPLogicalVolume_GetSolid(GPVPhysicalVolume_GetLogicalVolume(samplePhysical));
 
      G4double sampleSafety =
	 GPVSolid_DistanceToIn(sampleSolid,samplePoint);

      if ( sampleSafety<ourSafety )
      {
        ourSafety=sampleSafety;
      }

      if ( sampleSafety<=ourStep )
      {
        sampleDirection = GPAffineTransform_TransformAxis(&sampleTf,localDirection);
  
      G4double sampleStep =
	  GPVSolid_DistanceToIn2(sampleSolid,samplePoint,sampleDirection);

        if ( sampleStep<=ourStep )
        {
          ourStep  = sampleStep;
          *entering = true;
          *exiting  = false;
          *pBlockedPhysical = samplePhysical;
          *blockedReplicaNo  = -1;
        }
      }
    }
  }

  if ( currentProposedStepLength<ourSafety )
  {
    // Guaranteed physics limited
    //
    *entering = false;
    *exiting  = false;
    *pBlockedPhysical = GEOMETRYNULL;
    ourStep = kInfinity;
  }
  else
  {
    // Compute mother intersection if required
    //
    if ( motherSafety<=ourStep )
    {
      G4double motherStep = GPVSolid_DistanceToOut2(
                            motherSolid,
			    localPoint,
			    localDirection,
			    true,
			    validExitNormal,
			    exitNormal);
      if ( motherStep<=ourStep )
      {
        ourStep  = motherStep;
        *exiting  = true;
        *entering = false;
        if ( *validExitNormal )
        {
          GPRotationMatrix *rot = GPVPhysicalVolume_GetRotation(motherPhysical);
          if (rot)
          {
	    //            exitNormal *= rot->inverse();
	    GPRotationMatrix inv = GPRotationMatrix_inverse(rot);
            *exitNormal = GPRotationMatrix_apply(&inv, *exitNormal);
          }
        }
      }
      else
      {
        *validExitNormal = false;
      }
    }
  }
  *newSafety = ourSafety;

  return ourStep;
}

// ********************************************************************
// ComputeSafety
// ********************************************************************
//
FQUALIFIER
G4double GPNormalNavigation_ComputeSafety(GPNormalNavigation *This,
					  GPThreeVector localPoint,
					  GPNavigationHistory *history,
					  const G4double)
{
  GEOMETRYLOC GPVPhysicalVolume *motherPhysical, *samplePhysical;
  GEOMETRYLOC GPLogicalVolume *motherLogical;
  GEOMETRYLOC GPVSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;

  motherPhysical = GPNavigationHistory_GetTopVolume(history);
  motherLogical  = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid    = GPLogicalVolume_GetSolid(motherLogical);

  // Compute mother safety
  //
  motherSafety = GPVSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety; // Working isotropic safety

  // Compute daughter safeties 
  //
  localNoDaughters = GPLogicalVolume_GetNoDaughters(motherLogical);
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical =  GPLogicalVolume_GetDaughter(motherLogical,sampleNo);

    GPAffineTransform sampleTf;
    GPAffineTransform_Constructor3(&sampleTf,
				   GPVPhysicalVolume_GetRotation(samplePhysical),
				   GPVPhysicalVolume_GetTranslation(samplePhysical));
    
    GPAffineTransform_Invert(&sampleTf);

    const GPThreeVector samplePoint =
      GPAffineTransform_TransformPoint(&sampleTf,localPoint);
    GEOMETRYLOC const GPVSolid *sampleSolid =
      GPLogicalVolume_GetSolid(GPVPhysicalVolume_GetLogicalVolume(samplePhysical));
    const G4double sampleSafety =
      GPVSolid_DistanceToIn(sampleSolid,samplePoint);
    if ( sampleSafety<ourSafety )
    {
      ourSafety = sampleSafety;
    }
  }
  return ourSafety;
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
// $Id: GPNormalNavigation.icc,v 1.5 2010-11-04 08:57:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GPNormalNavigation Inline Implementation
//
// --------------------------------------------------------------------

// ********************************************************************
// LevelLocate
// ********************************************************************
//
FQUALIFIER
G4bool
GPNormalNavigation_LevelLocate( GPNormalNavigation *This,
				GPNavigationHistory *history,
				GEOMETRYLOC GPVPhysicalVolume* blockedVol,
				// G4int dummy,
				GPThreeVector  globalPoint,
				GPThreeVector* globalDirection,
				G4bool  pLocatedOnEdge, 
				GPThreeVector* localPoint )
{
  GEOMETRYLOC GPVPhysicalVolume *targetPhysical, *samplePhysical;
  GEOMETRYLOC GPLogicalVolume *targetLogical;
  GEOMETRYLOC GPVSolid *sampleSolid;
  GPThreeVector samplePoint;
  G4int targetNoDaughters;
  
  targetPhysical = GPNavigationHistory_GetTopVolume(history);
  targetLogical = GPVPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetNoDaughters = GPLogicalVolume_GetNoDaughters(targetLogical);

  G4bool found = false;

  if (targetNoDaughters!=0)
  {
    //
    // Search daughters in volume
    //
    for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
    {
      samplePhysical = GPLogicalVolume_GetDaughter(targetLogical,sampleNo);
      if ( samplePhysical!=blockedVol )
      {
        // Setup history
        //
	GPNavigationHistory_NewLevel(history, samplePhysical, kNormal,0);  
	//		    GPVPhysicalVolume_GetCopyNo(samplePhysical) );

        sampleSolid = GPLogicalVolume_GetSolid(
                      GPVPhysicalVolume_GetLogicalVolume(samplePhysical));

	GPAffineTransform aT = GPNavigationHistory_GetTopTransform(history);

        samplePoint = GPAffineTransform_TransformPoint(&aT, globalPoint);

        if( GPAuxiliaryNavServices_CheckPointOnSurface(sampleSolid, 
						       samplePoint, 
						       globalDirection, 
						       aT, pLocatedOnEdge) )
        {
          // Enter this daughter
          //
          *localPoint = samplePoint;
          found = true;
          break;
        }
        else
        {
          GPNavigationHistory_BackLevel(history);
	}
      }
    }
  }

  return found;
}

// ********************************************************************
// GetVerboseLevel
// ********************************************************************
//
FQUALIFIER
G4int GPNormalNavigation_GetVerboseLevel(GPNormalNavigation *This)
{
  return 0;
  //  return fLogger->GetVerboseLevel();
}

// ********************************************************************
// SetVerboseLevel
// ********************************************************************
//
/*
FQUALIFIER
void GPNormalNavigation_SetVerboseLevel(GPNormalNavigation *This,
					G4int level)
{
  ;
  //  fLogger->SetVerboseLevel(level);
}
*/

// ********************************************************************
// CheckMode
// ********************************************************************
//
FQUALIFIER
void GPNormalNavigation_CheckMode(GPNormalNavigation *This, 
				  G4bool mode)
{
  This->fCheck = mode;
}
