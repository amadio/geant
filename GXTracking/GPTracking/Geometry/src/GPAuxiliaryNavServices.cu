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
// $Id: G4AuxiliaryNavServices.icc,v 1.4 2007-05-22 07:48:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4AuxiliaryNavServices __Host__ __Device__ implementation
//
// --------------------------------------------------------------------

#include "GPAuxiliaryNavServices.h"

FQUALIFIER 
G4bool
GPAuxiliaryNavServices_CheckPointOnSurface( GEOMETRYLOC const GPVSolid* sampleSolid, 
                     GPThreeVector localPoint, 
                     GPThreeVector *globalDirection, 
                     GPAffineTransform sampleTransform,
                     const G4bool locatedOnEdge)
{
  GPThreeVector localDirection, sampleNormal;
  G4bool        enter = false;

  EInside insideSolid = GPVSolid_Inside(sampleSolid,localPoint); 
  
  if ( insideSolid!=kOutside ) 
  {
    G4bool checkDirection= locatedOnEdge && (globalDirection!=0);
    if( (insideSolid==kSurface) && checkDirection)
    {
      // We are probably located on an edge.
      //
      localDirection= GPAffineTransform_TransformAxis(&sampleTransform,*globalDirection); 

      // Check whether we enter the volume
      // 
      sampleNormal = GPVSolid_SurfaceNormal(sampleSolid,localPoint);

      if ( GPThreeVector_dot(sampleNormal,localDirection) <= 0 )
      {
        if( GPThreeVector_dot(sampleNormal,localDirection) == 0 )
        {
          // We can't decide yet, let's make sure we're entering the solid.
          // If by a confusion we entered the next solid we find out now
          // whether to leave or to enter.
          // This happens when we're on the surface or edge shared by two
          // solids
          //
          G4double distanceToIn =
	    GPVSolid_DistanceToIn2(sampleSolid, localPoint, localDirection );
          if( distanceToIn != kInfinity )
          {
            enter = true;
          } 
        }
        else
        {
          enter = true;
        }
      }
    }
    else
    {
      enter = true;
    }
  }
  return enter;
}

// --------------------------------------------------------------------

FQUALIFIER 
G4bool
GPAuxiliaryNavServices_CheckPointExiting( GEOMETRYLOC const GPVSolid* sampleSolid, 
                   const GPThreeVector localPoint, 
                   const GPThreeVector* globalDirection, 
                   const GPAffineTransform sampleTransform )
{
  if( !globalDirection )  { return false; }

  GPThreeVector localDirection, sampleNormal;
  G4bool        exiting = false;

  EInside insideSolid = GPVSolid_Inside(sampleSolid,localPoint); 
  if( insideSolid==kSurface )
  {
    localDirection= GPAffineTransform_TransformAxis(&sampleTransform,*globalDirection); 

    // Check whether we are exiting the volume
    // 
    sampleNormal = GPVSolid_SurfaceNormal(sampleSolid,localPoint);
    if ( GPThreeVector_dot(sampleNormal,localDirection) >= 0 )
    {
      if(  GPThreeVector_dot(sampleNormal,localDirection) == 0 )
      {
        // We can't decide yet, let's make sure we're entering the solid.
        // If by a confusion we entered the next solid we find out now
        // whether to leave or to exiting.
        // This happens when we're on the surface or edge shared by two
        // solids
        //
        G4double distanceToIn =
	  GPVSolid_DistanceToIn2(sampleSolid, localPoint, localDirection );
        if( distanceToIn != kInfinity )
        {
          exiting = true;
        } 
      }
      else
      {
        exiting = true;
      }
    }
  }
  return exiting;
}
