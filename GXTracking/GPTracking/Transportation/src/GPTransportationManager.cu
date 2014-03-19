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
// $Id: G4TransportationManager.cc,v 1.16 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4TransportationManager 
//
// Created : J.Apostolakis, 1997
// Reviewed: G.Cosmo, 2006
//  10.04.07 V.Ivanchenko  Use unique G4SafetyHelper
//
// --------------------------------------------------------------------

#include "GPTransportationManager.h"
#include "GPPropagatorInField.h"
#include "GPFieldManager.h"
#include "GPFieldMap.h"

//#include <algorithm>

//#include "G4GeometryMessenger.hh"
//#include "G4LogicalVolume.hh"
//#include "G4PVPlacement.hh"

// Initialise the static instance of the singleton
//
//G4TransportationManager* G4TransportationManager::fTransportationManager=0;

// ----------------------------------------------------------------------------
// Constructor
//

FQUALIFIER
void GPTransportationManager_Constructor(GPTransportationManager *This,
					 GPFieldMap *fieldMapArray) 
{ 
  //  if (fTransportationManager)
  //  {
  //    G4Exception("G4TransportationManager::G4TransportationManager()",
  //                "GeomNav0002", FatalException,
  //                "Only ONE instance of G4TransportationManager is allowed!");
  //  }

  // Create the navigator for tracking and activate it; add to collections
  //
  //  G4Navigator* trackingNavigator = new G4Navigator();
  //  trackingNavigator->Activate(true);
  //  fNavigators.push_back(trackingNavigator);
  //  fActiveNavigators.push_back(trackingNavigator);
  //  fWorlds.push_back(trackingNavigator->GetWorldVolume()); // NULL registered

  //  fGeomMessenger    = new G4GeometryMessenger(this);
  //  fSafetyHelper     = new G4SafetyHelper();

  //G4FWP: navigators (so geometry) are not available on the device yet
  //       only the field manager and the propagator are created

  //@@@G4FWP: all hooks to start transportation
  //1.construct magnetic field with the field map - this is the only input

  GPMagneticField magField;    
  GPMagneticField_Constructor(&magField,fieldMapArray);

  GPEquationOfMotion equaOfMotion;
  GPEquationOfMotion_Constructor(&equaOfMotion,&magField);

  GPClassicalRK4 rk4;
  GPClassicalRK4_Constructor(&rk4,&equaOfMotion,12);

  GPMagInt_Driver magIntDriver;
  GPMagInt_Driver_Constructor(&magIntDriver,1.0,&rk4,12,0);

  GPChordFinder chordFinder;
  GPChordFinder_Constructor(&chordFinder,&magIntDriver);	
  //GPChordFinder_Constructor2(&chordFinder,&magField,1.0,&rk4);	

  //2.construt the field manager and the propagator

  //fFieldManager     = new G4FieldManager();
  GPFieldManager aFieldManager;
  GPFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);
  This->fFieldManager = &aFieldManager;

  //fPropagatorInField=new G4PropagatorInField(trackingNavigator,fFieldManager);
  This->trackingNavigator = NULL;
  GPPropagatorInField propagatorInField;
  GPPropagatorInField_Constructor(&propagatorInField,0,
    				  &aFieldManager,0); 
  This->fPropagatorInField = &propagatorInField;

} 

// ----------------------------------------------------------------------------
// SetFieldManager(): Set the associated field manager.
//
FQUALIFIER
void GPTransportationManager_SetFieldManager( GPTransportationManager *This,
					      GPFieldManager* newFieldManager)
{
  This->fFieldManager = newFieldManager; 
  
  // Message the PropagatorInField, 
  // which also maintains this information (to be reviewed)
  //
  if( This->fPropagatorInField )
  {
    GPPropagatorInField_SetDetectorFieldManager( This->fPropagatorInField, 
						 newFieldManager );
  }
}

// ----------------------------------------------------------------------------
// GetFieldManager(): Return the associated field manager.
//
FQUALIFIER
GPFieldManager* GPTransportationManager_GetFieldManager( GPTransportationManager *This )
{
  return This->fFieldManager;
}

// ----------------------------------------------------------------------------
// SetPropagatorInField(): Set the associated propagator in field.
//

FQUALIFIER
void 
GPTransportationManager_SetPropagatorInField( GPTransportationManager *This,
					      GPPropagatorInField* newFieldPropagator )
{
  This->fPropagatorInField = newFieldPropagator;
}

// ----------------------------------------------------------------------------
// GetPropagatorInField(): Return the associated propagator in field.
//
FQUALIFIER
GPPropagatorInField* 
GPTransportationManager_GetPropagatorInField( GPTransportationManager *This )
{
  return This->fPropagatorInField;
}

