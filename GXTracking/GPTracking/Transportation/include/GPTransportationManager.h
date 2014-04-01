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
// $Id: G4TransportationManager.hh,v 1.12 2007-04-20 15:28:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4TransportationManager
//
// Class description:
//
// A singleton class which stores the (volume) navigator used by 
// the transportation process to do the geometrical tracking.
// It also stores a pointer to the propagator used in a (magnetic) 
// field and to the field manager.
// The class instance is created before main() is called, and
// in turn creates the navigator and the rest.

// Created:  10 March 1997, J. Apostolakis
// Reviewed: 26 April 2006, G. Cosmo
// --------------------------------------------------------------------

#ifndef  GPTransportationManager_hh
#define  GPTransportationManager_hh

//#include "G4SafetyHelper.hh"
//#include <vector>

#include "GPNavigator.h"
#include "GPPropagatorInField.h"
#include "GPFieldManager.h"
#include "GPMultiLevelLocator.h"

//class G4GeometryMessenger;
//class G4VPhysicalVolume;

struct GPTransportationManager 
{
  //std::vector<G4Navigator*> fNavigators; //@@@G4FWP no G4Navigator yet
  // The collection of all navigators registered 
  //std::vector<G4Navigator*> fActiveNavigators; //@@@G4FWP G4Navigator yet
  // The collection of only active navigators
  //std::vector<G4VPhysicalVolume*> fWorlds; //@@@G4FWP G4Navigator yet
  // The collection of worlds associated to the registered navigators
  
  GPNavigator* trackingNavigator; //@@@G4FWP - temporary NULL pointer
  GPMultiLevelLocator* multiLevelLocatpr; //@@@G4FWP - temporary NULL pointer

  GPPropagatorInField*    fPropagatorInField;
  GPFieldManager*         fFieldManager;
  //  G4GeometryMessenger*    fGeomMessenger;
  //  G4SafetyHelper*         fSafetyHelper;

  //  static G4TransportationManager*  fTransportationManager;
};

extern "C" {

FQUALIFIER
void GPTransportationManager_Constructor(GPTransportationManager *This,
					 GPFieldMap *fieldMapArray) ;

FQUALIFIER
void GPTransportationManager_SetFieldManager( GPTransportationManager *This,
					      GPFieldManager* newFieldManager);

FQUALIFIER
GPFieldManager* GPTransportationManager_GetFieldManager( GPTransportationManager *This );

FQUALIFIER
void 
GPTransportationManager_SetPropagatorInField( GPTransportationManager *This,
					      GPPropagatorInField* newFieldPropagator );
FQUALIFIER
GPPropagatorInField* 
GPTransportationManager_GetPropagatorInField( GPTransportationManager *This );

}

#endif 
