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
/// \file persistency/gdml/G01/include/G01DetectorConstruction.hh
/// \brief Definition of the G01DetectorConstruction class
//
//
// $Id: G01DetectorConstruction.hh 68025 2013-03-13 13:43:46Z gcosmo $
//
//
// Original Geant4 was calss: 
//  examples/extended/persistency/gdml/G01/include/G01DetectorConstruction.hh


#ifndef DETECTORCONSTRUCTION_H
#define DETECTORCONSTRUCTION_H 1

#include "G4VUserDetectorConstruction.hh"

class TGeoManager;
class G4UniformMagField;

/// Detector construction allowing to use the geometry read from the GDML file

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
 
    DetectorConstruction(G4VPhysicalVolume *setWorld = 0);
    ~DetectorConstruction();

    virtual G4VPhysicalVolume *Construct()
    {
      return fWorld;
    }

   void SetMagField(G4double);


  private:

    G4VPhysicalVolume  *fWorld;
    static TGeoManager *fgGeomMgrRoot;

    G4UniformMagField* magField;      //pointer to the magnetic field

};

#endif
