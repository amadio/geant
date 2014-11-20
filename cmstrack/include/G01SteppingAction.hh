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
// $Id$
//
/// \file G01SteppingAction.hh
/// \brief Definition of the G01SteppingAction class

#ifndef G01SteppingAction_h
#define G01SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Navigator;

/// Stepping action class
/// 
/// It holds data member fEnergy for accumulating the energy deposit
/// in a selected volume step by step.
/// The selected volume is set from  the detector construction via the  
/// SetVolume() function. The accumulated energy deposit is reset for each 
/// new event via the Reset() function from the event action.

class G01SteppingAction : public G4UserSteppingAction
{
  public:
    G01SteppingAction();
    virtual ~G01SteppingAction();

    // static access method
    static G01SteppingAction* Instance();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

    // reset accumulated energy
    void Reset();

    // get methods
    G4double GetEnergy() const { return fEnergy; }

    G4Navigator* GetNavigator() const {return fNavigator;}
    void SetNavigator(G4Navigator *nav) {fNavigator=nav;}
   
  private:
   G01SteppingAction(const G01SteppingAction&); // Not implemented
   G01SteppingAction& operator=(const G01SteppingAction&); // Not implemented
   
    static G01SteppingAction* fgInstance;  

    G4Navigator *fNavigator;
  
    G4double  fEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
