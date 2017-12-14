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
// $Id$
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "SteppingAnalysis.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
//#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "TPartIndex.h"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
    : detector((DetectorConstruction *)G4RunManager::GetRunManager()->GetUserDetectorConstruction()),
      eventaction((EventAction *)G4RunManager::GetRunManager()->GetUserEventAction()),
      doAnalysis(getenv("DO_STEPPING_ANALYSIS")), fSteppingAnalysis(0) {
  if (doAnalysis)
    fSteppingAnalysis = new SteppingAnalysis();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {
  if (doAnalysis)
    delete fSteppingAnalysis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction
#ifdef MAKESTAT
                  (const G4Step *aStep) {
  // get volume of the current step
  G4VPhysicalVolume *volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();

  G4double stepl = 0.;
  /*  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();

    if (volume == detector->GetAbsorber()) eventaction->AddAbs(edep,stepl);
    if (volume == detector->GetGap())      eventaction->AddGap(edep,stepl);
  */

  // determine which process happend
  G4int gVprocIndex = -1;
  G4VProcess const *g4proc = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if (!RunAction::isTabPhys) { // running with G4 phys. list. convert G4ProcName to GV
    if (G4String("UserSpecialCut") != g4proc->GetProcessName()) {
      G4int g4procCode = g4proc->GetProcessType() * 1000 + g4proc->GetProcessSubType();
      gVprocIndex = TPartIndex::I()->ProcIndex(g4procCode);
    } else {
      gVprocIndex = 18; // max -> userCut
    }
  } else { // running with tabulated phys.
    if (G4String("UserSpecialCut") != g4proc->GetProcessName()) {
      if (G4String("Transportation") == g4proc->GetProcessName())
        gVprocIndex = 0;
      else
        gVprocIndex = TPartIndex::I()->ProcIndex(g4proc->GetProcessName());
    } else {
      gVprocIndex = 18; // max -> userCut
    }
  }

  if (gVprocIndex < 0) {
    std::cout << "Unknown process [ " << g4proc->GetProcessName() << " ] in SteppingAction::UserSteppingAction !"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  //...
  stepl = aStep->GetStepLength() / CLHEP::cm; // change to cm
                                              // edep  /=GeV;                      // change to GeV
  if (volume == detector->GetGap()) {
    G4double pos = aStep->GetPreStepPoint()->GetPosition().x() / CLHEP::mm;
    G4int layer = (pos + 75.) / 15.;
    eventaction->FillPerStep(1, layer, edep, stepl, gVprocIndex);
  }
  if (volume == detector->GetAbsorber()) {
    G4double pos = aStep->GetPreStepPoint()->GetPosition().x() / CLHEP::mm;
    G4int layer = (pos + 75.) / 15.;
    eventaction->FillPerStep(0, layer, edep, stepl, gVprocIndex);
  }

  eventaction->AddOneStep();
#else
         (const G4Step *) {
#endif

  // example of saving random number seed of this event, under condition
  //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();

  //   Don't do !!
  //  if(doAnalysis) fSteppingAnalysis->DoIt(aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
