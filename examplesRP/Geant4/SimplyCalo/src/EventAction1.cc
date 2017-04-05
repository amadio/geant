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
/// \file electromagnetic/TestEm3/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 98762 2016-08-09 14:08:07Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction1.hh"

#include "Run1.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction1::EventAction1(DetectorConstruction* det)
:G4UserEventAction(),fDetector(det)
{
  fNbOfLayers = fDetector->GetNbOfLayers();
  fNbOfAbsor  = fDetector->GetNbOfAbsor()+1;
  fSumEdeps.resize(fNbOfLayers);
  fSumLengths.resize(fNbOfLayers);
  for (G4int i=0; i<fNbOfLayers; ++i) {
    fSumEdeps[i].resize(fNbOfAbsor,0.0);
    fSumLengths[i].resize(fNbOfAbsor,0.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction1::~EventAction1()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction1::BeginOfEventAction(const G4Event*)
{
  fNbOfLayers = fDetector->GetNbOfLayers();
  fNbOfAbsor  = fDetector->GetNbOfAbsor()+1;
  for (G4int i=0; i<fNbOfLayers; ++i) {
    for (G4int j=0; j<fNbOfAbsor; ++j) {
      fSumEdeps[i][j]   = 0.0;
      fSumLengths[i][j] = 0.0;
    }
  }
}

void EventAction1::FillPerStep(G4int layerNum, G4int absorNum, G4double edep, G4double stepl) {
  fSumEdeps[layerNum][absorNum]   += edep;
  fSumLengths[layerNum][absorNum] += stepl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction1::EndOfEventAction(const G4Event*)
{
  //get Run
  Run1* run = static_cast<Run1*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  for (G4int i=0; i<fNbOfLayers; ++i) {
    for (G4int j=0; j<fNbOfAbsor; ++j) {
      run->FillPerEvent(i,j,fSumEdeps[i][j],fSumLengths[i][j]);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
