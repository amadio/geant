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
/// \file electromagnetic/TestEm3/src/RunAction1.cc
/// \brief Implementation of the RunAction1 class
//
// $Id: RunAction1.cc 90969 2015-06-12 08:12:57Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction1.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
//#include "RunAction1Messenger.hh"
//#include "HistoManager.hh"
#include "Run1.hh"
#include "G4Timer.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction1::RunAction1(DetectorConstruction* det)
:G4UserRunAction(), fDetector(det), fRun(0), fTimer(0)
{
//  fRunMessenger = new RunAction1Messenger(this);
//  fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction1::~RunAction1()
{
//  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction1::GenerateRun()
{
  fRun = new Run1(fDetector);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction1::BeginOfRunAction(const G4Run*)
{
  if (isMaster) {
    //    G4Random::showEngineStatus();
    fTimer = new G4Timer();
    fTimer->Start();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction1::EndOfRunAction(const G4Run*)
{
  // compute and print statistic
  if (isMaster) {
    fTimer->Stop();
    if(!((G4RunManager::GetRunManager()->GetRunManagerType() ==
          G4RunManager::sequentialRM))) {
      G4cout << "\n" << "Total number of events:  "
             << fRun->GetNumberOfEvent() << G4endl;
      G4cout << "Master thread time:  "  << *fTimer << G4endl;
    }
    delete fTimer;
    fRun->EndOfRun();
  }
  //save histograms
//  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
//  if (analysis->IsActive()) {
//    analysis->Write();
//    analysis->CloseFile();
//  }

  // show Rndm status
  //  if (isMaster)  G4Random::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
