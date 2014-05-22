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

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "SimplePhysicsList.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "TabulatedDataManager.hh"
#include "TotalPhysicsProcess.hh"

#include <unistd.h>

void usage();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  G4bool isBachMode           = FALSE;  // do we run in bach mode ?  ; 
  char mFileName[512]         = "";     // macro filename (in case of bach mode) 
  G4int tabPhysVerboseLevel   = 0;      // verbosity level in tabulated physics 
  G4double tabPhysEnergyLimit = 3.0e-6; // low energy cut in tab. physics [GeV] 
  
  // parsing the arguments
  int c; 
  while ((c = getopt (argc, argv, "m:v:l:")) != -1)
    switch (c) {
      case 'm':
        isBachMode = TRUE;
        strncpy(mFileName,optarg,strlen(optarg));
        mFileName[strlen(optarg)] = '\0';
        break;
      case 'v':
        tabPhysVerboseLevel = atoi(optarg);
        break;
      case 'l':
        tabPhysEnergyLimit = atof(optarg);
        break;
      case '?':
        usage();
        return 1;
      default:
        abort ();
    }


  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction);
  //
  // runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new SimplePhysicsList);
  
  // Set user action classes
  //
  runManager->SetUserAction(new PrimaryGeneratorAction);
  //
  runManager->SetUserAction(new RunAction);
  //
  runManager->SetUserAction(new EventAction);
  //
  runManager->SetUserAction(new SteppingAction);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();

  // MaterialConverter::Instance()->CreateRootMaterials();
    
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // set some optional parameters in tabulated physics
  TabulatedDataManager::SetVerboseLevel( tabPhysVerboseLevel );
  TotalPhysicsProcess::SetEnergyLimit( tabPhysEnergyLimit ); 

  if (isBachMode) {  // batch mode
      G4String command = "/control/execute ";
      G4String fileName = mFileName;
      UImanager->ApplyCommand(command+fileName);
  } else {  // interactive mode : define UI session
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac"); 
#endif
      if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void usage()
{
  G4cout <<
  "NAME" << G4endl <<
  "    exampleN03 -- Geant4 with tabulated physics data " << G4endl << G4endl <<
  "SYNOPSIS" << G4endl <<
  "    exampleN03 [-l, -v, -m [file]] " << G4endl << G4endl <<
  "DESCRIPTION" << G4endl <<
  "    Run Geant4 /examples/novice/N03/exampleN03 with tabulated physics data " 
        << G4endl << G4endl <<
  "OPTIONS" << G4endl <<
  "    -l low energy cut in tabulated physics [GeV] (default 3.0e-6)" 
          << G4endl <<
  "    -v verbosity level tabulated physics (only >=2 is used at the moment) " 
          << G4endl <<
  "    -m [file] if we run it in bach mode where file is the Geant4 macro file" 
  << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....


