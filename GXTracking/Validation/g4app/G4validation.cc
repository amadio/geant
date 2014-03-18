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
#ifdef GPUNONRANDOM
  #warning GPUNONRANDOM is defined!
  #include "CLHEP/Random/NonRandomEngine.h"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsListBrem.hh"
#include "PhysicsListIoni.hh"
#include "PhysicsListMsc.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingActionBrem.hh"
#include "SteppingActionIoni.hh"
#include "SteppingActionMsc.hh"
#include "SteppingVerbose.hh"
#include "GPHistoManager.hh"

#include "cmsExpMagneticField.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include <iostream>
using std::cout;
using std::endl;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  G4int emModel = atoi(argv[2]);
  G4double fixedEnergy = 0;
  if(argc>=4) fixedEnergy = atof(argv[3]);
  // cout<<" emModel = "<< emModel
  //     <<" energy = "<< fixedEnergy
  //     <<" argv[0]="<< argv[0]  // command line
  //     <<" argv[1]="<< argv[1]
  //     <<" argv[2]="<< argv[2]
  //     <<" argv[3]="<< argv[3]
  //     << endl;

  //=======  Booking histograms =================
  GPHistoManager& hmgr = GPHistoManager::getInstance();

  hmgr.book1F("hstepBrem",100,0,10);
  hmgr.book1F("hstepIoni",100,0,1);
  hmgr.book1F("hstepMsc",100,0,10);

  if(emModel==1) hmgr.bookBremHistos();
  if(emModel==2) hmgr.bookIoniHistos();
  if(emModel==3) hmgr.bookMscHistos();

  //=======  end of histogram booking ============
  // Choose the Random engine
  //
#ifndef GPUNONRANDOM
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
#else
  CLHEP::NonRandomEngine* nonRandomEngine = new CLHEP::NonRandomEngine();
  nonRandomEngine->setNextRandom(0.123456);
  CLHEP::HepRandom::setTheEngine(nonRandomEngine);
  for(int i=0; i<10; ++i) {
    cout<<"Random: "<< nonRandomEngine->flat() << endl;
  }
#endif
  
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction(4,3,8.28*g/cm3));
  //
  if(emModel==1) runManager->SetUserInitialization(new PhysicsListBrem);
  if(emModel==2) runManager->SetUserInitialization(new PhysicsListIoni);
  if(emModel==3) runManager->SetUserInitialization(new PhysicsListMsc);

  // Set user action classes
  //
  G4int nParticlesPerEvent = 20000;
#ifdef GPUDEBUG
  nParticlesPerEvent = 2;
#endif
  PrimaryGeneratorAction* ppga = new PrimaryGeneratorAction(nParticlesPerEvent);
  ppga->SetEnergy(fixedEnergy);
  runManager->SetUserAction( ppga );
  //
  runManager->SetUserAction(new RunAction);
  //
  runManager->SetUserAction(new EventAction);
  //
  if(emModel==1) runManager->SetUserAction(new SteppingActionBrem);
  if(emModel==2) runManager->SetUserAction(new SteppingActionIoni);
  if(emModel==3) runManager->SetUserAction(new SteppingActionMsc);
  //
  runManager->SetUserAction(new TrackingAction);

  // Initialize G4 kernel
  //
  runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {  // interactive mode : define UI session
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

  hmgr.destroy();

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
