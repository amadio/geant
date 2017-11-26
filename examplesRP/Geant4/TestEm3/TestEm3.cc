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
/// \file electromagnetic/TestEm3/TestEm3.cc
/// \brief Main program of the electromagnetic/TestEm3 example
//
// $Id: TestEm3.cc 92914 2015-09-21 15:00:48Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//
//   This is an extended version of the TestEm3 standard Geant4 test used for simplified calorimeter tests.
//   The orginal Geant4 TestEm3 application was taken from:
//   == Geant4 version Geant4.10.02.ref08
//   == examples/extended/electromagnetic/TestEm3
//


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


#include <getopt.h>
#include <err.h>
#include <iostream>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static bool isPerformance = false;
static std::string  macrofile="";

static struct option options[] = {
   {"flag to run the application in performance mode (default FALSE)", no_argument, 0, 'p'},
   {"standard Geant4 macro file", required_argument, 0, 'm'},
   {0, 0, 0, 0}
 };

void help() {
  std::cout<<"\n "<<std::setw(100)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  Extended version of the standard Geant4 TestEM3 simplified calorimeter test.    \n"
           <<"  The corresponding GV application is caloAppRP. \n"
           << std::endl
           <<"  **** Parameters: \n"
           <<"      -p :   flag  ==> run the application in performance mode i.e. no scoring \n"
           <<"         :   -     ==> run the application in NON performance mode i.e. with scoring (default) \n"
           <<"      -m :   REQUIRED : the standard Geamt4 macro file\n"
           << std::endl;

  std::cout<<"\nUsage: SimplyCalo [OPTIONS] INPUT_FILE\n\n"<<std::endl;

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  std::cout<<"\n "<<std::setw(100)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}




int main(int argc,char** argv) {
 //
 // process arguments
  if (argc == 1) {
    help();
    exit(0);
  }
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "pm:", options, &optidx);
    if (c == -1)
      break;
    //
    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'p':
      isPerformance = true;
      break;
    case 'm':
      macrofile = optarg;
      break;
    default:
      help();
      errx(1, "unknown option %c", c);
    }
  }
 //
 //
  //choose the Random engine: set to CLHEP::MTwistEngine => we use vecgeom RNG that is the standard c++11 mtwist
  G4Random::setTheEngine(new CLHEP::MTwistEngine);
  G4Random::setTheSeed(12345678);
  //G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);

  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  // Number of threads can be defined via from macro
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList);

  // set user action classes
  runManager->SetUserInitialization(new ActionInitialization(detector,isPerformance));

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

//  if (argc!=1) {  // batch mode
    G4String command = "/control/execute ";
    G4String fileName = macrofile;
    UI->ApplyCommand(command+fileName);
//  } else {          //define visualization and UI terminal for interactive mode
//#ifdef G4VIS_USE
//      G4VisManager* visManager = new G4VisExecutive;
//      visManager->Initialize();
//#endif

//#ifdef G4UI_USE
//      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
//      ui->SessionStart();
//      delete ui;
//#endif

//#ifdef G4VIS_USE
//      delete visManager;
//#endif
//  }

  // job termination
  //
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
