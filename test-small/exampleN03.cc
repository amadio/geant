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

#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"

// incs. for setting the production cuts
#include "G4RegionStore.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4MaterialCutsCouple.hh"

// incs. for setting the tracking cuts
#include "G4UserSpecialCuts.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ProcessManager.hh"

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
void setG4ProductionCut();
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  G4bool isBachMode           = FALSE; // do we run in bach mode ?  ; 
  char mFileName[512]         = "";    // macro filename (in case of bach mode) 
  char physListName[512]      = "TABPHYS";  // physics list; default tabulated  
  G4int tabPhysVerboseLevel   = 0;      // verbosity level in tabulated physics 
  G4double tabPhysEnergyLimit = 3.0e-6; // low energy cut in tab. physics [GeV] 
  
  // parsing the arguments
  int c; 
  while ((c = getopt (argc, argv, "m:v:l:p:r")) != -1)
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
      case 'p':
        strncpy(physListName,optarg,strlen(optarg));
        physListName[strlen(optarg)] = '\0';
        break;
      case 'r':
        TabulatedDataManager::fgIsUseRange = TRUE;
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
  // set tracking cut energy value in case of G4 physics list 
  // (the energy value is the parameter given with -l flag in GeV)
  if(strcmp(physListName,"TABPHYS")) 
    DetectorConstruction::fTrackingCutInEnergy = tabPhysEnergyLimit*GeV; 

  runManager->SetUserInitialization(new DetectorConstruction);
  //

  // Set physics list
  if(!strcmp(physListName,"QBBC")) {
    runManager->SetUserInitialization(new QBBC);
  } else if(!strcmp(physListName,"FTFP_BERT")) {
    runManager->SetUserInitialization(new FTFP_BERT);
  } else if(!strcmp(physListName,"FTFP_BERT_HP")) {
    runManager->SetUserInitialization(new FTFP_BERT_HP);
  } else if(!strcmp(physListName,"TABPHYS")) {
    runManager->SetUserInitialization(new SimplePhysicsList);
    RunAction::isTabPhys = TRUE;  
  } else {
    G4cout << "Unknown physics list " << physListName << G4endl;
    exit(1);
  }

  G4cout << "Physics List:  " << physListName << G4endl;

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



  // set some optional parameters in tabulated physics in case of TABPHYS 
  if(!strcmp(physListName,"TABPHYS")) {
    TabulatedDataManager::SetVerboseLevel( tabPhysVerboseLevel );
    TotalPhysicsProcess::SetEnergyLimit( tabPhysEnergyLimit ); 
  }  

  // set production and tracking cuts in case of G4 physics lists 
  if(strcmp(physListName,"TABPHYS")) {
    // value is fixed now to 1 keV for consistency with the tabulated data
    setG4ProductionCut();  
    // setting the tracking cut energy for all registred particles 
    // (value has already been set in Det.Construction)
    G4ParticleDefinition* particle;
    G4ProcessManager *pm;
    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();     
    G4int np=theParticleTable->size();
    for(G4int i=0;i<np;++i)
     theParticleTable->GetParticle(i)->GetProcessManager()->AddProcess(new G4UserSpecialCuts(),-1,-1,1);     
  }

   

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
  "============================================================================"
  << G4endl <<
  "NAME" << G4endl <<
  "    exampleN03 ----------------- Geant4 example --------------------------- " 
       << G4endl << G4endl <<
  "SYNOPSIS" << G4endl <<
  "    exampleN03 [-l <value>, -v <value>, -m <FILE>, -p <NAME>, -r] "
        << G4endl << G4endl <<
  "DESCRIPTION" << G4endl <<
  "    Run Geant4 /examples/novice/N03/exampleN03 with optional physics list " 
        << G4endl << G4endl <<
  "OPTIONS" << G4endl <<
  "    -l low energy cut in tabulated physics [GeV] (default 3.0e-6)" 
          << G4endl <<
  "    -v verbosity level in tabulated physics (only >=2 is used at the moment)" 
          << G4endl <<
  "    -m <FILE> if we run it in bach mode; file is the Geant4 macro file" 
          << G4endl <<
  "    -p <NAME> physics list: TABPHYS (default), QBBC, FTFP_BERT, FTFP_BERT_HP"
          << G4endl <<
  "    -r use ranges, inverse ranges, continuous step limitation and proper way"
          << G4endl <<
  "       of along step energy loss computation with 1% linear loss limit"
          << G4endl <<
  "       in the tabulated physics. The pre-setp energy will be used to sample"
          << G4endl <<
  "       the type of the discrete interaction at the post-step point that "
          << G4endl <<
  "       corresponds to G4 without using integral approach. (This flag need to"
          << G4endl <<
  "       be used if the tracking cut value [-l <value>] < 3 keV)."
  << G4endl <<  
  "============================================================================"
  << G4endl << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Setting production cut (for gamma, e-, e+, p+) to a value given in energy.
// Will be done only if one of the original G4 physics list is used.
void setG4ProductionCut(){
  G4double ethresh = 1.0*keV; // production cut in energy to be set in range
   
  G4cout <<
  "============================================================================"
  << G4endl;
  printf("\t\tSetting Geant4 production cuts in energy START\n");
  G4MaterialTable *theMaterialTable = (G4MaterialTable*)G4Material::GetMaterialTable();
  //printf("Number of materials = %d\n",theMaterialTable->size());

  G4VRangeToEnergyConverter *converter[4];
  converter[0] = new G4RToEConvForGamma();
  converter[1] = new G4RToEConvForElectron();
  converter[2] = new G4RToEConvForPositron();
  converter[3] = new G4RToEConvForProton();
      
  G4RegionStore *theRegionStore = (G4RegionStore*) G4RegionStore::GetInstance();
      
  G4int nRegions = theRegionStore->size();
  printf("\t\tWe found %d regions:\n",nRegions);

  for(G4int i=0;i<nRegions;++i) {
    G4Region *reg = (*theRegionStore)[i];
    std::cout<<"\t\tName of "<< i <<"-th region is: "<< reg->GetName() << std::endl;

    for(G4int imat=0; imat<theMaterialTable->size(); ++imat) {
      G4MaterialCutsCouple *mcCouple = reg->FindCouple((*theMaterialTable)[imat]);
      if(mcCouple!=0) {
         G4ProductionCuts* pcuts = reg->GetProductionCuts();
         const char* parcuts[4]={"gamma","e-","e+","proton"};
         printf("Production Cuts (cm) : ");
         for(G4int ic=0; ic<4; ++ic) 
           printf("%s=%12.2g  [cm]  ",parcuts[ic],pcuts->GetProductionCut(ic)/cm);
         printf("\n");
         G4Material *mat = (*theMaterialTable)[imat];  
         printf("Production Cuts (GeV): ");
         for(G4int ic=0; ic<4; ++ic)
           printf("%s=%12.2g  [GeV]  ", parcuts[ic],
                converter[ic]->Convert(pcuts->GetProductionCuts()[ic],mat)/GeV);
         printf("\n");

         for(G4int ic=0; ic<4; ++ic) {
           G4double lmin = 1*nm;
           G4double lmax = 10*km;
           while(std::abs(lmin-lmax)>0.5*(lmin+lmax)*1e-10) {
             G4double lmid = 0.5*(lmin+lmax);
             if(converter[ic]->Convert(lmid,mat) > ethresh) lmax=lmid;
             else lmin = lmid;
           }
           reg->GetProductionCuts()->SetProductionCut(0.5*(lmin+lmax),ic);
         }

         //After setting
         pcuts = reg->GetProductionCuts();
         printf("Production Cuts (cm) : ");
         for(G4int ic=0; ic<4; ++ic) 
           printf("%s=%12.2g  [cm]  ",parcuts[ic],pcuts->GetProductionCut(ic)/cm);
         printf("\n");
         printf("Production Cuts (GeV): ");
         for(G4int ic=0; ic<4; ++ic)
           printf("%s=%12.2g  [GeV]  ", parcuts[ic],
         converter[ic]->Convert(pcuts->GetProductionCuts()[ic],mat)/GeV);
         printf("\n");
      }
    }
  }

  printf("\t\tSetting Geant4 production cuts in energy END\n");
  G4cout <<
  "============================================================================"
  << G4endl;

  // delete converters     
  for(G4int i=0; i<4; ++i) delete converter[i]; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

