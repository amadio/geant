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
/// \file persistency/gdml/G01/load_gdml.cc
/// \brief Main program of the persistency/gdml/G01 example
//
//
// $Id: load_gdml.cc 68025 2013-03-13 13:43:46Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 - load_gdml
//
// --------------------------------------------------------------

#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "VTfileio.h"

#include <vector>

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4RegionStore.hh"

#include "G01PrimaryGeneratorAction.hh"
#include "G01DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G01EventAction.hh"
#include "G01SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4GDMLParser.hh"

int main(int argc,char **argv)
{
   G4cout << G4endl;
   G4cout << "Usage: load_gdml <intput_gdml_file:mandatory>"
          << " <output_gdml_file:optional>" << G4endl;
   G4cout << G4endl;

   if (argc<2)
   {
      G4cout << "Error! Mandatory input file is not specified!" << G4endl;
      G4cout << G4endl;
      return -1;
   }

   G4GDMLParser parser;

// Uncomment the following if wish to avoid names stripping
// parser.SetStripFlag(false);

   parser.Read(argv[1]);
   
   if (argc>4)
   {
      G4cout << "Error! Too many arguments!" << G4endl;
      G4cout << G4endl;
      return -1;
   }

   G4RunManager* runManager = new G4RunManager;

   runManager->SetUserInitialization(new G01DetectorConstruction(
                                     parser.GetWorldVolume()));
   //   runManager->SetUserInitialization(new G01PhysicsList);
   runManager->SetUserInitialization(new FTFP_BERT);
   runManager->SetUserAction(new G01PrimaryGeneratorAction);
   runManager->SetUserAction(new G01EventAction);
   runManager->SetUserAction(new G01SteppingAction);

   runManager->Initialize();

   if (argc>=3)
   {
      parser.Write(argv[2], G4TransportationManager::GetTransportationManager()->
                   GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
   }

   G4UImanager* UImanager = G4UImanager::GetUIpointer();
 
   ///////////////////////////////////////////////////////////////////////
   //
   // Example how to retrieve Auxiliary Information
   //
   const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
   std::vector<G4LogicalVolume*>::const_iterator lvciter;

   VTfileio * io = VTfileio::I();

   io->OutFile(new TFile("cmstrack.root","recreate"));

   TObjArray *lv = io->GetVolumeDictionary();

   G4Region *region = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");
   G4ProductionCuts *wcuts = new G4ProductionCuts;
   wcuts->SetProductionCut(1*cm); // same cuts for gamma, e- and e+
   region->SetProductionCuts(wcuts);

   const int nlv = lvs->size();
   G4ProductionCuts *rcuts = new G4ProductionCuts[nlv];

   for( int ilv=0; ilv<nlv; ++ilv) 
   {
     G4LogicalVolume *glv = (*lvs)[ilv];
     lv->Add(new TObjString(glv->GetName()));
     double radl = glv->GetMaterial()->GetRadlen();
     region = new G4Region(G4String("Region_")+glv->GetName());
     region->AddRootLogicalVolume(glv);
     if(ilv==1623) rcuts[ilv].SetProductionCut(0.5*radl);
     else rcuts[ilv].SetProductionCut(0.25*radl);
     region->SetProductionCuts(&rcuts[ilv]);
    
     G4GDMLAuxListType auxInfo = parser.GetVolumeAuxiliaryInformation(glv);
     std::vector<G4GDMLAuxPairType>::const_iterator ipair = auxInfo.begin();
     for( ipair = auxInfo.begin(); ipair != auxInfo.end(); ipair++ )
     {
       G4String str=ipair->type;
       G4String val=ipair->value;
       G4cout << " Auxiliary Information is found for Logical Volume :  "
              << glv->GetName() << G4endl;
       G4cout << " Name of Auxiliary type is     :  " << str << G4endl;
       G4cout << " Associated Auxiliary value is :  " << val << G4endl;
     }
   }
   
   io->WriteDictionaries();

   //
   // End of Auxiliary Information block
   //
   ////////////////////////////////////////////////////////////////////////

   if (argc==4)   // batch mode  
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[3];
     UImanager->ApplyCommand(command+fileName);
   }
   else           // interactive mode
   {
#ifdef G4UI_USE
     G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
     G4VisManager* visManager = new G4VisExecutive;
     visManager->Initialize();
     UImanager->ApplyCommand("/control/execute vis.mac");
#endif
     ui->SessionStart();
#ifdef G4VIS_USE
     delete visManager;
#endif
     delete ui;
#endif
   }

   delete runManager;
   delete io;

   return 0;
}
