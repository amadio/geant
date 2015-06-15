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
#include "TGeoManager.h"
#include "TFile.h"
#include "VTfileio.h"

#include <vector>
#include <sys/stat.h>

#include "G4RunManager.hh"
#include "G4FieldManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G01PrimaryGeneratorAction.hh"
#include "G01DetectorConstruction.hh"
#include "G01Field.hh"
#include "FTFP_BERT.hh"
#include "G01EventAction.hh"
#include "G01SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "MyGDMLWriteStructure.hh"

int main(int argc, char **argv) {
  using CLHEP::s;
  struct stat sb;

  G4cout << G4endl;
  G4cout << "Usage: load_gdml <intput_gdml_file:mandatory>"
         << " <output_gdml_file:optional>" << G4endl;
  G4cout << G4endl;

  if (argc < 2) {
    G4cout << "Error: mandatory input file is not specified!" << G4endl;
    G4cout << G4endl;
    return -1;
  }

  G4GDMLParser parser(new G4GDMLReadStructure(), new MyGDMLWriteStructure());
  //   parser.SetStripFlag(false);

  // Uncomment the following if wish to avoid names stripping
  // parser.SetStripFlag(false);

  parser.Read(argv[1], FALSE);

  if (argc > 4) {
    G4cout << "Error! Too many arguments!" << G4endl;
    G4cout << G4endl;
    return -1;
  }

  G4RunManager *runManager = new G4RunManager;

  runManager->SetUserInitialization(new G01DetectorConstruction(parser.GetWorldVolume()));

  //   runManager->SetUserInitialization(new G01PhysicsList);
  runManager->SetUserInitialization(new FTFP_BERT);
  runManager->SetUserAction(new G01PrimaryGeneratorAction);
  runManager->SetUserAction(new G01EventAction);
  runManager->SetUserAction(new G01SteppingAction);

  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------

  G01Field *myField = new G01Field;
  G4FieldManager *fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(myField);
  fieldMgr->CreateChordFinder(myField);

  runManager->Initialize();

  G4VPhysicalVolume *world =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();

  G4Navigator *secondNavigator = new G4Navigator();
  secondNavigator->SetWorldVolume(world);
  G01SteppingAction::Instance()->SetNavigator(secondNavigator);

  if (!stat("random_start.rand", &sb)) {
    printf("Setting the random number from file\n");
    G4Random::restoreEngineStatus("random_start.rand");
  }

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  ///////////////////////////////////////////////////////////////////////
  //
  // Example how to retrieve Auxiliary Information
  //
  const G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume *>::const_iterator lvciter;

  VTfileio *io = VTfileio::I();

  io->OutFile(new TFile("cmstrack.root", "recreate"));

  // Create root file
  TGeoManager::Import(argv[1]);
  gGeoManager->Write(gGeoManager->GetName());

  G4Region *region = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");
  G4ProductionCuts *wcuts = new G4ProductionCuts;
  wcuts->SetProductionCut(1 * CLHEP::cm); // same cuts for gamma, e- and e+
  region->SetProductionCuts(wcuts);

  const int nlv = lvs->size();
  G4ProductionCuts *rcuts = new G4ProductionCuts[nlv];
  std::vector<G4double> dcuts(4);

  printf("Setting cuts for %d logical volumes\n", nlv);

  G4double tlim = 1e-6 * s;
  char regname[14];
  for (int ilv = 0; ilv < nlv; ++ilv) {
    G4LogicalVolume *glv = (*lvs)[ilv];
    io->AddVolume(glv->GetName());
    io->AddShape(glv->GetSolid()->GetEntityType().c_str());
    double radl = glv->GetMaterial()->GetRadlen();
    snprintf(regname, 13, "Region_%4.4d_", ilv);
    region = new G4Region(G4String(regname) + glv->GetName());
    region->AddRootLogicalVolume(glv);
    G4double cfact = 0.25;
    // Lead Crystals are parametrised in reality
    //     if(ilv==1623) cfact=2;
    // This cone of copper very far needs high cuts
    // if(ilv==2891) cfact=4;
    dcuts[0] = 4 * cfact * radl;
    dcuts[1] = cfact * radl;
    dcuts[2] = cfact * radl;
    dcuts[3] = cfact * radl;
    rcuts[ilv].SetProductionCuts(dcuts);
    region->SetProductionCuts(&rcuts[ilv]);
    glv->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, tlim, 0., cfact * radl));

    G4GDMLAuxListType auxInfo = parser.GetVolumeAuxiliaryInformation(glv);
    if (!auxInfo.empty()) {
      //      G4cout << " Auxiliary Information is found for Logical Volume :  "
      //         << glv->GetName() << G4endl;
      //  std::vector<G4GDMLAuxPairType>::const_iterator ipair = auxInfo.begin();
      //  for (ipair = auxInfo.begin(); ipair != auxInfo.end(); ipair++) {
      //    G4String str = ipair->type;
      //    G4String val = ipair->value;
      //    G4cout << str << " = " << val << G4endl;
      // }
    }
  }

  io->WriteDictionaries();

#if NAVTEST
// check overlaps
//   const G4PhysicalVolumeStore *pvs = G4PhysicalVolumeStore::GetInstance();
// for (G4PhysicalVolumeStore::const_iterator it = pvs->begin() ; it !=
// pvs->end(); ++it) {
//   (*it)->CheckOverlaps(1000,1e-3,FALSE,1);
//}
#endif
  //
  // End of Auxiliary Information block
  //
  ////////////////////////////////////////////////////////////////////////

  if (argc >= 3) {
    parser.Write(argv[2], world->GetLogicalVolume());
  }

  if (argc == 4) // batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[3];
    UImanager->ApplyCommand(command + fileName);
  } else // interactive mode
  {
#ifdef G4UI_USE
    G4UIExecutive *ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    G4VisManager *visManager = new G4VisExecutive;
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

  G4Random::saveEngineStatus("random_end.rand");

  parser.Write("finalgdml.gdml", G4TransportationManager::GetTransportationManager()
                                     ->GetNavigatorForTracking()
                                     ->GetWorldVolume()
                                     ->GetLogicalVolume());

  delete runManager;
  delete io;

  return 0;
}
