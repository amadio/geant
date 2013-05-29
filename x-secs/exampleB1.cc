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
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "B1SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "QBBC.hh"
#include "G4HadronicProcess.hh"
#include "G4VEmProcess.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4ProductionCutsTable.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"

//---- ROOT stuff
#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new B1DetectorConstruction());

  // Physics list
  G4VModularPhysicsList* physicsList = new QBBC;
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);

  // Primary generator action
  runManager->SetUserAction(new B1PrimaryGeneratorAction());

  // Set user action classes
  //
  // Stepping action
  runManager->SetUserAction(new B1SteppingAction());     

  // Event action
  runManager->SetUserAction(new B1EventAction());

  // Run action
  runManager->SetUserAction(new B1RunAction());
     
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

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  // ---------------------------------------------------------------------
  //
  // Here x-sections should be defined... let's see...
  const G4VModularPhysicsList* pl = dynamic_cast<const G4VModularPhysicsList*>(runManager->GetUserPhysicsList());
  G4cout << pl << G4endl;
  //
  const G4int nbins = 1000;
  const G4double emin = 1.e-6*GeV;
  const G4double emax = 1.e3*GeV;
  const G4double delta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
  const G4double lemin = TMath::Log10(emin/GeV);
  const G4double lemax = TMath::Log10(emax/GeV);
  
  printf("%f %f \n",lemin, lemax);
  G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator *theParticleIterator = theParticleTable->GetIterator();
  theParticleIterator->reset();
  G4ParticleDefinition* particle;
  TFile *fh = new TFile("xsec.root","recreate");
  while((*theParticleIterator)()) {
    char string[1024]="\0";
    G4int ip=0;
    particle = theParticleIterator->value();
    sprintf(string,"%-20s:",(const char *)particle->GetParticleName());
    G4ProcessManager* pManager = particle->GetProcessManager();
    if ( pManager == 0) {
      G4cout << "No process manager" << G4endl;
    } else {
      G4ProcessVector* pList = pManager->GetProcessList();
      for (G4int idx=0; idx<pList->size(); idx++) {
	const G4VProcess* p = (*pList)[idx];
	//	if(p->GetProcessType()==1||p->GetProcessType()==6) continue;
	//	if(p->GetProcessType()==2) 
	//	  if(p->GetProcessSubType()==10||p->GetProcessSubType()==2) continue;
	sprintf(&string[strlen(string)]," [%s,%d,%d]",
		(const char *)p->GetProcessName(),p->GetProcessType(),
		p->GetProcessSubType());
	ip=1;
	if(!strcmp((const char*)particle->GetParticleName(),"pi+")&&
	   !strcmp((const char*)p->GetProcessName(),"pi+Inelastic")) {
	  printf("Found pi+ && p+Inelastic\n");
	  TH1F *fpiin = new TH1F("fpiin",(const char*)p->GetProcessName(),1000,lemin,lemax);
	  G4Material *iron =  G4Material::GetMaterial("G4_Fe");
	  G4double en = emin;
	  for(int i=0;i<1000;++i) {
	    fpiin->SetBinContent(i,((G4HadronicProcess*)p)->GetElementCrossSection(
	   new G4DynamicParticle(particle,G4ThreeVector(0,0,1),en),iron->GetElement(0))/barn);
	    en*=delta;
	    //	    printf("%f %f\n",TMath::Power(10,-3.+0.006*i)*GeV,fpiin->GetBinContent(i));
	  }
	  fpiin->Write();
	} else 	if(!strcmp((const char*)particle->GetParticleName(),"e-")&&
			   !strcmp((const char*)p->GetProcessName(),"eIoni")) {
	   printf("Found %s && %s\n",(const char*)particle->GetParticleName(),(const char*)p->GetProcessName());
	   TH1F *feioni = new TH1F("feioni",(const char*)p->GetProcessName(),1000,lemin,lemax);
	   G4Material *iron =  G4Material::GetMaterial("G4_Fe");

	   G4ProductionCutsTable* theCoupleTable =
	      G4ProductionCutsTable::GetProductionCutsTable();
	   
	   size_t numOfCouples = theCoupleTable->GetTableSize();
	   const G4MaterialCutsCouple* couple = 0;
	   for (size_t i=0; i<numOfCouples; i++) {
	      couple = theCoupleTable->GetMaterialCutsCouple(i);
	      if (couple->GetMaterial() == iron) break;
	   }
	   printf("Iron Z %f, A %f\n",iron->GetZ(), iron->GetA()*mole/g);
	   G4double natomscm3 = (Avogadro*iron->GetDensity()*cm3)/(iron->GetA()*mole);
	   printf("Density = %g, natoms/cm3 = %g\n",iron->GetDensity()*cm3/g,natomscm3);
	   G4VEnergyLossProcess *pt = (G4VEnergyLossProcess*)p;
	   G4double en=emin;
	   for(int i=0;i<1000;++i) {
	      G4double xsec = pt->CrossSectionPerVolume(en,couple)*cm/natomscm3/barn;
	      feioni->SetBinContent(i,xsec);
	      en*=delta;
	   }
	   feioni->Write();
	}
      }
    }
    if(ip) printf("%s\n",string);
  }
  fh->Write();
  fh->Close();

//---------------------------------------------------------------------------
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
