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
#include "B1materials.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "QBBC.hh"
#include "G4HadronicProcess.hh"
#include "G4VEmProcess.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialTable.hh"

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
#include <TMXsec.h>
#include <TRandom.h>

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

  // Set user action classes4
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
    //
    // Definitions
    //
    G4ParticleDefinition* particle;

    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();

    G4int np=theParticleTable->size();

    printf("Found %d particles\n",np);

    //
    const G4int maxproc = 10;  // max of 10 proc per particle
    const G4int nbins = 1000;
    const G4double emin = 1.e-6*GeV;
    const G4double emax = 1.e3*GeV;
    const G4double delta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
    const G4double lemin = TMath::Log10(emin/GeV);
    const G4double lemax = TMath::Log10(emax/GeV);
    
    printf("Log emin %f Log emax %f \n",lemin, lemax);
    

    TFile *fh = new TFile("xsec.root","recreate");
    /*
      snippets

      const G4VModularPhysicsList* pl = dynamic_cast<const G4VModularPhysicsList*>(runManager->GetUserPhysicsList());
    
      G4MaterialTable *theMaterialTable = (G4MaterialTable*)G4Material::GetMaterialTable();
      G4ElementTable  *theElementTable = (G4ElementTable*)G4Element::GetElementTable();
      
      for(i=0;i<matTable->size();i++)
      { delete (*(matTable))[i]; }
      matTable->clear();
      for(i=0;i<elemTable->size();i++)
      { delete (*(elemTable))[i]; }
      elemTable->clear();
      }    
      
      G4ParticleTable::G4PTblDicIterator *theParticleIterator = theParticleTable->GetIterator();
      theParticleIterator->reset();
      while((*theParticleIterator)()) {
      particle = theParticleIterator->value(); 
    */

    //
    // Let's build the material table
    G4MaterialTable *theMaterialTable = (G4MaterialTable*)G4Material::GetMaterialTable();
    TMXsec **secTable = new TMXsec*[nmaterials];
    //
    // now let's go for cross sections
    char string[1024]="\0";
    // Array of cross sections, say max 10 processes per particle
    G4float *pxsec = new G4float[maxproc*nbins];
    G4float *msang = new G4float[nbins];
    G4float *mslen = new G4float[nbins];
    G4float *dedx  = new G4float[nbins];
    G4int   *pdic  = new G4int[maxproc];

    // Print all processes for all particles for reference
    G4double maxwidth=0;
    G4double minlife =1e9;
    for(G4int i=0; i<np; ++i) {
       particle = theParticleTable->GetParticle(i);
       char sl[3]="ll";
       if(particle->IsShortLived()) strcpy(sl,"sl");
       G4double life = particle->GetPDGLifeTime()/s;
       G4double width = particle->GetPDGWidth();
       sprintf(string,"%-20s: %s %12.5g %12.5g %12.5g",(const char *)particle->GetParticleName(),sl,life,width,life*width/hbar_Planck*s);
       G4int nproc=0;
       G4ProcessManager* pManager = particle->GetProcessManager();
       if ( pManager == 0) {
	  G4cout << "No process manager (no X-sec!) for " << particle->GetParticleName() << G4endl;
       } else {
	  // Here we get the process list for this particle
	  G4ProcessVector* pList = pManager->GetProcessList();

	  // loop over all processes
	  for (G4int idx=0; idx<pList->size(); idx++) {
	     const G4VProcess* p = (*pList)[idx];
	     G4int ptype = p->GetProcessType();
	     if(ptype !=1 && ptype !=6) ++nproc;
	     sprintf(&string[strlen(string)]," [%s,%d,%d]",
		     (const char *)p->GetProcessName(),p->GetProcessType(),
		     p->GetProcessSubType());
	  }
       }
       if(nproc) {
	  if(life>0 && life<minlife) minlife=life;
	  if(width>maxwidth) maxwidth=width;
       }
       printf("%s\n",string);
    }
    G4double maxwidthl=hbar_Planck/(s*minlife);
    G4double minlifew =hbar_Planck/(s*maxwidth);
    printf("Reactions for life above %g(%g) and width below %g(%g) (%g)\n",
	   minlife,maxwidthl,maxwidth,minlifew,minlife*maxwidth/hbar_Planck*s);
	     
    G4ProductionCutsTable* theCoupleTable =
       G4ProductionCutsTable::GetProductionCutsTable();
    
    size_t numOfCouples = theCoupleTable->GetTableSize();

    // loop over all materials
    G4double totsize = 0;
    for(G4int imat=0; imat<nmaterials; ++imat) {
       // Just a check that we are finding out the same thing...
       const G4Material *matt = G4Material::GetMaterial(material[imat]);
       const G4Material *mat = (*theMaterialTable)[imat+1];  // skip G4_galactic
       if(matt!=mat) printf("Funny %s %s!\n",(const char*) mat->GetName(),(const char*) mat->GetName());
       secTable[imat] = new TMXsec(mat->GetZ(),mat->GetA()*mole/g,emin,emax,nbins,np);

       G4double natomscm3 = (Avogadro*mat->GetDensity()*cm3)/(mat->GetA()*mole);
       printf("Material = %s density = %g, natoms/cm3 = %g\n",
	      (const char*) mat->GetName(),mat->GetDensity()*cm3/g,natomscm3);

       // get the couple table for em processes
       const G4MaterialCutsCouple* couple = 0;
       for (size_t i=0; i<numOfCouples; i++) {
	  couple = theCoupleTable->GetMaterialCutsCouple(i);
	  if (couple->GetMaterial() == mat) break;
       }

       // loop over all particles
       for(G4int i=0; i<np; ++i) {
	  particle = theParticleTable->GetParticle(i);
	  sprintf(string,"%-20s:",(const char *)particle->GetParticleName());
	  G4ProcessManager* pManager = particle->GetProcessManager();
	  if ( pManager == 0) {
	     G4cout << "No process manager (no X-sec!) for " << particle->GetParticleName() << G4endl;
	  } else {
	     // Here we get the process list for this particle
	     G4ProcessVector* pList = pManager->GetProcessList();
	     G4int nproc=0;

	     // loop over all processes
	     for (G4int idx=0; idx<pList->size(); idx++) {
		G4VProcess* p = (*pList)[idx];
		sprintf(&string[strlen(string)]," [%s,%d,%d]",
			(const char *)p->GetProcessName(),p->GetProcessType(),
			p->GetProcessSubType());

		// if it is transportation or decay we bail out
		if(p->GetProcessType()==1||p->GetProcessType()==6) continue;

		// Let's start with what we know
		if(p->GetProcessType() == 4) {
		   // no parametrization for Z > 92 and inhelastic (but why a crash??)
		   if(mat->GetZ() > 92 && 
		      ( i == 382 || i == 383 ) && 
		      p->GetProcessSubType() == 121) continue; 
		   totsize += nbins;
		   if(p->GetProcessSubType() == 151) continue; // Capture at rest, will see later

		   // Hadronic interaction -- just store x-sec
		   printf("Mat %s Part [%3d] %s adding process %s [%d,%d]\n",
			  (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
			  (const char*) p->GetProcessName(),
			  p->GetProcessType(),p->GetProcessSubType());
		   G4double en=emin;
		   G4HadronicProcess *ph = (G4HadronicProcess*)p;
		   for(G4int j=0; j<nbins; ++j) {
		      G4DynamicParticle *dp = new G4DynamicParticle(particle,G4ThreeVector(0,0,1),en);
		      pxsec[nproc*nbins+j] =  ph->GetElementCrossSection(dp,mat->GetElement(0))/barn;
		      en*=delta;
		      delete dp;
		   }
		   pdic[nproc]=ph->GetProcessType()*1000+ph->GetProcessSubType();
		   ++nproc;
		} else if (p->GetProcessType() == 2) {
		   
		   // Em process. A bit more complicated
		   
		   totsize += nbins;
		   G4VEnergyLossProcess *pt = dynamic_cast<G4VEnergyLossProcess*>(p);
		   if(pt) {
		      printf("Mat %s Part [%3d] %s adding process %s [%d,%d]\n",
			     (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
			     (const char*) p->GetProcessName(),
			     p->GetProcessType(),p->GetProcessSubType());
		      G4double en=emin;
		      for(G4int j=0; j<nbins; ++j) {
			 pxsec[nproc*nbins+j] =  pt->CrossSectionPerVolume(en,couple)*cm/natomscm3/barn;
			 en*=delta;
		      }
		      pdic[nproc]=pt->GetProcessType()*1000+pt->GetProcessSubType();
		      ++nproc;
		   } else {
		      printf("%s: %s[%d,%d] not a GVEnergyLossProcess\n",
			     (const char*) particle->GetParticleName(),
			     (const char*)p->GetProcessName(),
			     p->GetProcessType(),p->GetProcessSubType());
		   }
		}
	     }
	  }
       }
    }
    totsize += nbins*nmaterials;
    totsize /= 0.25*1024*1024;
    printf("Tot size = %11.4gMB\n",totsize);
    fh->Write();
    fh->Close();
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
