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
#include "G4LossTableManager.hh"
#include "G4VMultipleScattering.hh"
#include "G4VMscModel.hh"
#include "G4Step.hh"

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
#include <TPartIndex.h>

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
    // G4float *msang = new G4float[nbins];
    // G4float *mslen = new G4float[nbins];
    G4float *dedx  = new G4float[nbins];
    G4int   npdic = 0;
    short   *pdic  = new short[maxproc];
    short   *ndic = new short[100];
    char    (*cdic)[DICLEN+1] = new char[100][DICLEN+1];

    // Print all processes for all particles for reference
    G4double maxwidth=0;
    G4double minlife =1e9;
    G4int *pindex = new G4int[np];
    G4int *preac  = new G4int[np];
    G4int npreac = 0;
    char **parttab = new char*[np];
    for(G4int i=0; i<np; ++i) {
       particle = theParticleTable->GetParticle(i);
       char sl[3]="ll";
       if(particle->IsShortLived()) strcpy(sl,"sl");
       G4double life = particle->GetPDGLifeTime()/s;
       G4double width = particle->GetPDGWidth();
       sprintf(string,"%-20s (%10d): %s %12.5g %12.5g %12.5g",(const char *)particle->GetParticleName(),
	       particle->GetPDGEncoding(),sl,life,width,life*width/hbar_Planck*s);

       const char *pnam = (const char*) particle->GetParticleName();
       G4int pnl = strlen(pnam);
       parttab[i] = new char[pnl+1];
       strcpy(parttab[i],pnam);
       parttab[i][pnl]='\0';

       pindex[i]=particle->GetPDGEncoding();
       G4bool nproc=FALSE;
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
	     if(ptype !=1 && ptype !=6) nproc=TRUE;
	     sprintf(&string[strlen(string)]," [%s,%d,%d]",
		     (const char *)p->GetProcessName(),p->GetProcessType(),
		     p->GetProcessSubType());
	  }
       }
       if(nproc) {
	  preac[npreac++]=pindex[i];
	  if(life>0 && life<minlife) minlife=life;
	  if(width>maxwidth) maxwidth=width;
       }
       printf("%s\n",string);
    }
    G4double maxwidthl=hbar_Planck/(s*minlife);
    G4double minlifew =hbar_Planck/(s*maxwidth);
    printf("Reactions for life above %g(%g) and width below %g(%g) (%g)\n",
	   minlife,maxwidthl,maxwidth,minlifew,minlife*maxwidth/hbar_Planck*s);
    printf("Found %d particles with reactions\n",npreac);
	     
    G4ProductionCutsTable* theCoupleTable =
       G4ProductionCutsTable::GetProductionCutsTable();
    
    size_t numOfCouples = theCoupleTable->GetTableSize();

    // Order the particle table
    for(G4int i=0; i<np+1; ++i) 
       for(G4int j=i; j<np; ++j) {
	  if(pindex[i]>pindex[j]) {
	     G4int itmp = pindex[i];
	     pindex[i] = pindex[j];
	     pindex[j] = itmp;
	     char *ctmp = parttab[i];
	     parttab[i] = parttab[j];
	     parttab[j] = ctmp;
	  }
       }
    FILE *fout = fopen("file.txt","w");
    char line[120];
    line[0]='\0';
    for(G4int i=0; i<np; ++i) {
       if(strlen(line)+strlen(parttab[i])+3<120) {
	  strcat(line,"\"");
	  strcat(line,parttab[i]);
	  strcat(line,"\",");
       } else {
	  fprintf(fout,line);
	  fprintf(fout,"\n");
	  line[0]='\0';
       }
    }
    if(strlen(line)>0) {
       fprintf(fout,line);
       fprintf(fout,"\n");
    }

    line[0]='\0';
    for(G4int i=0; i<np; ++i) {
       char cnum[30];
       sprintf(cnum,"%d,",pindex[i]);
       if(strlen(line)+strlen(cnum)<120) {
	  strcat(line,cnum);
       } else {
	  fprintf(fout,line);
	  fprintf(fout,"\n");
	  line[0]='\0';
       }
    }
    if(strlen(line)>0) {
       fprintf(fout,line);
       fprintf(fout,"\n");
    }
    fclose(fout);

    TPartIndex::I()->SetPartTable(parttab,pindex,np);
    printf("code %d part %s\n",211,TPartIndex::I()->PartName(TPartIndex::I()->PartIndex(211)));

    // Push particle table into the TPartIndex

    

    // loop over all materials
    const G4ThreeVector  dirz(0,0,1);
    G4double totsize = 0;
    G4int npr=0;
    for(G4int imat=0; imat<nmaterials; ++imat) {
       //       printf("Material position %f %f %f\n",MaterialPosition[imat][0],MaterialPosition[imat][1],MaterialPosition[imat][2]);
       G4ThreeVector pos(MaterialPosition[imat][0],MaterialPosition[imat][1],MaterialPosition[imat][2]);
       // Just a check that we are finding out the same thing...
       G4Material *matt = G4Material::GetMaterial(material[imat]);
       const G4Material *mat = (*theMaterialTable)[imat+1];  // skip G4_galactic
       if(matt!=mat) printf("Funny %s %s!\n",(const char*) mat->GetName(),(const char*) mat->GetName());
       TMXsec *mxsec = secTable[imat] = new TMXsec(mat->GetZ(),
						   mat->GetA()*mole/g,emin,emax,nbins,npreac);

       G4double natomscm3 = (Avogadro*mat->GetDensity()*cm3)/(mat->GetA()*mole);
       printf("Material = %s density = %g, natoms/cm3 = %g\n",
	      (const char*) mat->GetName(),mat->GetDensity()*cm3/g,natomscm3);

       // get the couple table for em processes
       const G4MaterialCutsCouple* couple = 0;
       for (size_t i=0; i<numOfCouples; i++) {
	  couple = theCoupleTable->GetMaterialCutsCouple(i);
	  if (couple->GetMaterial() == mat) break;
       }
       if(couple == 0) { 
          G4cerr << "ERROR> Not found couple for material " << mat->GetName() << G4endl;
       }

       //       G4DynamicParticle dynParticle(particle, G4ThreeVector(1,0,0), 100*MeV);
       //G4Track aTrack( &dynParticle, 0.0, G4ThreeVector(0,0,0));
       G4Step step;
       //       aTrack.SetStep( &step );  

       G4StepPoint *preStepPt= step.GetPreStepPoint();
       preStepPt->SetMaterialCutsCouple(couple); 
       preStepPt->SetMaterial(matt); // const_cast<G4Material *>(mat) ); 

       // loop over all particles
       G4int kpreac=0; // number of particles with reaction, should always be the same
       for(G4int i=0; i<np; ++i) {
	  particle = theParticleTable->GetParticle(i);
	  sprintf(string,"%-20s:",(const char *)particle->GetParticleName());
	  G4ProcessManager* pManager = particle->GetProcessManager();
	  G4bool nproc=FALSE; // whether we have a process for this particle
	  G4bool bdedx=FALSE; // whether we have dedx for this particle
	  G4bool bmulsc=FALSE; // whether we have multiple scattering for this particle
	  G4int nprxs=0; // number of processes with x-secs
	  if ( pManager == 0) {
	     G4cout << "No process manager (no X-sec!) for " << particle->GetParticleName() << G4endl;
	  } else {
	     // Here we get the process list for this particle
	     G4ProcessVector* pList = pManager->GetProcessList();

	     // loop over all processes
	     for (G4int idx=0; idx<pList->size(); idx++) {
		G4VProcess* p = (*pList)[idx];
		sprintf(&string[strlen(string)]," [%s,%d,%d]",
			(const char *)p->GetProcessName(),p->GetProcessType(),
			p->GetProcessSubType());

		// Add information to the process dictionary

		G4int pcode = p->GetProcessType()*1000+p->GetProcessSubType();
		G4int pf = 1;
		for(G4int ipt=0; ipt<npdic; ++ipt) if(ndic[ipt]==pcode) {
		      pf=0;
		      break;}
		if(pf) {
		   printf("Adding process %d %d %s\n",p->GetProcessType(),p->GetProcessSubType(),(const char*) p->GetProcessName());
		   ndic[npdic]=pcode;
		   strncpy(cdic[npdic],(const char *)p->GetProcessName(),DICLEN);
		   cdic[npdic][DICLEN]='\0';
		   ++npdic;
		}
		
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
		   /*		   printf("Mat %s Part [%3d] %s adding process %s [%d,%d]\n",
			  (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
			  (const char*) p->GetProcessName(),
			  p->GetProcessType(),p->GetProcessSubType()); */
		   G4double en=emin;

		   G4HadronicProcess *ph = (G4HadronicProcess*)p;
		   for(G4int j=0; j<nbins; ++j) {
		      G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
		      pxsec[nprxs*nbins+j] =  ph->GetElementCrossSection(dp,mat->GetElement(0))/barn;
		      en*=delta;
		      delete dp;
		   }
		   pdic[nprxs]=ph->GetProcessType()*1000+ph->GetProcessSubType();
		   ++nprxs;
		   nproc=TRUE;
		} else if (p->GetProcessType() == 2) {
		   
		   // Em process. A bit more complicated
		   
		   totsize += nbins;
		   if(G4VEnergyLossProcess *pt = dynamic_cast<G4VEnergyLossProcess*>(p)) {
		      /*		      printf("Mat %s Part [%3d] %s adding process %s [%d,%d]\n",
			     (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
			     (const char*) p->GetProcessName(),
			     p->GetProcessType(),p->GetProcessSubType()); */
		      G4double en=emin;
		      for(G4int j=0; j<nbins; ++j) {
			 pxsec[nprxs*nbins+j] =  pt->CrossSectionPerVolume(en,couple)*cm/natomscm3/barn;
			 en*=delta;
		      }
		      pdic[nprxs]=pt->GetProcessType()*1000+pt->GetProcessSubType();
		      ++nprxs;
		      nproc=TRUE;
		   } else if(G4VEmProcess *pt = dynamic_cast<G4VEmProcess*>(p)) {
		      /* printf("Mat %s Part [%3d] %s adding process %s [%d,%d]\n",
			     (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
			     (const char*) p->GetProcessName(),
			     p->GetProcessType(),p->GetProcessSubType()); */
		      G4double en=emin;
		      for(G4int j=0; j<nbins; ++j) {
			 pxsec[nprxs*nbins+j] =  pt->CrossSectionPerVolume(en,couple)*cm/natomscm3/barn;
			 en*=delta;
		      }
		      pdic[nprxs]=pt->GetProcessType()*1000+pt->GetProcessSubType();
		      ++nprxs;
		      nproc=TRUE;
		   } else if(p->GetProcessSubType() == 10) {
		      // Multiple scattering, let's see what we can do here
		      // just a check...
		      if(!particle->GetPDGCharge()) {
			 printf("ms for a neutral particle %s\n",(const char *) particle->GetParticleName());
			 exit(1);
		      }
		      G4VMultipleScattering *pms = dynamic_cast<G4VMultipleScattering*>(p);
		      if(!pms) printf("Cannot cast to MS!!!\n");
		      G4double en=emin;
		      for(G4int j=0; j<nbins; ++j) {
			 G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
			 G4Track *track = new G4Track(dp,0.,pos);
                         track->SetStep( &step ); 
			 pms->StartTracking(track);
			 //			 printf("%p\n",cc);
			 //			 printf("Material %s\n",(const char*)track->GetMaterialCutsCouple()->GetMaterial()->GetName());
                         // G4VMscModel *msmod = dynamic_cast<G4VMscModel*>(pms->SelectModel(en,0));

                         // Full call to AlongStepGPIL
			 G4int nrep=1;
			 if(!(strcmp((const char*) particle->GetParticleName(),"e-") ||
			      strcmp((const char*) mat->GetName(),"G4_Fe")) && 
			    fabs(en-0.05)<1e-3) {
			    printf("Found e- in iron and energy %fMeV!!!\n",en*1000);
			    nrep=10;
			    }
                         const G4double previousStep= 0.0;
                         const G4double proposedStep= 10.0*mm;
                         G4double currentSafety= 9.0*mm;
  		         G4GPILSelection selection;
			 G4double aveang = 0;
			 G4double sigang = 0;
			 G4double avecor = 0;
			 G4double sigcor = 0;
			 for(G4int is=0; is<nrep; ++is) {
			    G4double stepSize= 
			       pms->AlongStepGetPhysicalInteractionLength(
									  *track,
									  previousStep,
									  proposedStep,
									  currentSafety,
									  &selection);
			    
                         // Else mimic call to AlongStepGPIL
                         // G4double geomLen= 100.0*mm; 
                         //   G4double trueLimit = 
                         // msmod->ComputeTruePathLengthLimit( *track, geomLen); // Called from AlongGPIL
			    
                         // Sample Scattering
			    G4ThreeVector  dirnew(0,0,0), displacement;
			    G4ParticleChangeForMSC* particleChng; 
			    
			    particleChng= dynamic_cast<G4ParticleChangeForMSC *>(pms->AlongStepDoIt( *track, step)); 
			    dirnew= *(particleChng->GetMomentumDirection()); 
			    G4double angle = dirnew.angle(dirz);
//			    printf("Correction %f, angle %f, en %f\n",proposedStep/stepSize,180*angle/pi,track->GetKineticEnergy());
			 }

			 en*=delta;
			 nproc=TRUE;
			 bmulsc=TRUE;

                         delete track;
                         delete dp;
		      }
		   } else {
		      printf("%s: %s[%d,%d] Cannot handle yet\n",
			     (const char*) particle->GetParticleName(),
			     (const char*)p->GetProcessName(),
			     p->GetProcessType(),p->GetProcessSubType());
		   }
		}
	     }
	  }
	  // Here we build the dedx tables
	  if(nproc && particle->GetPDGCharge()) {
	     // we only consider ionisation for particles that have other processes
	     // This is overshooting, because even with a proton at 100TeV and a lifetime
	     // of 7.4e-20 (which is the shortest lifetime for which we have reactions)
	     // we still have a c*tau*gamma of 2.5e-3 microns!!
	     G4double en=emin;
	     for(G4int j=0; j<nbins; ++j) {
		dedx[j] = G4LossTableManager::Instance()
		   ->GetDEDX(particle,en,couple);
		en*=delta;
	     }
	     bdedx=TRUE;
	  }
	  if(nproc) {
	     if(preac[npr%npreac]!=particle->GetPDGEncoding() &&
		preac[kpreac]!=particle->GetPDGEncoding()) {
		// check that the index table is consistent
		printf("==========================>>>> Error %d %d %d\n",
		       npr%npreac,preac[npr%npreac],particle->GetPDGEncoding());
	     }
	     // we have a material with x-sec for this particle. we add it
	     mxsec->AddPart(kpreac,particle->GetPDGEncoding(),nbins,nprxs,emin,emax);
	     if(nprxs) mxsec->AddPartXS(kpreac,pxsec,pdic);
	     if(bdedx) mxsec->AddPartIon(kpreac,dedx);
	     ++npr; // Total number of processes to calculate size
	     ++kpreac; // number of processes for this particle
	  } // end of "if we have processes" for this particle
       } // end of particle loop
       mxsec->Finalise();
       if(kpreac!=npreac) { // number of processes should not change with time!
	  printf("Error !!! kpreac(%d) != npreac(%d)\n",kpreac,npreac);
	  exit(1);
       }
       //       mxsec->Dump();
       mxsec->Write(material[imat]);
    } // end of material loop
    TPartIndex::I()->Write("PartIndex");
    fh->Write();
    fh->Close();
    totsize += nbins*nmaterials;
    totsize /= 0.25*1024*1024;
    printf("Particles with reactions = %d, tot size = %11.4gMB\n",npr/nmaterials, totsize);
      // Print dictionary
    for(G4int id=0; id<npdic; ++id) printf("Reac #%d code %-6d %s\n",id,ndic[id],cdic[id]);
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
