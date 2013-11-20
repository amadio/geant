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

#include "SampleDisInt.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
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
#include "G4TransportationManager.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4RegionStore.hh"

#include "G4Proton.hh"

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
#include <TEXsec.h>
#include <TEFstate.h>
#include <TRandom.h>
#include <TPartIndex.h>
#include <TClass.h>
#include <TVirtualStreamerInfo.h>
#include <TDatabasePDG.h>
#include <TMap.h>
#include <TObjString.h>
#include <TFinState.h>
#include "TPDecay.h"

#include <unistd.h>

extern char *optarg;
extern int optind;
extern int optopt;
extern int opterr;
extern int optreset;

using namespace CLHEP;

void usage()
{
  G4cout <<
  "NAME" << G4endl <<
  "     tabxsec -- generate physics tables for the G5 prototype from G4" << G4endl << G4endl <<
  "SYNOPSIS" << G4endl <<
  "     tabxsec [-Eesvx] [file]" << G4endl << G4endl <<
  "DESCRIPTION" << G4endl <<
  "     Run one of more G4 events for all the known elements and extract cross sections and secondary particles" << G4endl <<
  "     generated during the interaction, tabulating them for use in the G5 prototype. file is the G4 input" << G4endl <<
  "     macro. This code has been derived from the G4 example B1" << G4endl << G4endl <<
  "     The following options are available:" << G4endl <<
  "     -E ene  energy of the particles generated in GeV" << G4endl <<
  "     -e num  number of G4 events to generate" << G4endl <<
  "     -s num  number of samples of the final state to generate and store in fstat.root (default 0, no file generated)" << G4endl <<
  "     -v num  verbosity level (at the moment only 0, 1 and 2 are used)" << G4endl <<
  "     -x      generate the cross sections and store them in the file xsec.root" << G4endl <<
  "     -z nun  minimum Z of the material to treat (default 1)" << G4endl <<
  "     -Z num  maximum Z of the material to treat (default 100)" << G4endl <<
  "     -k ene  minimum kinetic energy in GeV (default 1e-8)" << G4endl <<
  "     -K ene  maximum kinetic energy in GeV (default 1e4)" << G4endl <<
  "     -n num  number of energy bins (default 100)" << G4endl <<
  "     -p pyl  physics list (default FTFP_BERT)" << G4endl <<
  "     -i      interactive graphics mode" << G4endl;
}

void CloseTiming();

#include <fenv.h>
#include "fp_exception_glibc_extension.h"

TParticlePDG *AddParticleToPdgDatabase(const G4String& name,
                                       G4ParticleDefinition* particleDefinition);

void DefineParticles();
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static TMap partDict;

int main(int argc,char** argv)
{
  
  feenableexcept(FE_OVERFLOW|FE_DIVBYZERO|FE_INVALID);
  
  G4int nsample=0;
  G4int ngener=0;
  G4int verbose=0;
  G4bool interact=FALSE;
  G4double evEnergy=1*GeV;
  G4bool xsecs=FALSE;
  G4int zmin = 1;
  G4int zmax = 100;
  G4float emin = 1e-6;
  G4float emax = 1e3;
  G4int nbins = 100;
  char physlist[80]="FTFP_BERT";
  G4int sl=0;
  /* getopt stuff */
  int c;
  int opterr = 0;
  /* end of getopt vars */
  
  /* getopt processing */
  while ((c = getopt (argc, argv, "s:e:v:E:ixz:Z:n:k:K:p:")) != -1)
    switch (c)
  {
    case 's':
      nsample=atoi(optarg);
      break;
    case 'e':
      ngener=atoi(optarg);
      break;
    case 'z':
      zmin=atoi(optarg);
      break;
    case 'Z':
      zmax=atoi(optarg);
      break;
    case 'n':
      nbins=atoi(optarg);
      break;
    case 'v':
      verbose=atoi(optarg);
      break;
    case 'E':
      sscanf(optarg,"%lfl",&evEnergy);
      break;
    case 'k':
      sscanf(optarg,"%fl",&emin);
      break;
    case 'K':
      sscanf(optarg,"%fl",&emax);
      break;
    case 'x':
      xsecs=TRUE;
      break;
    case 'i':
      interact=TRUE;
      break;
    case 'p':
      sl=strlen(optarg);
      strncpy(physlist,optarg,sl);
      physlist[sl]='\0';
      break;
    case '?':
      usage();
      return 1;
    default:
      abort ();
  }
  
  argc -= optind-1;
  argv += optind-1;
  /*  for (int index = 0; index < argc; index++)
   printf ("Non-option argument %s\n", argv[index]);*/
  if(nsample && ! ngener) {
    G4cout << "Sampling will fail if at least one event is not generated!"
    << " Settting ngener = 1" << G4endl;
    ngener = 1;
  }
  
  /* end of getopt stuff */
  
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
  G4VModularPhysicsList* physicsList=0;
  if(!strcmp(physlist,"QBBC")) {
    physicsList = new QBBC;
  } else if(!strcmp(physlist,"FTFP_BERT")) {
    physicsList = new FTFP_BERT;
  } else if(!strcmp(physlist,"FTFP_BERT_HP")) {
    physicsList = new FTFP_BERT_HP;
  } else {
    G4cout << "Unknown physics list " << physlist << G4endl;
    exit(1);
  }

  physicsList->SetVerboseLevel(verbose);
  runManager->SetUserInitialization(physicsList);
  
  // Primary generator action
  B1PrimaryGeneratorAction *gen = new B1PrimaryGeneratorAction();
  gen->SetVerbose(verbose);
  runManager->SetUserAction(gen);
  
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
  
  if (!interact) {
    if(argc!=1) {
      // batch mode
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    } else {
      UImanager->ApplyCommand("/control/verbose 2");
      UImanager->ApplyCommand("/control/saveHistory");
      UImanager->ApplyCommand("/run/verbose 2");
    }
    if(!nsample && !xsecs) {
      // only event generation is requested
      runManager->BeamOn( ngener );
    } else {
      // ---------------------------------------------------------------------
      //
      // Here we extract x-sections and final states
      //
      // Definitions
      //
      G4ParticleDefinition* particle;
      
      G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
      
      G4int np=theParticleTable->size();
      
      printf("Found %d particles\n",np);
      
      // These quantities could be defined via program options
      const G4int maxproc = 10;  // max of 10 proc per particle
      emin = emin*GeV;
      emax = emax*GeV;
      const G4double delta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
      const G4double lemin = TMath::Log10(emin);
      const G4double lemax = TMath::Log10(emax);
      
      printf("Log emin %f Log emax %f \n",lemin, lemax);
      
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
      // We might need to run a BeamOn here, not sure
      // for the moment we do not
      //    runManager->BeamOn( 0 );
      
      //
      // Let's build the material table
      G4MaterialTable *theMaterialTable = (G4MaterialTable*)G4Material::GetMaterialTable();
      //
      // now let's go for cross sections
      char string[1024]="\0";
      // Array of cross sections, say max 10 processes per particle
      G4float *pxsec = new G4float[maxproc*nbins];
      G4float *msang = new G4float[nbins];
      G4float *mslen = new G4float[nbins];
      G4float *msasig = new G4float[nbins];
      G4float *mslsig = new G4float[nbins];
      G4float *dedx  = new G4float[nbins];
      G4int   *pdic  = new G4int[maxproc];
      
      TDatabasePDG *ipdg = TDatabasePDG::Instance();
      DefineParticles(); // Add ions and other G4 particles
      // Print all processes for all particles for reference
      G4double maxwidth=0;
      G4double minlife =1e9;
      G4int *pPDG = new G4int[np];
      G4int *preac  = new G4int[np];
      TParticlePDG **pdpdg = new TParticlePDG*[np];
      G4int npreac = 0;
      G4int npchar = 0;
      G4int ntotproc = 0;
      
      // Loop over all the particles defined by G4 in order to add them
      // to the PDG database table
      for(G4int i=0; i<np; ++i) {
        particle = theParticleTable->GetParticle(i);
        pdpdg[i] = AddParticleToPdgDatabase(particle->GetParticleName(),particle);
        char sl[3]="ll";
        if(particle->IsShortLived()) strcpy(sl,"sl");
        G4double life = particle->GetPDGLifeTime()/s;
        G4double width = particle->GetPDGWidth();
        sprintf(string,"%-20s (%10d): %s %11.5g %11.5g %11.5g",(const char *)particle->GetParticleName(),
                particle->GetPDGEncoding(),sl,life,width,life*width/hbar_Planck*s);
        
        pPDG[i]=pdpdg[i]->PdgCode();
        if(!pPDG[i]) {
          // the only 0 pdg should be the rootino
          TObjString *conversion = (TObjString*) partDict.GetValue((const char *)particle->GetParticleName());
          const char *partnam = conversion->String().Data();
          if(strcmp("Rootino",partnam)) {
            printf("Fatal, PDG for particle name %s not found!!\n",partnam);
            exit(1);
          }
        }
        
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
            if(ptype !=1 && ptype !=6) {
              nproc=TRUE;
              ++ntotproc;
            }
            sprintf(&string[strlen(string)]," [%s,%d,%d]",
                    (const char *)p->GetProcessName(),p->GetProcessType(),
                    p->GetProcessSubType());
          }
        }
        TParticlePDG *proot = ipdg->GetParticle(particle->GetPDGEncoding());
        if(!proot)
          printf("Problem!! Particle %s has no ROOT equivalent\n",(const char*)particle->GetParticleName());
        if(nproc && strcmp("GenericIon",(const char *)particle->GetParticleName())) {
          preac[npreac++]=pPDG[i];
          if(life>0 && life<minlife) minlife=life;
          if(width>maxwidth) maxwidth=width;
          if(proot->Charge()) ++npchar;
        }
        printf("%s\n",string);
      }
      
      TPartIndex::I()->SetNPartReac(npreac);
      TPartIndex::I()->SetNPartCharge(npchar);
      TPartIndex::I()->SetEnergyGrid(emin/GeV,emax/GeV,nbins);
      
      printf("Particles with reactions %d, of which %d charged\n",npreac,npchar);
      printf("Average number of proc per particle %10.2g\n",(1.*ntotproc)/npreac);
      
      // Push particle table into the TPartIndex
      // Put all the particles with reaction and charge first
      // then all the particle with reactions
      // then all the rest
      
      for(G4int i=0; i<np-1; ++i) {
        G4bool inti = FALSE;
        G4bool chgi = (TDatabasePDG::Instance()->GetParticle(pPDG[i])->Charge()!=0);
        for(G4int jr=0; jr<npreac; ++jr)
          if(preac[jr]==pPDG[i]) {inti=TRUE; break;}
        for(G4int j=i; j<np; ++j) {
          G4bool intj = FALSE;
          G4bool chgj = (TDatabasePDG::Instance()->GetParticle(pPDG[j])->Charge()!=0);
          for(G4int jr=0; jr<npreac; ++jr)
            if(preac[jr]==pPDG[j]) {intj=TRUE; break;}
          if((!inti && intj) || (inti && intj && (!chgi && chgj))) {
            G4int itmp = pPDG[i];
            pPDG[i] = pPDG[j];
            pPDG[j] = itmp;
            G4bool btmp = inti;
            inti = intj;
            intj = btmp;
            btmp = chgi;
            chgi = chgj;
            chgj = btmp;
          }
        }
      }
      
      //----------------------------------------------------------------------------------------------------------
      //
      // test that the particles with reactions are contiguous...
      // to verify that our ordering worked
      G4bool endint=FALSE;
      G4bool endchg=FALSE;
      for(G4int i=0; i<np; ++i) {
        G4bool inter = FALSE;
        G4bool charg = FALSE;
        for(G4int j=0; j<npreac; ++j)
          if(preac[j]==pPDG[i]) {inter=TRUE; break;}
        charg = (TDatabasePDG::Instance()->GetParticle(pPDG[i])->Charge()!=0);
        if(!inter) {
          if(!endint) {
            printf("First particle without interaction #%d %s PDG %d\n",
                   i,TDatabasePDG::Instance()->GetParticle(pPDG[i])->GetName(),pPDG[i]);
            endint=TRUE;
          }
        } else if(endint) {
          printf("Fatal: particles with reactions not contiguous i=%d  npreac=%d\n",i,npreac);
        }
        if(!charg && !endint) {
          if(!endchg) {
            printf("First interacting particle without charge #%d %s PDG %d\n",
                   i,TDatabasePDG::Instance()->GetParticle(pPDG[i])->GetName(),pPDG[i]);
            endchg=TRUE;
          }
        } else if(!endint && endchg) {
          printf("Fatal: interacting particles with charge not contiguous i=%d  npchar=%d\n",i,npchar);
        }
      }
      //----------------------------------------------------------------------------------------------------------
      
      // Store the particle table into PartIndex
      TPartIndex::I()->SetPartTable(pPDG,np);
      // Check that "pions will be pions"
      printf("Code %d part %s\n",211,TPartIndex::I()->PartName(TPartIndex::I()->PartIndex(211)));
      
      // See what is the range of lifetimes for particles with reactions
      G4double maxwidthl=hbar_Planck/(s*minlife);
      G4double minlifew =hbar_Planck/(s*maxwidth);
      printf("Reactions for life above %g(%g) and width below %g(%g) (%g)\n",
             minlife,maxwidthl,maxwidth,minlifew,minlife*maxwidth/hbar_Planck*s);
      printf("Found %d particles with reactions\n",npreac);
      
      // Print out available reactions
      printf("Available reactions :");
      
      char line[120];
      line[0]='\0';
      char chpdf[20];
      for(G4int i=0; i<npreac; ++i) {
        snprintf(chpdf,19,"%d,",preac[i]);
        if(strlen(line)+strlen(chpdf)<120) {
          strcat(line,chpdf);
        } else {
          printf("%s\n", line);
          snprintf(line,19,"%s",chpdf);
        }
      }
      if(strlen(line)>0) printf("%s\n",line);
      
      // Write out particles
      // not very useful indeed, so it is disabled
      Bool_t outpart=kFALSE;
      if(outpart) {
        FILE *fout = fopen("file.txt","w");
        
        line[0]='\0';
        char cnum[30];
        for(G4int i=0; i<np; ++i) {
          snprintf(cnum,29,"%d,",pPDG[i]);
          if(strlen(line)+strlen(cnum)<120) {
            strcat(line,cnum);
          } else {
            fprintf(fout,"%s\n",line);
            snprintf(line,119,"%s",cnum);
          }
        }
        if(strlen(line)>0) fprintf(fout,"%s\n",line);
        
        fclose(fout);
      }
      
      //
      // Range converters to determine the cuts in range that
      // provide a fixed cut in energy
      G4VRangeToEnergyConverter *converter[4];
      converter[0] = new G4RToEConvForGamma();
      converter[1] = new G4RToEConvForElectron();
      converter[2] = new G4RToEConvForPositron();
      converter[3] = new G4RToEConvForProton();
      
      G4RegionStore *theRegionStore = (G4RegionStore*) G4RegionStore::GetInstance();
      
      G4int nRegions = theRegionStore->size();
      printf("We found %d regions\n",nRegions);
      
      G4double curxs=0;
      
      for(G4int imat=0; imat<nmaterials; ++imat) {
        G4Material *mat = (*theMaterialTable)[imat+1];  // skip G4_galactic
        G4Region *reg = theRegionStore->GetRegion(materialVec[imat]);
        
        /*
         // Printing the cuts just to see whether we understand them
         // not necessary any more
         G4ProductionCuts* pcuts = reg->GetProductionCuts();
         const char* parcuts[4]={"gamma","e-","e+","proton"};
         printf("Production Cuts (cm) : ");
         for(G4int ic=0; ic<4; ++ic) printf("%s=%12.2g  ",parcuts[ic],pcuts->GetProductionCut(ic)/cm);
         printf("\n");
         printf("Production Cuts (GeV): ");
         for(G4int ic=0; ic<4; ++ic)
         printf("%s=%12.2g  ",parcuts[ic],converter[ic]->Convert(pcuts->GetProductionCuts()[ic],mat)/GeV);
         printf("\n");
         */
        
        // Now set cuts to specified value via a simple bisection algorithm
        // The 1keV is arbitrary, however it is reasonable
        const G4double ethresh = 1*keV;
        for(G4int ic=0; ic<4; ++ic) {
          //	  G4double cut = reg->GetProductionCuts()->GetProductionCut(ic);
          G4double lmin = 1*nm;
          G4double lmax = 10*km;
          while(std::abs(lmin-lmax)>0.5*(lmin+lmax)*1e-10) {
            G4double lmid = 0.5*(lmin+lmax);
            if(converter[ic]->Convert(lmid,mat) > ethresh) lmax=lmid;
            else lmin = lmid;
          }
          reg->GetProductionCuts()->SetProductionCut(0.5*(lmin+lmax),ic);
        }
      }
      
      for(G4int i=0; i<4; ++i) delete converter[i];
      
      // Create our own vector of particles - since the particle table is not const in Geant4
      // New particles can be added during interaction generation and then everything is
      // messed up in the particle table.
      //
      // We did not decide yet how to deal with the generated particles that are not in the
      // initial list. For the moment we ignore them and assign the energy to the kerma
      // A better solution may be found later
      
      std::vector<G4ParticleDefinition*> particleVector;
      for(G4int i=0; i<np; ++i) {
        // particleVector[i] = theParticleTable->GetParticle(i);
        particleVector.push_back( theParticleTable->GetParticle(i) );
      }
      printf("Particle Vector has %lu contents.\n", particleVector.size() );
      if( particleVector.size() == 0 ) { printf("Cannot work without particles."); exit(1); }
      
      // Here we call BeamOn to be sure to initalise cross sections.
      // This is still not enough but it is necessary
      
      runManager->BeamOn( ngener );
      
      const G4ThreeVector  dirz(0,0,1);
      G4ThreeVector *pos = new G4ThreeVector(0,0,0);
      
      TFile *fh = 0;
      
      if(nsample) {
        // ------------------------------------------ Sample decays a la Geant4 ----------------------------------------
        
        fh = new TFile("fstate.root","recreate");
        fh->SetCompressionLevel(0);
        
        TFinState decayfs[np];
        for(G4int i=0; i<np; ++i) {
          particle = particleVector[i];
          if(particle->GetDecayTable()) {
            
            // if the particle is a generic ion we bail out
            if(!strcmp("GenericIon",(const char *)particle->GetParticleName())) continue;
            
            G4int partindex = TPartIndex::I()->PartIndex(pdpdg[i]->PdgCode());
            // printf("partindex %d pdg %d\n",partindex,pdpdg[i]->PdgCode());
            if(partindex<0) {
              printf("Error, unknown PDG %d for %s\n",particle->GetPDGEncoding(),
                     (const char *)particle->GetParticleName());
              exit(1);
            }
            sprintf(string,"%-20s:",(const char *)particle->GetParticleName());
            G4ProcessManager* pManager = particle->GetProcessManager();
            
            // Loop over all the process active for this particle
            if ( pManager == 0) {
              G4cout << "No process manager (no X-sec!) for " << particle->GetParticleName() << G4endl;
            } else {
              // Here we get the process list for this particle
              G4ProcessVector* pList = pManager->GetProcessList();
              
              // loop over all processes defined for this particle
              //
              for (G4int idx=0; idx<pList->size(); idx++) {
                G4VProcess* p = (*pList)[idx];
                sprintf(&string[strlen(string)]," [%s,%d,%d]",
                        (const char *)p->GetProcessName(),p->GetProcessType(),
                        p->GetProcessSubType());
                
                // Add information to the process dictionary
                
                G4int pcode = p->GetProcessType()*1000+p->GetProcessSubType();
                G4int pindex = TPartIndex::I()->ProcIndex(pcode);
                if(pindex<0) {
                  printf("Error: Found unknown process %d %d %s\n",
                         p->GetProcessType(),p->GetProcessSubType(),(const char*) p->GetProcessName());
                  exit(1);
                }
                
                if(p->GetProcessType()==6) {
                  // ----------------------------------- Decay ---------------------------------------
                  
                  if(verbose > 2)
                    printf("Part [%3d] %s adding process %s [%d,%d]\n",
                           i, (const char*) particle->GetParticleName(),
                           (const char*) p->GetProcessName(),
                           p->GetProcessType(),p->GetProcessSubType());
                  G4double en=0; // Let's try decay at rest
                  
                  G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
                  
                  if(nsample) {
                    G4Material *matt = (*theMaterialTable)[0];
                    printf("-------------------------------------------------  Sampling %s %s on %s @ %11.4e GeV ---------------------------------------------\n",
                           (const char *) particle->GetParticleName(),
                           (const char *) p->GetProcessName(),
                           (const char*) matt->GetName(),
                           en/GeV);
                    SampDisInt(matt, pos, dp, p, nsample, verbose, decayfs[partindex]);
                  }
                  delete dp;
                  // ----------------------------------- Decay ---------------------------------------
                }
              }
            }
          }
        }
        // ------------------------------------------ Sample decays a la Geant4 ----------------------------------------
        TPDecay *decayTable = new TPDecay(nsample,np,decayfs);
        fh->WriteObject(decayTable,"DecayTable");
      }
      
      // From here on we tabulate the cross sections and sample the interactions
      // These two functions may be split later
      
      G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
      size_t numOfCouples = theCoupleTable->GetTableSize();
      
      // loop over all materials defined
      //
      G4double totsize = 0;
      G4int npr=0;
      G4Navigator *nav = G4TransportationManager::GetTransportationManager()->
      GetNavigatorForTracking();
      TList *allElements = new TList();
      TEXsec *mxsec=0;
      TEFstate *mfstate=0;
      Int_t totfs=0;
      Int_t curfs=0;
      TFinState *vecfs=0;
      
      TFinState rcaptfs[np];
      for(G4int imat=0; imat<nmaterials; ++imat) {
        if(verbose) printf("Material position %f %f %f\n",MaterialPosition[imat][0],MaterialPosition[imat][1],MaterialPosition[imat][2]);
        pos->set(MaterialPosition[imat][0],MaterialPosition[imat][1],MaterialPosition[imat][2]);
        // It seems that this is not necessary
        //       nav->LocateGlobalPointAndUpdateTouchableHandle(*pos,dirz,fTouchableHandle,false);
        
        // Just a check that we are finding out the same thing...
        G4Material *matt = G4Material::GetMaterial(materialVec[imat]);
        const G4Material *mat = (*theMaterialTable)[imat+1];  // skip G4_galactic
        if(matt!=mat) printf("Funny %s %s!\n",(const char*) mat->GetName(),(const char*) mat->GetName());
        
        if(mat->GetZ()<zmin || mat->GetZ()>zmax) continue;
        
        G4int amat = mat->GetA()*mole/g;
        amat = 0; // set a = 0 to indicate natural element.
        G4double dens = mat->GetDensity()*cm3/g;
        
        // Create container class for the x-sections of this material
        allElements->Add(mxsec =
                         new TEXsec(mat->GetZ(),amat,dens,npreac));
        if(nsample)
          mfstate = new TEFstate(mat->GetZ(),amat,dens);
        
        G4double natomscm3 = (Avogadro*mat->GetDensity()*cm3)/(mat->GetA()*mole);
        printf("Material = %s density = %g, natoms/cm3 = %g\n",
               (const char*) mat->GetName(),mat->GetDensity()*cm3/g,natomscm3);
        
        // get the couple table for em processes
        const G4MaterialCutsCouple* couple = 0;
        const G4MaterialCutsCouple* ctmp = 0;
        
        for (size_t i=0; i<numOfCouples; i++) {
          ctmp = theCoupleTable->GetMaterialCutsCouple(i);
          if (ctmp->GetMaterial() == mat) {
            // this is just a sanity check
            if(couple) {
              printf("Found second couple for %s]n",
                     (const char*) mat->GetName());
              exit(1);
            }
            couple = ctmp;
            //	     break;
          }
        }
        
        if(couple == 0) {
          G4cerr << "ERROR> Not found couple for material " << mat->GetName() << G4endl;
          exit(1);
        }
        
        // Print production cuts for this material (is already done by G4, but just to see whether we understand them)
        // and make sure that they are what they should be (1keV)
        
        G4ProductionCuts *pcuts = couple->GetProductionCuts();
        const char* parcuts[4]={"gamma","e-","e+","proton"};
        printf("Production Cuts (cm) : ");
        for(G4int ic=0; ic<4; ++ic) printf("%s=%12.2g  ",parcuts[ic],pcuts->GetProductionCuts()[ic]/cm);
        printf("\n");
        printf("Production Cuts (GeV): ");
        Double_t cuts[4];
        for(G4int ic=0; ic<4; ++ic) {
          cuts[ic] = (*G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(ic))[couple->GetIndex()]/GeV;
          printf("%s=%12.2g  ",parcuts[ic],cuts[ic]);
        }
        printf("\n");
        mxsec->SetCuts(cuts);
        
        //       G4DynamicParticle dynParticle(particle, G4ThreeVector(1,0,0), 100*MeV);
        //G4Track aTrack( &dynParticle, 0.0, G4ThreeVector(0,0,0));
        
        // loop over all particles to store cross sections
        G4int kpreac=0; // number of particles with reaction, should always be the same
        
        for(G4int i=0; i<np; ++i) {
          particle = particleVector[i];
          
          // if the particle is a generic ion we bail out
          if(!strcmp("GenericIon",(const char *)particle->GetParticleName())) continue;
          
          G4int partindex = TPartIndex::I()->PartIndex(pdpdg[i]->PdgCode());
          // printf("partindex %d pdg %d\n",partindex,pdpdg[i]->PdgCode());
          if(partindex<0) {
            printf("Error, unknown PDG %d for %s\n",particle->GetPDGEncoding(),
                   (const char *)particle->GetParticleName());
            exit(1);
          }
          sprintf(string,"%-20s:",(const char *)particle->GetParticleName());
          G4ProcessManager* pManager = particle->GetProcessManager();
          G4bool nproc=FALSE; // whether we have a process for this particle
          G4bool bdedx=FALSE; // whether we have dedx for this particle
          G4bool bmulsc=FALSE; // whether we have multiple scattering for this particle
          G4int nprxs=0; // number of processes with x-secs
          
          if(nsample) {
            // In case we are requested to sample events, we sample all events for all
            // the channels and all the energies for a given particle, so we have to
            // allocate "enough" final states
            totfs = TPartIndex::I()->NProc()*nbins;
            vecfs = new TFinState[totfs];
          }
          
          // Loop over all the process active for this particle
          if ( pManager == 0) {
            G4cout << "No process manager (no X-sec!) for " << particle->GetParticleName() << G4endl;
          } else {
            // Here we get the process list for this particle
            G4ProcessVector* pList = pManager->GetProcessList();
            
            // loop over all processes defined for this particle
            //
            curfs = 0;
            for (G4int idx=0; idx<pList->size(); idx++) {
              G4VProcess* p = (*pList)[idx];
              sprintf(&string[strlen(string)]," [%s,%d,%d]",
                      (const char *)p->GetProcessName(),p->GetProcessType(),
                      p->GetProcessSubType());
              
              // Add information to the process dictionary
              
              G4int pcode = p->GetProcessType()*1000+p->GetProcessSubType();
              G4int pindex = TPartIndex::I()->ProcIndex(pcode);
              if(pindex<0) {
                printf("Error: Found unknown process %d %d %s\n",
                       p->GetProcessType(),p->GetProcessSubType(),(const char*) p->GetProcessName());
                exit(1);
              }
              
              // if it is transportation or decay we bail out
              if(p->GetProcessType()==1) {
                continue;
                // Decay we already did
              } else if(p->GetProcessType()==6) {
                continue;
                // From here on the processes we want to store
              } else if(p->GetProcessType() == fHadronic ) {  // 4
                // no parametrization for Z > 92 and inelastic (but why a crash??)
                if(mat->GetZ() > 92 &&
                   ( i == 382 || i == 383 ) &&
                   p->GetProcessSubType() == 121) continue;
                if(p->GetProcessSubType() == 151) {
                  // Capture at rest
                  continue;
                  if(verbose > 2)
                    printf("Part [%3d] %s adding process %s [%d,%d]\n",
                           i, (const char*) particle->GetParticleName(),
                           (const char*) p->GetProcessName(),
                           p->GetProcessType(),p->GetProcessSubType());
                  G4double en=0; // Let's try capture at rest
                  
                  G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
                  
                  if(nsample) {
                    printf("-------------------------------------------------  Sampling %s %s on %s @ %11.4e GeV ---------------------------------------------\n",
                           (const char *) particle->GetParticleName(),
                           (const char *) p->GetProcessName(),
                           (const char*) matt->GetName(),
                           en/GeV);
                    SampDisInt(matt, pos, dp, p, nsample, verbose, rcaptfs[partindex]);
                    //                     printf("vecfs[%d*%d+%d=%d].Print(): ",nbins,nprxs,j,nbins*nprxs+j); vecfs[nbins*nprxs+j].Print();
                  }
                  delete dp;
                  //                     printf("vecfs[%d*%d+%d=%d].Print(): ",nbins,nprxs,j,nbins*nprxs+j); vecfs[nbins*nprxs+j].Print();
                } else {
                  
                  // Hadronic interaction -- just store x-sec
                  if(verbose > 2)
                    printf("Mat %s Part [%3d] %s adding process %s [%d,%d]\n",
                           (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
                           (const char*) p->GetProcessName(),
                           p->GetProcessType(),p->GetProcessSubType());
                  G4double en=emin;
                  
                  G4HadronicProcess *ph = (G4HadronicProcess*)p;
                  G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
                  G4bool iszero=TRUE;
                  for(G4int j=0; j<nbins; ++j) {
                    curxs = ph->GetElementCrossSection(dp,mat->GetElement(0),mat);
                    if(curxs < 0 || curxs > 1e10) {
                      printf("%s %s on %s @ %f GeV: xs %14.7g\n",
                             (const char*) particle->GetParticleName(),
                             (const char*) p->GetProcessName(),
                             (const char*) mat->GetName(),
                             en/GeV,curxs);
                      if(curxs>0) curxs = 1e10;
                      else curxs=0;
                    }
                    pxsec[nprxs*nbins+j] = curxs/barn;
                    iszero &= (curxs<=0.);
                    if( /*particle == G4Proton::Proton() &&*/ nsample) {
                      // Here we sample proton interactions only when xsec>0 -- just for a test
                      if(curxs) {
                        printf("-------------------------------------------------  Sampling %s %s on %s @ %11.4e GeV ---------------------------------------\n",
                               (const char *) particle->GetParticleName(),
                               (const char *) ph->GetProcessName(),
                               (const char *) mat->GetName(),
                               en/GeV);
                        SampDisInt(matt, pos, dp, ph, nsample, verbose, vecfs[nbins*nprxs+j]);
                        //                     printf("vecfs[%d*%d+%d=%d].Print(): ",nbins,nprxs,j,nbins*nprxs+j); vecfs[nbins*nprxs+j].Print();
                      }
                      if(curfs != nbins*nprxs+j) {
                        printf("Eeeeeek! %d %d\n",curfs,nbins*nprxs+j);
                        exit(1);
                      }
                      ++curfs;
                    }
                    en*=delta;
                    dp->SetKineticEnergy(en);
                  }
                  delete dp;
                  if(pcode != ph->GetProcessType()*1000+ph->GetProcessSubType()) {
                    printf("Error: process code mismatch 1\n");
                    exit(1);
                  }
                  if(iszero) {
                    printf("%s for %s on %s has 0 x-sec\n",
                           (const char *) p->GetProcessName(),
                           (const char *) particle->GetParticleName(),
                           (const char *) mat->GetName());
                    
                    // Roll back the final states
                    if(nsample) curfs-=nbins;
                  } else {
                    // Add process only if there is at least one bin with nonzero x-sec
                    totsize += nbins;
                    pdic[nprxs]=pindex;
                    ++nprxs;
                    nproc=TRUE;
                  }
                }
              } else if (p->GetProcessType() == 2) {
                
                // Em process. A bit more complicated
                
                
                if(G4VEnergyLossProcess *ptEloss = dynamic_cast<G4VEnergyLossProcess*>(p)) {
                  if(verbose>2)
                    printf("Mat %s Part [%3d] %s adding ELoss process %s [%d,%d]\n",
                           (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
                           (const char*) p->GetProcessName(),
                           p->GetProcessType(),p->GetProcessSubType());
                  G4double en=emin;
                  
                  G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
                  G4bool iszero=TRUE;
                  for(G4int j=0; j<nbins; ++j) {
                    curxs = ptEloss->CrossSectionPerVolume(en,couple);
                    if(curxs < 0 || curxs > 1e10) {
                      printf("%s %s on %s @ %f GeV: xs %14.7g\n",
                             (const char*) particle->GetParticleName(),
                             (const char*) p->GetProcessName(),
                             (const char*) mat->GetName(),
                             en/GeV,curxs);
                      if(curxs>0) curxs = 1e10;
                      else curxs=0;
                    }
                    pxsec[nprxs*nbins+j] = curxs*cm/natomscm3/barn;
                    iszero &= (curxs<=0.);
                    if( /*particle == G4Electron::Electron() && */nsample) {
                      if(curxs) {
                        printf("-------------------------------------------------  Sampling %s %s on %s @ %11.4e GeV ---------------------------------------\n",
                               (const char *) particle->GetParticleName(),
                               (const char *) ptEloss->GetProcessName(),
                               (const char *) mat->GetName(),
                               en/GeV);
                        SampDisInt(matt, pos, dp, ptEloss, nsample, verbose, vecfs[nbins*nprxs+j]);
                      }
                      if(curfs != nbins*nprxs+j) {
                        printf("Eeeeeek! %d %d\n",curfs,nbins*nprxs+j);
                        exit(1);
                      }
                      ++curfs;
                    }
                    en*=delta;
                    dp->SetKineticEnergy(en);
                  }
                  delete dp;
                  if(pcode != ptEloss->GetProcessType()*1000+ptEloss->GetProcessSubType()) {
                    printf("Error: process code mismatch 1\n");
                    exit(1);
                  }
                  if(iszero) {
                    printf("%s for %s on %s has 0 x-sec\n",
                           (const char *) p->GetProcessName(),
                           (const char *) particle->GetParticleName(),
                           (const char *) mat->GetName());
                    // Roll back the final states
                    if(nsample) curfs-=nbins;
                  } else {
                    // Add process only if there is at least one bin with nonzero x-sec
                    totsize += nbins;
                    pdic[nprxs]=pindex;
                    ++nprxs;
                    nproc=TRUE;
                  }
                } else if(G4VEmProcess *ptEm = dynamic_cast<G4VEmProcess*>(p)) {
                  if(verbose > 2)
                    printf("Mat %s Part [%3d] %s adding EM process %s [%d,%d]\n",
                           (const char*) mat->GetName(), i, (const char*) particle->GetParticleName(),
                           (const char*) p->GetProcessName(),
                           p->GetProcessType(),p->GetProcessSubType());
                  
                  G4double en=emin;
                  G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
                  G4bool iszero=TRUE;
                  for(G4int j=0; j<nbins; ++j) {
                    curxs = ptEm->CrossSectionPerVolume(en,couple);
                    if(curxs < 0 || curxs > 1e10) {
                      printf("%s %s on %s @ %f GeV: xs %14.7g\n",
                             (const char*) particle->GetParticleName(),
                             (const char*) p->GetProcessName(),
                             (const char*) mat->GetName(),
                             en/GeV,curxs);
                      if(curxs>0) curxs = 1e10;
                      else curxs=0;
                    }
                    pxsec[nprxs*nbins+j] = curxs*cm/natomscm3/barn;
                    iszero &= (curxs<=0.);
                    if( /*particle == G4Positron::Positron() &&*/ nsample) {
                      if(curxs) {
                        printf("-------------------------------------------------  Sampling %s %s on %s @ %11.4e GeV ---------------------------------------\n",
                               (const char *) particle->GetParticleName(),
                               (const char *) ptEm->GetProcessName(),
                               (const char *) mat->GetName(),
                               en/GeV);
                        SampDisInt(matt, pos, dp, ptEm, nsample, verbose, vecfs[nbins*nprxs+j]);
                      }
                      if(curfs != nbins*nprxs+j) {
                        printf("Eeeeeek! %d %d\n",curfs,nbins*nprxs+j);
                        exit(1);
                      }
                      ++curfs;
                    }
                    en*=delta;
                    dp->SetKineticEnergy(en);
                  }
                  delete dp;
                  if(pcode != ptEm->GetProcessType()*1000+ptEm->GetProcessSubType()) {
                    printf("Error: process code mismatch 1\n");
                    exit(1);
                  }
                  if(iszero) {
                    printf("%s for %s on %s has 0 x-sec\n",
                           (const char *) p->GetProcessName(),
                           (const char *) particle->GetParticleName(),
                           (const char *) mat->GetName());
                    // Roll back the final states
                    if(nsample) curfs-=nbins;
                  } else {
                    // Add process only if there is at least one bin with nonzero x-sec
                    totsize += nbins;
                    pdic[nprxs]=pindex;
                    ++nprxs;
                    nproc=TRUE;
                  }
                  
                } else if(p->GetProcessSubType() == 10) {
                  // Multiple scattering, let's see what we can do here
                  
                  // First check particle and process
                  if(!particle->GetPDGCharge()) {
                    printf("ms for a neutral particle %s\n",(const char *) particle->GetParticleName());
                    exit(1);
                  }
                  G4VMultipleScattering *pms = dynamic_cast<G4VMultipleScattering*>(p);
                  if(!pms) { printf("Cannot cast to MS!!!\n"); exit(1); }
                  
                  G4double en=emin;
                  memset(msang,0,nbins*sizeof(G4float));
                  memset(mslen,0,nbins*sizeof(G4float));
                  memset(msasig,0,nbins*sizeof(G4float));
                  memset(mslsig,0,nbins*sizeof(G4float));
                  for(G4int j=0; j<nbins; ++j) {
                    
                    // Define the step size - as ftrR * range
                    G4double maxStep= 1.0*cm;
                    G4double ftrR= 0.1;
                    G4VMscModel *msModel = dynamic_cast<G4VMscModel*>(pms->SelectModel(en,0));
                    G4double range = msModel->GetRange(particle, en, couple);
                    const G4double proposedStep= std::min( ftrR*range, maxStep );
                    
                    //	printf("%p\n",cc);
                    //	printf("Material %s\n",(const char*)track->GetMaterialCutsCouple()->GetMaterial()->GetName());
                    
                    const G4int nrep=10;
                    
                    for(G4int is=0; is<nrep; ++is) {
                      const G4double previousStep= 0.0;
                      G4DynamicParticle *dp = new G4DynamicParticle(particle,dirz,en);
                      G4Track *track = new G4Track(dp,0.,*pos);
                      G4Step   step;
                      step.SetStepLength( proposedStep );
                      
                      G4StepPoint *preStepPt= step.GetPreStepPoint();
                      preStepPt->SetMaterialCutsCouple(couple);
                      preStepPt->SetMaterial(matt); // const_cast<G4Material *>(mat) );
                      G4double currentSafety= 1.1 * proposedStep;
                      preStepPt->SetSafety( currentSafety );
                      
                      track->SetStep( &step );
                      step.SetTrack( track );
                      G4StepPoint *postStepPt= step.GetPostStepPoint();
                      postStepPt->SetMomentumDirection( dirz );
                      postStepPt->SetSafety( currentSafety );
                      
                      pms->StartTracking(track);
                      // Value for return - safety also may be used
                      G4GPILSelection selection;
                      
                      // Option 1. Mimic call to AlongStepGPIL
                      G4double stepSize=
                      pms->AlongStepGetPhysicalInteractionLength(
                                                                 *track,
                                                                 previousStep,
                                                                 proposedStep,
                                                                 currentSafety,
                                                                 &selection);
                      
                      // Option 2. Get path correction from MSc model's method
                      // G4double geomLen= 100.0*mm;
                      //   G4double trueLimit =
                      // msmod->ComputeTruePathLengthLimit( *track, geomLen); // Called from AlongGPIL
                      
                      step.SetStepLength( stepSize );  // Set Geometry step size - mimics tracking
                      G4ThreeVector  dirnew(0,0,0), displacement;
                      
                      // Sample Scattering by calling the full DoIt
                      G4VParticleChange *pc= pms->AlongStepDoIt( *track, step);
                      // *************
                      G4ParticleChangeForMSC* particleChng= dynamic_cast<G4ParticleChangeForMSC *>(pc);
                      if( particleChng == 0) { G4cout << "ERROR> Incorrect type of Particle Change" << G4endl;  exit(1); }
                      
                      dirnew= *(particleChng->GetMomentumDirection());
                      G4double angle = dirnew.angle(dirz);
                      /*			    printf("Correction %f, angle %f, en %f\n",
                       proposedStep/stepSize,  180*angle/pi, track->GetKineticEnergy() );*/
                      
                      msang[j] += angle;
                      msasig[j] += angle*angle;
                      mslen[j] += proposedStep/stepSize;
                      mslsig[j] += (proposedStep/stepSize)*(proposedStep/stepSize);
                      
                      pc->UpdateStepForAlongStep(&step);
                      // step->UpdateTrack();
                      
                      // Post Step Do It
                      // G4VParticleChange *pc2= pms->PostStepDoIt( *track, step);
                      // ************
                      // pc2->UpdateStepForPostStep(&step);
                      // step->UpdateTrack();
                      
                      delete track;
                      delete dp;
                    }
                    msang[j] /= nrep;
                    msasig[j] = sqrt(TMath::Max(msasig[j]/nrep-msang[j]*msang[j],0.f));
                    mslen[j] /=nrep;
                    mslsig[j] = sqrt(TMath::Max(mslsig[j]/nrep-mslen[j]*mslen[j],0.f));
                    
                    en*=delta;
                  }
                  nproc=TRUE;
                  bmulsc=TRUE;
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
            // This is overshooting, because even with a particle at 100TeV and a lifetime
            // of 7.4e-20 (which is the shortest lifetime for which we have reactions)
            // we still have a c*tau*gamma of 2.5e-3 microns!!
            G4double en=emin;
            for(G4int j=0; j<nbins; ++j) {
              // Store dE/dx in GeV*g/cm^2
              dedx[j] = (G4LossTableManager::Instance()
                         ->GetDEDX(particle,en,couple)/GeV) / (mat->GetDensity()*cm3/g);
              en*=delta;
            }
            bdedx=TRUE;
          }
          if(nproc) {
            if(partindex >= TPartIndex::I()->NPartReac()) {
              printf("Error, particle %s PDG %d with reaction %d > %d\n",
                     (const char *)particle->GetParticleName(),
                     particle->GetPDGEncoding(),
                     partindex,TPartIndex::I()->NPartReac());
              exit(1);
            }
            // Add particle container
            mxsec->AddPart(partindex,pPDG[partindex],nprxs);
            // if there are reactions, add particle cross sections
            if(nprxs) mxsec->AddPartXS(partindex,pxsec,pdic);
            // if there is dedx, add it
            if(bdedx) mxsec->AddPartIon(partindex,dedx);
            // if there is multipls scattering, add it
            if(bmulsc) mxsec->AddPartMS(partindex,msang,msasig,mslen,mslsig);
            ++npr; // Total number of processes to calculate size
            ++kpreac; // number particle with processes
            if(nsample) { // We have been requested for sampling final states
              /*             if(partindex==32) {
               printf("Storing for protons ns %d nprxs %d\n", nsample, nprxs);
               for(G4int i=0; i<nprxs; ++i) {
               printf("Reaction %d %s\n",i,TPartIndex::I()->ProcName(pdic[i]));
               }
               for(G4int j=0; j<nprxs; ++j) {
               for(G4int i=0; i<nbins; ++i) {
               printf("%s %f ",TPartIndex::I()->ProcName(pdic[j]),TPartIndex::I()->EGrid()[i]);
               vecfs[nbins*j+i].Print();
               }
               }
               }
               */
              // Trim vecfs -- Ugly till I find a better way...
              TFinState * newvfs = new TFinState[nprxs*nbins];
              for(G4int ifs=0; ifs<nprxs*nbins; ++ifs) {
                newvfs[ifs]=vecfs[ifs];
              }
              mfstate->AddPart(partindex, pPDG[partindex], nsample, nprxs, pdic,newvfs);
              // Set to 0 the pointer to vecfs because now it is owned by the class
            }
          } // end of "if we have processes" for this particle
          if(nsample) delete [] vecfs;
        } // end of particle loop
        if(kpreac!=npreac) { // number of processes should not change with time!
          printf("Error !!! kpreac(%d) != npreac(%d)\n",kpreac,npreac);
          exit(1);
        }
        if(nsample) {
          fh->WriteObject(mfstate,TPartIndex::I()->EleSymb(mat->GetZ()));
          delete mfstate;
        }
      } // end of material loop
      
      // Write all cross sections
      if(xsecs) {
        TFile *fh = new TFile("xsec.root","recreate");
        fh->SetCompressionLevel(0);
        //allElements->Add(TPartIndex::I());
        TPartIndex::I()->Write("PartIndex");
        allElements->Write();
        fh->Write();
        fh->Close();
      }
      if(nsample) {
        fh->Write();
        fh->Close();
      }
      delete allElements;

      totsize += nbins*nmaterials;
      totsize /= 0.25*1024*1024;
      printf("Particles with reactions = %d, tot size = %11.4gMB\n",npr/nmaterials, totsize);
      // Print dictionary
      //    for(G4int id=0; id<npdic; ++id) printf("Reac #%d code %-6d %s\n",id,ndic[id],cdic[id]);

      for(G4int i=0; i<np; ++i) delete pdpdg[i];
      delete [] pdpdg;
      delete ipdg;
      
    }
    CloseTiming();
    
  } else {
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
//_____________________________________________________________________________
TParticlePDG *AddParticleToPdgDatabase(const G4String& name,
                                       G4ParticleDefinition* particleDefinition)
{
  // Add the particle definition in TDatabasePDG
  // Code by Ivana Hrivnacova
  
  // Return if particle was already added
  G4int pdgEncoding = particleDefinition->GetPDGEncoding();
  TParticlePDG *particlePDG=0;
  if(pdgEncoding)
    particlePDG
    = TDatabasePDG::Instance()->GetParticle(pdgEncoding);
  else {
    TObjString *conversion = (TObjString*) partDict.GetValue((const char *)particleDefinition->GetParticleName());
    if( conversion )
      particlePDG
	    = TDatabasePDG::Instance()->GetParticle(conversion->String().Data());
    else
      particlePDG
	    = TDatabasePDG::Instance()->GetParticle((const char*) name);
  }
  if ( particlePDG )  return particlePDG;
  
  // Get particle data
  G4String g4Name = particleDefinition->GetParticleName();
  G4int pdgQ = G4int(particleDefinition->GetPDGCharge()/eplus);
  // !! here we do not save dynamic charge but the static one
  G4String g4Type = particleDefinition->GetParticleType();
  G4String rootType = g4Type;
  if ( g4Type == "nucleus" ||  g4Type == "anti_nucleus") rootType="Ion";
  
  G4bool verbose = FALSE;
  if (verbose) {
    G4cout << "Adding particle to TDatabasePDG " << G4endl;
    G4cout << "   name:   " << g4Name << G4endl;
    G4cout << "   g4name: " << name << G4endl;
    G4cout << "   PDG:    " << pdgEncoding << G4endl;
    G4cout << "   pdgQ:   " << pdgQ << G4endl;
    G4cout << "   type:   " << rootType << G4endl;
  }
  
  // Add particle to TDatabasePDG
  return TDatabasePDG::Instance()
  ->AddParticle(name, g4Name,
                particleDefinition->GetPDGMass()/GeV,
                particleDefinition->GetPDGStable(),
                particleDefinition->GetPDGWidth()/GeV,
                pdgQ*3, rootType, pdgEncoding);
}


//
// public methods This code is copied from G4virtual MC
// written by I.Hrivnacova
//

//_____________________________________________________________________________
void DefineParticles()
{
  /// Add special particles with standard PDG = 0 to TDatabasePDG
  /// and map them to G4 particles objects.
  
  const Int_t kspe=50000000;
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  
  // optical phothon
  if ( !pdgDB->GetParticle(kspe+50) )
    pdgDB->AddParticle("Cherenkov", "Cherenkov", 0, kFALSE,
                       0,0,"Special",kspe+50);
  partDict.Add(new TObjString("opticalphoton"),new TObjString("Cherenkov"));
  
  // feedback phothon
  if ( !pdgDB->GetParticle(kspe+51) )
    pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton", 0, kFALSE,
                       0, 0, "Special", kspe+51);
  
  // charged rootino
  if ( !pdgDB->GetParticle(kspe+52) )
    pdgDB->AddParticle("ChargedRootino","ChargedRootino", 0, kFALSE,
                       0, -1, "Special", kspe+52);
  partDict.Add(new TObjString("chargedgeantino"), new TObjString("ChargedRootino"));
  
  // generic ion
  // This particle should not appear in tracking (as it is commented
  // in class G4GenericIon), but as it does, we map it anyway
  if ( !pdgDB->GetParticle(kspe+60) )
    pdgDB->AddParticle("GenericIon", "GenericIon",  0.938272, kTRUE,
                       0, 1, "Special", kspe+60);
  partDict.Add(new TObjString("GenericIon"),new TObjString("GenericIon"));
  
  
  // Light ions
  // Get PDG codes from Geant4 the rest as in TGeant3
  const Double_t kGeV=0.9314943228;
  const Double_t kHslash = 1.0545726663e-27;
  const Double_t kErgGeV = 1/1.6021773349e-3;
  const Double_t kHshGeV = kHslash*kErgGeV;
  const Double_t kYearsToSec = 3600*24*365.25;
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle;
  
  particle = particleTable->FindParticle("deuteron");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("Deuteron","Deuteron",2*kGeV+8.071e-3, kTRUE,
                       0, 3, "Ion", particle->GetPDGEncoding());
  }
  
  particle = particleTable->FindParticle("triton");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("Triton","Triton",3*kGeV+14.931e-3,kFALSE,
                       kHshGeV/(12.33*kYearsToSec), 3, "Ion",
                       particle->GetPDGEncoding());
  }
  
  particle = particleTable->FindParticle("alpha");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("Alpha","Alpha",4*kGeV+2.424e-3, kTRUE,
                       kHshGeV/(12.33*kYearsToSec), 6, "Ion",
                       particle->GetPDGEncoding());
  }
  
  particle = particleTable->FindParticle("He3");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("HE3", "HE3", 3*kGeV+14.931e-3, kFALSE,
                       0, 6, "Ion", particle->GetPDGEncoding());
  }
  
  // Light anti-ions
  // Get parameters from Geant4
  
  particle = particleTable->FindParticle("anti_deuteron");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("AntiDeuteron", "AntiDeuteron", 1.875613, kTRUE,
                       0, -3, "Ion", particle->GetPDGEncoding());
  }
  
  particle = particleTable->FindParticle("anti_triton");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("AntiTriton", "AntiTriton", 2.808921, kTRUE,
                       0, -3, "Ion", particle->GetPDGEncoding());
  }
  
  particle = particleTable->FindParticle("anti_alpha");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("AntiAlpha","AntiAlpha", 3.727379, kTRUE,
                       0, -6, "Ion", particle->GetPDGEncoding());
  }
  
  particle = particleTable->FindParticle("anti_He3");
  if ( particle && ! pdgDB->GetParticle(particle->GetPDGEncoding()) ) {
    pdgDB->AddParticle("AntiHE3", "AntiHE3", 2.808391, kTRUE,
                       0, -6, "Ion", particle->GetPDGEncoding());
  }
  partDict.Add(new TObjString("geantino"), new TObjString("Rootino"));
  
}


