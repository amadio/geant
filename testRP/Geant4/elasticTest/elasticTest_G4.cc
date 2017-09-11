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
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>


#include <getopt.h>
#include <err.h>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

#include "G4Timer.hh"

#include "G4TouchableHistory.hh"


#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4DecayPhysics.hh"
#include "G4ProcessManager.hh"


#include "G4HadronElastic.hh"

#include "Hist.h"
using userapplication::Hist;

//
// default values of the input parameters
static std::string   materialName("G4_C");               // material is lead
static int           numHistBins       = 100;             // number of histogram bins between min/max values
static double        numSamples        = 1.e+7;           // number of required final state samples
static double        primaryEnergy     = 1000.0;           // primary particle energy in [MeV]

static struct option options[] = {
  {"material-name     (with a G4_ prefix i.e. NIST material)      - default: G4_Pb"  , required_argument, 0, 'm'},
  {"primary-energy    (in internal energy units i.e. [MeV])       - default: 100"    , required_argument, 0, 'E'},
  {"number-of-samples (number of required final state samples)    - default: 1.e+7"  , required_argument, 0, 'f'},
  {"number-of-bins    (number of bins in the histogram)           - default: 100"    , required_argument, 0, 'n'},
  {"help"                                                                            , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
void help();



int main(int argc, char** argv) {
  //
  // Get input parameters
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "hm:E:f:n:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'm':
       materialName = optarg;
       break;
    case 'E':
      primaryEnergy = (double)strtof(optarg, NULL);
      if (primaryEnergy<=0)
        errx(1, "primary particle energy must be positive");
      break;
    case 'f':
      numSamples = (double)strtof(optarg, NULL);
      if (numSamples<=0)
        errx(1, "number of final state samples must be positive");
      break;
    case 'n':
      numHistBins = (int)strtof(optarg, NULL);
      if (numHistBins<=0)
        errx(1, "number of histogram bins must be positive");
      break;
    case 'h':
       help();
       return 0;
       break;
    default:
      help();
      errx(1, "unknown option %c", c);
    }
  }



  G4String mname(materialName);         // material
  G4double   energy  = primaryEnergy;   // primary energy of the gamma photon
  G4double   stat    = numSamples;      // number of samples
  //G4int      verbose = 1;
  G4Material *mat    = G4NistManager::Instance()->FindOrBuildMaterial(mname);


  if(!mat) { exit(1); }

  // Set random engine to MTwist: the same that we use in GeantV (we use the std c++11 inmp.)
  CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);

  //
  G4Proton::Proton();
  G4cout.setf( std::ios::scientific, std::ios::floatfield );
  
  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));
  G4DecayPhysics decays;
  decays.ConstructParticle();
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ions = partTable->GetIonTable();
  partTable->SetReadiness();
  ions->CreateAllIon();
  ions->CreateAllIsomer();
  

  //--------- Geometry definition
  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,mat,"Box",0,0,0);
  G4PVPlacement*   pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",lFrame,0,false,0);

  // -------------------------------------------------------------------
  // -------- Start run processing
  G4StateManager* g4State=G4StateManager::GetStateManager();
  if (! g4State->SetNewState(G4State_Init)) {
    G4cout << "error changing G4state"<< G4endl;;
  }

  //
  
  // Instantiate models
  G4HadronElastic elastic;

  // -------- Track
  G4ThreeVector aPosition  = G4ThreeVector(0.,0.,0.);
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  const G4ParticleDefinition* part = G4Proton::Proton();
  G4DynamicParticle dParticle(part,aDirection,energy);

  G4Track* track = new G4Track(&dParticle,0.0,aPosition);
  //
  G4TouchableHandle fpTouchable(new G4TouchableHistory());
  track->SetTouchableHandle(fpTouchable);
  
  // -------- Step
  if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle)) {
    G4cout << "G4StateManager PROBLEM! " << G4endl;
  }
  
  // set the EmModel
  G4HadronicInteraction *model = &elastic;

  // print outs
  G4cout<< mat;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Particle       =  " << part->GetParticleName()                                    << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Kinetic energy =  " << energy/MeV << "  [MeV] "                                   << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;


  double xMin =  0.0;
  double xMax =  1.0;
  Hist *histo_t    = new Hist(xMin, xMax, numHistBins);
  
  // ------- Histograms name
  G4String nname = mname; // material name
  G4String ss    = "G4_"; // strip the G4_ prefix
  G4String::size_type ipos = nname.find(ss);
  if (ipos!=std::string::npos) {
     nname.erase(ipos, ss.length());
  }
  //

  // start sampling
  G4cout<< "   -------------------------------------------------------------------------------- " << G4endl;
  G4cout<< "   Sampling is running : .........................................................  " << G4endl;
  // Sampling
  
  std::vector<G4DynamicParticle*> vdp;
  dParticle.SetKineticEnergy(energy);
  dParticle.SetMomentumDirection(0.,0.,1.);
  G4Timer *timer;

  // hadron projectile
  G4HadProjectile hadpro(dParticle);
  G4Nucleus targetNucleus(mat);
  G4HadFinalState* finalstate;
  
  timer = new G4Timer();
  timer->Start();
  for (long int iter=0; iter<stat; ++iter) {
    
    //    fParticleChange->InitializeForPostStep(*track);   
    //    model->SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    
    finalstate =   model->ApplyYourself(hadpro, targetNucleus);
    
    // sample momentum transfer using Lab. momentum
    //  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
    //				    G4double plab,
    //				    G4int Z, G4int A)    

    histo_t->Fill(std::acos(finalstate->GetMomentumChange().z()));

  }
  timer->Stop();
  double timeInSec = timer->GetRealElapsed();
  delete timer;

  FILE *f     = fopen("test.dat","w");
  double norm = 1./numSamples;

  for (int i=0; i<histo_t->GetNumBins(); ++i) {
    fprintf(f,"%d\t%.8g\t%.8g\n",i,histo_t->GetX()[i]+0.5*histo_t->GetDelta(),histo_t->GetY()[i]*norm);
  }
  
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Time of sampling =  " << timeInSec << " [s]" << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;

  
  delete pFrame;
  delete lFrame;
  delete sFrame;
  partTable->DeleteAllParticles();

  G4cout << "###### End of test #####" << G4endl;
}


void help() {
  G4cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<G4endl;
  G4cout<<"  Geant4 conversionTest application for testing Geant4 Bethe-Heitler model for e-/e+"
        <<"  pair production by photons."
        << G4endl;
  G4cout<<"\n  Usage: conversionTest_G4 [OPTIONS] \n"<<G4endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  G4cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<G4endl;
}
