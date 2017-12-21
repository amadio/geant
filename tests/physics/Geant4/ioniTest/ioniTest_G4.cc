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
//
//      File name:     ioni.cc
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 September 2011 from test30.cc
//
//      Modifications:
//      M Novak  April 2017: special version for testing GeantV
//        ionisation models againts Geant4 models.
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include <getopt.h>
#include <err.h>

#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
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

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4TouchableHistory.hh"

#include "G4MollerBhabhaModel.hh"

#include "G4ParticleChangeForLoss.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"

//
// default values of the input parameters
static std::string   particleName("e-");                  // primary particle is electron
static std::string   materialName("G4_Pb");               // material is lead
static int           numHistBins       = 100;             // number of histogram bins between min/max values
static double        numSamples        = 1.e+7;           // number of required final state samples
static double        primaryEnergy     = 100.0;           // primary particle energy in [MeV]
static double        prodCutValue      = 1.0;             // by default in length and internal units i.e. [mm]

static struct option options[] = {
  {"particle-name     (possible particle names: e-, e+)           - default: e-"     , required_argument, 0, 'p'},
  {"material-name     (with a G4_ prefix i.e. NIST material)      - default: G4_Pb"  , required_argument, 0, 'm'},
  {"primary-energy    (in internal energy units i.e. [MeV])       - default: 100"    , required_argument, 0, 'E'},
  {"number-of-samples (number of required final state samples)    - default: 1.e+7"  , required_argument, 0, 'f'},
  {"number-of-bins    (number of bins in the histogram)           - default: 100"    , required_argument, 0, 'n'},
  {"cut-vale          (secondary production threshold [mm])       - default: 1.0"    , required_argument, 0, 'c'},
  {"help"                                                                            , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
void help();



int main(int argc, char** argv) {
  //
  // Get input parameters
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "hp:m:E:f:n:c:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'p':
       particleName = optarg;
        break;
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
    case 'c':
      prodCutValue = (double)strtof(optarg, NULL);
      if (prodCutValue<=0)
        errx(1, "production cut value must be positive");
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

  G4cout << "========================================================" << G4endl;
  G4cout << "======        Ioni Test Starts        ========" << G4endl;
  G4cout << "========================================================" << G4endl;


  G4String mname(materialName);       // material
  G4double   energy  = primaryEnergy;   // primary energy of the e-/e+
  G4double   stat    = numSamples;    // number of samples
  //G4int      verbose = 1;
  G4Material *mat    = G4NistManager::Instance()->FindOrBuildMaterial(mname);

  // Set random engine to MTwist: the same that we use in GeantV (we use the std c++11 inmp.)
  //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
  //CLHEP::HepRandom::setTheEngine(new CLHEP::Ranlux64Engine);

  //CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  //G4cout<< "  rndmEngine = " <<rndmEngine->name() <<G4endl;

  if(!mat) { exit(1); }

  // Track
  G4ThreeVector aPosition  = G4ThreeVector(0.,0.,0.);
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4Electron::Electron();
  G4Positron::Positron();
  G4ParticleDefinition* part = nullptr;
  G4String pname = "electron";
  if (particleName=="e-") {
    part  = G4Electron::Electron();
  } else if (particleName=="e+") {
    part  = G4Positron::Positron();
    pname = "positron";
  } else {
    G4cout<< "  *** unknown particle name = " << particleName << G4endl;
    help();
    return 0;
  }

  G4DynamicParticle dParticle(part,aDirection,energy);

  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  //--------- Geometry definition
  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box*           sFrame   = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame   = new G4LogicalVolume(sFrame,mat,"Box",0,0,0);
  G4PVPlacement*   pFrame   = new G4PVPlacement(0,G4ThreeVector(),"Box",lFrame,0,false,0);

  G4DataVector cuts;
  cuts.push_back(std::max(prodCutValue,990*eV));
  G4ProductionCuts* pcut = new G4ProductionCuts();
  pcut->SetProductionCut(cuts[0], 0); // set cut for e-
  pcut->SetProductionCut(cuts[0], 1); // set cut for e-
  pcut->SetProductionCut(cuts[0], 2); // set cut for e-
  pcut->SetProductionCut(cuts[0], 3); // set cut for e-

  G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(mat, pcut);
  couple->SetIndex(0);

  G4Region* reg = new G4Region("DefaultRegionForTheWorld");
  reg->AddRootLogicalVolume(lFrame);
  reg->UsedInMassGeometry(true);
  reg->SetProductionCuts(pcut);
  reg->RegisterMaterialCouplePair(mat, couple);
//  G4cout << reg->FindCouple(mat) << G4endl;

  G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  theCoupleTable->UpdateCoupleTable(pFrame);
  // get production cut in energy for gamma
  G4double elCutEnergy = (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[0];
  G4double minE = elCutEnergy;
  if (part==G4Electron::Electron()) {
    minE *= 2.;
  }
  if (energy<=minE) {
    G4cout<< " *** Primary energy = " << energy/MeV
          << " [MeV] is <= minimum energy = " << minE/MeV
          << " [MeV] so there is no secondary e- production at this energy!"
          << G4endl;
    return 0;
  }

  // -------------------------------------------------------------------
  // -------- Start run processing
  G4StateManager* g4State=G4StateManager::GetStateManager();
  if (! g4State->SetNewState(G4State_Init)) {
    G4cout << "error changing G4state"<< G4endl;
  }
  // Initilize models

  // create model
  G4MollerBhabhaModel mbm;

  G4ParticleChangeForLoss* fParticleChange = new G4ParticleChangeForLoss();
  mbm.SetParticleChange(fParticleChange, 0);

  mbm.Initialise(part, cuts);

  G4VEmModel *model  = &mbm;


  // -------- Track
  G4Track* gTrack;
  gTrack = new G4Track(&dParticle,0.0,aPosition);
  G4TouchableHandle fpTouchable(new G4TouchableHistory());
  gTrack->SetTouchableHandle(fpTouchable);
  const G4Track* track = const_cast<G4Track*>(gTrack);
  // -------- Step
  if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle)) {
    G4cout << "G4StateManager PROBLEM! " << G4endl;
  }
  // -------- Event loop
  std::vector<G4DynamicParticle*> vdp;
  vdp.reserve(1);
  dParticle.SetKineticEnergy(energy);


  // print outs
  G4cout<< mat;
  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Particle       =  " << part->GetParticleName()                                    << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Kinetic energy =  " << energy/MeV << "  [MeV] "                                   << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Model name     =  " << model->GetName()                                           << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  // compute some integrated quantities using the G4VEmModel interface methods
  // check if we compute atomic-cross section: only for single elemnt materials
  G4bool isSingleElementMaterial = false;
  if (mat->GetNumberOfElements()==1) {
    isSingleElementMaterial = true;
  }
  //
  // Note: restrictedDEDX and unrestrictedDEDX will be non-zero value only in case of EMModels that has the ComputeDEDX
  //       method implemented i.e. in case of energy loss intercations (i.e. ioni. and brem.)
  G4double restrictedDEDX          = 0.0;
  G4double unRestrictedDEDX        = 0.0;
  // Note: atomicCrossSection is computed only in case of materials that has single element
  G4double atomicCrossSection      = 0.0;
  G4double macroscopicCrossSection = 0.0;
  //
  // use the model to compute restricted stopping power
  restrictedDEDX   = model->ComputeDEDXPerVolume(mat, part, energy, elCutEnergy);
  // use the model to compute un-restricted stopping power
  unRestrictedDEDX = model->ComputeDEDXPerVolume(mat, part, energy, energy);
  // use the model to compute atomic cross section (only in case of single element materials)
  if (isSingleElementMaterial) {
    atomicCrossSection  = model->ComputeCrossSectionPerAtom(part, energy, mat->GetZ(), mat->GetA(), elCutEnergy);
  }
  // use the model to compute macroscopic cross section
  macroscopicCrossSection = model->CrossSectionPerVolume(mat, part, energy, elCutEnergy);
  //
  // print out integrated quantities:
  // -atomic cross section
  if (isSingleElementMaterial) {
    std::cout<< "   cross section per atom      :";
    std::cout<< std::setw(14) << std::scientific << std::right << atomicCrossSection/barn
             << std::setw(14) << std::left << "     [barn]";
    std::cout<<std::endl;
  }
  //
  // -macroscopic cross section
  std::cout<< "   cross section per volume    :";
  std::cout<< std::setw(14) << std::scientific << std::right << macroscopicCrossSection/(1./cm)
           << std::setw(14) << std::left << "     [1/cm]";
  std::cout<<std::endl;
  //
  // -restricted stopping power
  std::cout<< "   resricted dE/dx  (MeV/cm)   :";
  std::cout<< std::setw(14) << std::scientific << std::right << restrictedDEDX/(MeV/cm)
           << std::setw(14) << std::left << "   [MeV/cm]";
  std::cout<<std::endl;
  //
  // -unrestricted stopping power
  std::cout<< "   unresricted dE/dx (MeV/cm)  :";
  std::cout<< std::setw(14) << std::scientific << std::right << unRestrictedDEDX/(MeV/cm)
           << std::setw(14) << std::left << "   [MeV/cm]";
  std::cout<<std::endl;
  //===========================================================================================//

  // prepare histograms:
  // ------- Histograms name
  // energy distribution(t): is sampled in varibale t/primaryEnergy
  // angular distribution(theta): is sampled in variable cos(theta)
  //
  Histo    histo;
  // strip down the G4_ prefix from material name
  std::string nname = mname;
  std::string ss    = "G4_";
  std::string::size_type ipos = nname.find(ss);
  if (ipos!=std::string::npos) {
     nname.erase(ipos, ss.length());
  }
  // set histo name prefix
  G4String hname    = "ioni_G4_";
  //
  // set up a histogram for the secondary e- energy(t) : t/primaryEnergy
  G4int    nbins    = numHistBins;
  G4double xmin     = 0.;
  G4double xmax     = 1.;
  G4String theHName = hname + "secondary_energy_" + nname;
  histo.Add1D("1",theHName,nbins,xmin,xmax);
  //
  // set up histogram for the secondary gamma direction(theta) : cos(theta)
  xmin     = -1.;
  xmax     =  1.;
  theHName = hname + "secondary_angular_" + nname;;
  histo.Add1D("2",theHName,nbins,xmin,xmax);
  //
  // set up a histogram for the post interaction primary e-/e+ energy(E1) : E1/primaryEnergy
  xmin     = 0.0;
  xmax     = 1.0;
  theHName = hname + pname +"_energy_" + nname;;
  histo.Add1D("3",theHName,nbins,xmin,xmax);
  //
  // set up a histogram for the post interaction primary e-/e+ direction(theta) : cos(theta)
  xmin     = -1.0;
  xmax     =  1.0;
  theHName = hname + pname +"_angular_" + nname;;
  histo.Add1D("4",theHName,nbins,xmin,xmax);
  // book histos
  hname += nname;
  histo.SetFileName(hname);
  histo.Book();


  // start sampling
  G4cout<< "   -------------------------------------------------------------------------------- " << G4endl;
  G4cout<< "   Sampling is running : .........................................................  " << G4endl;
  // Sampling

  // SB ---------------------------
  G4double timeInSec = 0.0;
  G4Timer* timer = new G4Timer();
  timer->Start();
  for (long int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    model->SampleSecondaries(&vdp,couple,&dParticle,elCutEnergy,energy);
    // if there is any secondary gamma then get it
    if (vdp.size()>0) {
      // secondary
      G4double eSec    = vdp[0]->GetKineticEnergy()/energy;   // k/E_prim
      if (eSec>0.0) {
        histo.Fill(0,eSec,1.0);
      }
      G4double costSec = vdp[0]->GetMomentumDirection().z();

      histo.Fill(1,costSec,1.0);
      // go for the post interaction primary
      G4double ePrim    = fParticleChange->GetProposedKineticEnergy()/energy;
      if (ePrim>0.0) {
        histo.Fill(2,ePrim,1.0);
      }
      G4double costPrim = fParticleChange->GetProposedMomentumDirection().z();
      histo.Fill(3,costPrim,1.0);
      delete vdp[0];
      vdp.clear();
    }
  }
  timer->Stop();
  timeInSec = timer->GetRealElapsed();
  delete timer;


  //
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   Time of sampling =  " << timeInSec << " [s]" << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;

  // -------- Committing the transaction with the tree
  histo.ScaleH1(0, 0.25/stat);
  histo.ScaleH1(1, 1./stat);
  histo.ScaleH1(2, 0.25/stat);
  histo.ScaleH1(3, 1./stat);
  //
  histo.Save();
  G4cout<< "   Histogram is written  into file =  " << hname << G4endl;
  G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
  G4cout<< "   ================================================================================ "<< G4endl << G4endl;

  delete pFrame;
  delete lFrame;
  delete sFrame;
  partTable->DeleteAllParticles();

  return 0;
}



void help() {
  G4cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<G4endl;
  G4cout<<"  Geant4 ioniTest application for testing Geant4 e-/e+ model for ionization."
           << G4endl;
  G4cout<<"\n  Usage: ioniTest_G4 [OPTIONS] \n"<<G4endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  G4cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<G4endl;
}
