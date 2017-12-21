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
//      File name:     photo.cc
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 September 2011 from test30.cc
//
//      Modifications:
//
//      M Bandieramonte April 2017:
//      special version for testing GeantV
//      photoelectric model againts Geant4 model.
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
#include "G4ParticleChangeForGamma.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
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

#include "G4PEEffectFluoModel.hh"
//#include "G4LivermorePhotoElectricModel_new.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4DipBustGenerator.hh"

#include "G4ParticleChangeForLoss.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"

//
// default values of the input parameters
static std::string   particleName("gamma");                         // primary particle is gamma
static std::string   materialName("G4_Pb");                         // material is lead
static std::string   photoelectricModelName("photoelectricLiv");    // name of the photoelectric model to test
static int           numHistBins       = 100;                       // number of histogram bins between min/max values
static double        numSamples        = 1.e+07;                    // number of required final state samples 1.e+7
static double        primaryEnergy     = 10;                        // primary particle energy in [MeV]
static double        prodCutValue      = 1.0;                       // by default in length and internal units i.e. [mm]
static bool          isAngular         = false;                     // angular or energy distribution is required ?

static struct option options[] = {
    {"particle-name     (possible particle names: gamma)            - default: gamma"               , required_argument, 0, 'p'},
    {"material-name     (with a G4_ prefix i.e. NIST material)      - default: G4_Pb"               , required_argument, 0, 'm'},
    {"primary-energy    (in internal energy units i.e. [MeV])       - default: 10"                  , required_argument, 0, 'E'},
    {"number-of-samples (number of required final state samples)    - default: 1.e+7"               , required_argument, 0, 'f'},
    {"number-of-bins    (number of bins in the histogram)           - default: 100"                 , required_argument, 0, 'n'},
    {"model-name        (photoelectricFluo, photoelectricLiv)       - default: photoelectricLiv"    , required_argument, 0, 'b'},
    {"cut-vale          (secondary production threshold [mm])       - default: 1.0"                 , required_argument, 0, 'c'},
    {"isangular         (angular distribution is required ?)        - default: false"               , no_argument      , 0, 'a'},
    {"help"                                                                                         , no_argument      , 0, 'h'},
    {0, 0, 0, 0}
};
void help();



int main(int argc, char** argv) {
    //
    // Get input parameters
    while (true) {
        int c, optidx = 0;
        c = getopt_long(argc, argv, "ehap:m:E:f:n:b:c:", options, &optidx);
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
            case 'b':
                photoelectricModelName = optarg;
                break;
            case 'c':
                prodCutValue = (double)strtof(optarg, NULL);
                if (prodCutValue<=0)
                    errx(1, "production cut value must be positive");
                break;
            case 'a':
                isAngular = true;
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
    G4cout << "======        Photoelectric Test Starts         ========" << G4endl;
    G4cout << "========================================================" << G4endl;
    
    
    G4String mname(materialName);           // material
    G4double   energy  = primaryEnergy;     // primary energy of the gamma
    G4int      stat    = numSamples;        // number of samples
    //G4int      verbose = 1;
    G4Material *mat    = G4NistManager::Instance()->FindOrBuildMaterial(mname);
    
    // Set random engine to MTwist: the same that we use in GeantV (we use the std c++11 inmp.)
    //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
    //CLHEP::HepRandom::setTheEngine(new CLHEP::Ranlux64Engine);
    
    //CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
    //G4cout<< "  rndmEngine = " <<rndmEngine->name() <<G4endl;
    
    if(!mat) { exit(1); }
    
    //we have only one model for the moment so we don't need this option
    if (!(photoelectricModelName=="photoelectricFluo" || photoelectricModelName=="photoelectricLiv")) {
      G4cout << "  *** unknown photoelectric. model name = " << photoelectricModelName << G4endl;
      help();
      return 0;
    }
    
    
    // Track
    G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
    G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
    G4Gamma::Gamma();
    G4ParticleDefinition* part = G4Gamma::Gamma();
    
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
    cuts.push_back(prodCutValue);
    G4ProductionCuts* pcut = new G4ProductionCuts();
    pcut->SetProductionCut(cuts[0], 0); // set cut for gamma
    pcut->SetProductionCut(cuts[0], 1); // set cut for gamma
    pcut->SetProductionCut(cuts[0], 2); // set cut for gamma
    pcut->SetProductionCut(cuts[0], 3); // set cut for gamma
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
    //G4double gammaCutEnergy = (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[0];
    //G4cout<< " gammaCutEnergy = " << gammaCutEnergy << " [MeV]" << G4endl;
    
    // -------------------------------------------------------------------
    // -------- Start run processing
    G4StateManager* g4State=G4StateManager::GetStateManager();
    if (! g4State->SetNewState(G4State_Init)) {
        G4cout << "error changing G4state"<< G4endl;
    }
    // Initilize models
    
    //  G4EmParameters::Instance()->Dump();
    // create models
    G4PEEffectFluoModel photoelectricFluo;
    G4LivermorePhotoElectricModel photoelectricLiv;
    
    // Set particle change object
    G4ParticleChangeForGamma* fParticleChange = new G4ParticleChangeForGamma();
    photoelectricLiv.SetParticleChange(fParticleChange, 0);
    photoelectricFluo.SetParticleChange(fParticleChange, 0);
    
    // Initilise models
    photoelectricLiv.Initialise(part, cuts);
    photoelectricLiv.SetCurrentCouple(couple);
    photoelectricFluo.Initialise(part, cuts);
    photoelectricFluo.SetCurrentCouple(couple);
    
    G4VEmModel *model  = &photoelectricFluo;
    if (photoelectricModelName=="photoelectricLiv") {
        model  = &photoelectricLiv;
    }
    
    // ------- Histograms name
    // if energy distribution(k) is required     : log10(k/primaryEnergy) will be the variable
    // if angular distribution(theta) is required: log10(1-cos(theta)*0.5) will be the variable
    //
    Histo    histo;
    std::string nname = mname;
    std::string ss    = "G4_";
    std::string::size_type ipos = nname.find(ss);
    if (ipos!=std::string::npos) {
        nname.erase(ipos, ss.length());
    }
    std::ostringstream strs;
    strs << energy;
    std::string str = strs.str();
    
    G4String hname    = model->GetName()+"_energy_" + str + "MeV_"+ nname;
    G4int    nbins    = numHistBins;
    G4double xmin     = -3.0;
    G4double xmax     =  3.0;
    G4String theHName = "Photoelectron energy";
    G4String strDistribution = "energy-distribution";

    
    if (isAngular) {  // log10(1-cos(theta)*0.5) will be the variable
        xmin     = -12.;
        xmax     = 0.5;
        theHName = "Photoelectron angular";
        hname    =  model->GetName()+"_angular_" + str + "_"+ nname;
        strDistribution = "angular-distribution";
    }
    histo.Add1D("1",theHName,nbins,xmin,xmax);
    histo.SetFileName(hname);
    histo.Book();
    G4cout << "Histograms are booked output file <" << hname << "> "<< G4endl;
    
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
    G4double cost, z, e;
    
    
    G4cout<< mat;
    G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();
    G4cout<< "   -------------------------------------------------------------------------------- "<< G4endl;
    G4cout<< "   Particle       =  " << part->GetParticleName()                                    << G4endl;
    G4cout<< "   -------------------------------------------------------------------------------- "<< G4endl;
    G4cout<< "   Kinetic energy =  " << energy/MeV << "  [MeV] "                                   << G4endl;
    G4cout<< "   -------------------------------------------------------------------------------- "<< G4endl;
    G4cout<< "   Model name     =  " << model->GetName()                                           << G4endl;
    G4cout<< "   -------------------------------------------------------------------------------- "<< G4endl;
    
    
    // compute some integrated quantities using the G4VEmModel interface methods
    // check if we compute atomic-cross section: only for single elemnt materials
    G4bool isSingleElementMaterial = false;
    if (mat->GetNumberOfElements()==1) {
        isSingleElementMaterial = true;
    }
    //

    // Note: atomicCrossSection is computed only in case of materials that has single element
    G4double atomicCrossSection      = 0.0;
    G4double macroscopicCrossSection = 0.0;
    //

    // use the model to compute atomic cross section (only in case of single element materials)
    if (isSingleElementMaterial) {
        atomicCrossSection  = model->ComputeCrossSectionPerAtom(part, energy, mat->GetZ(), mat->GetA());
    }
   
    //UNCOMMENT TO TEST CrossSectionPerVolume method
    //clock_t  start = clock();
    //for (int i = 0; i<stat; i++)
    // use the model to compute macroscopic cross section
    macroscopicCrossSection = model->CrossSectionPerVolume(mat, part, energy);
    //clock_t  end = clock();
    //std::cout<<"CrossSectionPerVolume ex-time: "<<(end-start)/(double(CLOCKS_PER_SEC))<<std::endl;
    
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
    //===========================================================================================//
    
    

    
    G4cout<< "   Sampling is running : "<<  strDistribution << "................................  "<< G4endl;
    // Sampling
    
    // SB ---------------------------
    G4double timeInSec = 0.0;
    
    G4Timer* timer = new G4Timer();
    timer->Start();
    for (G4int iter=0; iter<stat; ++iter) {
        
        fParticleChange->InitializeForPostStep(*track);
        model->SampleSecondaries(&vdp,couple,&dParticle,energy); ///gammaCutEnergy not needed
        
        if(vdp.size() > 0) {
            e = vdp[0]->GetKineticEnergy();
            if (!isAngular && e>0.0) {
                e = log10((energy - e)*1000);
                histo.Fill(0,e,1.0);
                }
            cost = vdp[0]->GetMomentumDirection().z();
            z = 0.5*(1.0 - cost);
            
            if (isAngular && z>0.0) {
                G4double val = std::log10(z);
                if (val>-12.){
                    histo.Fill(0,val,1.0);
                }
            }
            delete vdp[0];
            vdp.clear();
        }
    }
    timer->Stop();
    timeInSec = timer->GetRealElapsed();
    //G4cout << "photoelectricSG:  "  << *timer << G4endl;
    delete timer;
    
    G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
    G4cout<< "   Time of sampling =  " << timeInSec << " [s]" << G4endl;
    G4cout<< "   -------------------------------------------------------------------------------- "<<G4endl;
    
    // -------- Committing the transaction with the tree
    if (!isAngular) {
        histo.ScaleH1(0, 0.25/stat);
    } else {
        histo.ScaleH1(0, 1./stat);
        //
        
        
    }
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
    G4cout<<"  Geant4 photoelectricTest application for testing some Geant4 gamma models for photoelectric."
    << G4endl;
    G4cout<<"\n  Usage: photoelectricTest_G4 [OPTIONS] \n"<<G4endl;
    for (int i = 0; options[i].name != NULL; i++) {
        printf("\t-%c  --%s\n", options[i].val, options[i].name);
    }
    G4cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<G4endl;
}
