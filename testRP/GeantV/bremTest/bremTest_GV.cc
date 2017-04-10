
/**
 * @brief  bremTest_GV: GeantV-Real-physics application to test models for e-/e+ bremsstrahlung.
 * @author M Novak
 * @date   April 2017
 *
 * A GeantV-Real-physics application to test models for e-/e+ bremsstrahlung available in the
 * given physics list: eSeltzerBergerBrems, eRelativisticBrems. The models are extracted from
 * the given physics list. Note, that the corresponding Geant4 application does not use a
 * physics-list to build-up and initialise the models but instead it creates and initialises
 * the models directly.
 * Both the energy and angular distributions of the emitted bremsstrahlung photons can be tested
 * with the application (see the available options below by running with --help option)
 *
 * Run ./bremTest_GV --help    for more details!
 * The corresponding Geant4 test application is bremTest_G4.
 */



#include <iostream>
#include <vector>
#include <cstdio>
#include <ctime>

#include <getopt.h>
#include <err.h>

#include "Material.h"
#include "Element.h"
#include "MaterialCuts.h"

// vecgeom includes
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"

#include "Region.h"
#include "PhysicsListManager.h"
#include "PhysicsList.h"
#include "PhysicsParameters.h"

// include the user defined physics list
#include "UserPhysicsList.h"

#include "ELossTableManager.h"


#include "PhysicsManagerPerParticle.h"
#include "PhysicsProcess.h"

#include "Particle.h"
#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "EMPhysicsProcess.h"
#include "EMModelManager.h"
#include "EMModel.h"


#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"

#include "LightTrack.h"
#include "PhysicsData.h"

// a simple histogram class: can be changed later
#include "Hist.h"

using geantphysics::Material;
using geantphysics::Element;
using geantphysics::MaterialCuts;  // this is just to print the table

using geantphysics::PhysicsListManager;
using geantphysics::PhysicsList;
using geantphysics::PhysicsParameters;


using geantphysics::PhysicsManagerPerParticle;
using geantphysics::PhysicsProcess;

using geantphysics::Particle;
using geantphysics::Electron;
using geantphysics::Positron;
using geantphysics::Gamma;

using geantphysics::EMPhysicsProcess;

using geantphysics::ELossTableManager;

using geantphysics::EMModelManager;
using geantphysics::EMModel;

using geantphysics::SeltzerBergerBremsModel;
using geantphysics::RelativisticBremsModel;

using geantphysics::LightTrack;
using geantphysics::PhysicsData;

using userapplication::Hist;

//
// default values of the input parameters
static std::string   particleName("e-");                  // primary particle is electron
static std::string   materialName("NIST_MAT_Pb");         // material is lead
static std::string   bremModelName("eSeltzerBergerBrems");// name of the bremsstrahlung model to test
static int           numHistBins       = 100;             // number of histogram bins between min/max values
static double        numSamples        = 1.e+7;           // number of required final state samples
static double        primaryEnergy     = 0.1;             // primary particle energy in [GeV]
static double        prodCutValue      = 0.1;             // by default in length and internal units i.e. [cm]
static bool          isProdCutInLength = true;            // is the production cut value given in length ?
static bool          isAngular         = false;           // angular or energy distribution is required ?

static struct option options[] = {
  {"particle-name     (possible particle names: e-, e+)                        - default: e-"                 , required_argument, 0, 'p'},
  {"material-name     (with a NIST_MAT_ prefix; see more in material doc.)     - default: NIST_MAT_Pb"        , required_argument, 0, 'm'},
  {"primary-energy    (in internal energy units i.e. [GeV])                    - default: 0.1"                , required_argument, 0, 'E'},
  {"number-of-samples (number of required final state samples)                 - default: 1.e+7"              , required_argument, 0, 'f'},
  {"number-of-bins    (number of bins in the histogram)                        - default: 100"                , required_argument, 0, 'n'},
  {"model-name        (eSeltzerBergerBrems, eRelativisticBrems)                - default: eSeltzerBergerBrems", required_argument, 0, 'b'},
  {"cut-vale          (secondary production threshold value for all particles) - default: 0.1"                , required_argument, 0, 'c'},
  {"cut-in-energy     (is the production cut value given in energy ? )         - default: false"              , no_argument      , 0, 'e'},
  {"isangular         (angular distribution is required ?)                     - default: false"              , no_argument      , 0, 'a'},
  {"help"                                                                                                     , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
void help();

// method to create photon energy or angular distribution histogram form final states
double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut, Particle *particle,
                          EMModel *emModel, Hist *histo, bool angular);

//===========================================================================================//
int main(int argc, char *argv[]) {
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
      bremModelName = optarg;
      break;
    case 'c':
      prodCutValue = (double)strtof(optarg, NULL);
      if (prodCutValue<=0)
        errx(1, "production cut value must be positive");
      break;
    case 'e':
      isProdCutInLength = false;
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

//============================= User defined input data =====================================//
//
// Create target material: which is supposed to be a NIST Material
Material *matDetector = Material::NISTMaterial(materialName);
//
// Set particle (Electron,Positron,Gamma)
Particle *particle = nullptr;
if (particleName=="e-") {
  particle = Electron::Definition();
} else if (particleName=="e+") {
  particle = Positron::Definition();
} else {
  std::cerr<< " ***** ERROR: unknown particle name: " << particleName << std::endl;
  help();
  return 0;
}
//
// Set particle kinetic energy
double kineticEnergy    = primaryEnergy;

//
// Set process name: the name of the bremsstrahlung process is "eBrem" in the UserPhysicsList
std::string processName = "eBrem";

// Set model name:
std::string modelName   = bremModelName;

//
// Set production cuts if needed
bool   iscutinlength    = isProdCutInLength;
double gcut             = prodCutValue;
double emcut            = prodCutValue;
double epcut            = prodCutValue;
//===========================================================================================//



//============= Initialization i.e. building up and init the physics ========================//
//
// Create a dummy vecgeom::geometry with only one volume i.e. world; create a region and set production cuts
//
// create a vecgeom::LogicalVolume
vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1.,1.,1.);
vecgeom::LogicalVolume worldl(&worldParams);
// create one region and assigne to the logical volume
vecgeom::Region *aRegion = new vecgeom::Region("ARegion",iscutinlength, gcut, emcut, epcut);
worldl.SetRegion(aRegion);
// create a dummy vecgeom::media for this volume:
// Only the vecgeom::material index (that will be 0) is used later to identify our material during the init.
worldl.SetTrackingMediumPtr(new vecgeom::Medium((matDetector->GetName()).c_str(),new vecgeom::Material(),new double()));
vecgeom::GeoManager::Instance().SetWorld(worldl.Place());
vecgeom::GeoManager::Instance().CloseGeometry();

// print the material table
//std::cerr<< Material::GetTheMaterialTable();


//
// Create all MaterialCuts
//
MaterialCuts::CreateAll();
// print all MaterialCuts
// std::cout<<MaterialCuts::GetTheMaterialCutsTable()<<std::endl;




//
// Set number of regions in the PhysicsListManager: during normal physics init., it is done
// automatically in the PhysicsProcessHandler::Initialize()
PhysicsListManager::Instance().SetNumberOfRegions(vecgeom::Region::GetNumberOfRegions());


//
// Create one user physics list.
//
//   THIS IS VERY SIMILAR To Geant4 STYLE
//   The same as above but if we have only one physics list and the active region vector is not provided then that one
//   phyics list will be used in all regions

    // this is what the user needs to do:
    // 1. create the physics list
    PhysicsList *thePhysicsList = new userapplication::UserPhysicsList("UserDefined-PhysicsList");
    // The user can get the PhysicsParameters form the physics list and change default values:
    //  - unlike in case of usual simulation, now we request to compute the CSDA range table
    thePhysicsList->GetPhysicsParameters()->SetIsComputeCSDARange(true);
    // 2. register the physics list:
    PhysicsListManager::Instance().RegisterPhysicsList(thePhysicsList);

//
// print out the PhysicsParameters obejct: we ha only one physics list so we have only one physics parameter object
 std::cout<<PhysicsParameters::GetThePhysicsParametersTable()[0];

// Build all registered physics lists: during normal physics init., it is done
// automatically in the PhysicsProcessHandler::Initialize()
PhysicsListManager::Instance().BuildPhysicsLists();
//===========================================================================================//



//======================== Getting the required information =================================//
//
//  Now get the list of EM processes assigned to the particle and use them.
//
std::cout<<std::endl<<std::endl;
std::cout<< "   ================================================================================ \n"
         << "   ================================    RESULTS    ================================= \n"
         << "   ================================================================================ \n";

// first get the MaterialCuts: we have only one
const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(aRegion->GetIndex(),matDetector->GetIndex());
std::cout<< "  "<< matCut->GetMaterial() << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   MaterialCuts: \n" << matCut;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Particle       =  " << particle->GetName() << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Kinetic energy =  " << kineticEnergy/geant::MeV << "  [MeV] " << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Physics list   =  " << thePhysicsList->GetName() << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Process name   =  " << processName << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Model name     =  " << modelName << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

//
// Get all processes, assigned to this partcile and active in the given region, that has discrete part.
PhysicsManagerPerParticle *thePhysManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
// check if the particle has any processes assined to (in the current region): if not than the thePhysManager = nullptr
size_t numPostStepCandidateProcesses = 0;
std::vector<PhysicsProcess*> thePostStepCandProcVect;
if (thePhysManager) {
  thePostStepCandProcVect       = thePhysManager->GetListPostStepCandidateProcesses();
  numPostStepCandidateProcesses = thePostStepCandProcVect.size();
}
std::cout<< "   The particle has " << numPostStepCandidateProcesses
         << " processes with discrete part assinged to." << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
// terminate if there is no any post-step candidate processes assigned to the particle in the given region.
if (!numPostStepCandidateProcesses) {
  std::cout << "   ================================================================================ "<< std::endl << std::endl;
  return 0;
}


// Try to get the required process from the list of post-step-candiates
PhysicsProcess *proc = nullptr;
for (size_t i=0; i<thePostStepCandProcVect.size(); ++i) {
  if (thePostStepCandProcVect[i]->GetName()==processName) {
    proc = thePostStepCandProcVect[i];
  }
}
// Check if required process was found and terminate otherwise
if (!proc) {
  std::cout<< "   The required process with name =  " << processName << " was not found in the physics-list: "
           << thePhysicsList->GetName()<< std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  return 0;
}
// Try to get the required EM model from the process
EMPhysicsProcess *emProc = nullptr;
if (proc->GetType()==geantphysics::ProcessType::kElectromagnetic
    || proc->GetType()==geantphysics::ProcessType::kEnergyLoss) {
  emProc = static_cast<EMPhysicsProcess*>(proc);
}
if (!emProc) {
  std::cout<< "   The required process with name =  " << processName << " is not an EM process! "   << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  return 0;
}
// If the process is an EM process then try to get the required model from the EMModelManager of the EMPhysicsProcess
EMModel *emModel = nullptr;
if (emProc) {
  EMPhysicsProcess            *emProc         = static_cast<EMPhysicsProcess*>(proc);
  EMModelManager              *emModelManager = emProc->GetModelManager();
  const std::vector<EMModel*> theEMModelVect  = emModelManager->GetModelListInRegion(matCut->GetRegionIndex());
  size_t numEMModels = theEMModelVect.size();
  for (size_t i=0; i<numEMModels; ++i) {
    if (theEMModelVect[i]->GetName()==modelName) {
      emModel = theEMModelVect[i];
    }
  }
}
// Check if the required model has been found in the process
if (!emModel) {
  std::cout<< "   The required model with name =  " << modelName << " was not found in the model list of process: "
           << processName << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  return 0;
}
// Everything is ok: the process and the model have been found in the physics list.
// So the final check: the particle kinetic energy is within the min/max usage limits of the model?
if (primaryEnergy<emModel->GetLowEnergyUsageLimit() || primaryEnergy>emModel->GetHighEnergyUsageLimit()) {
  std::cout<< "   The required primary energy =  " << primaryEnergy/geant::GeV << " [GeV] is out of the range of \n"
           << "   of the model usage limits in the given physics-list : \n"
           << "     - minimum kinetic energy:  " << emModel->GetLowEnergyUsageLimit() << " [GeV]\n"
           << "     - maximum kinetic energy:  " << emModel->GetHighEnergyUsageLimit() << " [GeV]\n"
           << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  return 0;
}

//
// Create a histogram to sample the emitted photon energy(k) distribution or the angular distribution(theta):
// if energy distribution(k) is required     : log10(k/primaryEnergy) will be the variable
// if angular distribution(theta) is required: log10(1-cos(theta)*0.5) will be the variable
double xMin = std::log10(matCut->GetProductionCutsInEnergy()[0]/primaryEnergy); // log10(gcut/primaryEnergy)
double xMax = 0.1;
std::string strDistribution = "energy-distribution";
if (isAngular) {
  xMin = -12.0;
  xMax =   0.5;
  strDistribution = "angular-distribution";
}
Hist *histo     = new Hist(xMin, xMax, numHistBins);
std::cout<< "   Sampling is running : "<<  strDistribution << "........................................  " << std::endl;
// Sampling
double timeInSec = sampleDistribution(numSamples, primaryEnergy, matCut, particle, emModel, histo, isAngular);
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Time of sampling =  " << timeInSec << " [s]" << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

// fileName
char fileName[512];
if (!isAngular) {
  sprintf(fileName,"brem_GV_energy_%s.ascii",(matCut->GetMaterial()->GetName()).c_str());
} else {
  sprintf(fileName,"brem_GV_angular_%s.ascii",(matCut->GetMaterial()->GetName()).c_str());
}
std::cout<< "   Histogram is written  into file =  " << fileName << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

FILE *f = fopen(fileName,"w");
// print out the histogram
double norm = 1./numSamples;
if (!isAngular) {
  norm = 0.25/numSamples;
}
for (int i=0; i<histo->GetNumBins(); ++i) {
 //std::cout << i << " " <<h->GetX()[i]+0.5*h->GetDelta()<<"  "<<std::setprecision(8)<<h->GetY()[i]*norm<<std::endl;
 fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
}
fclose(f);


std::cout << "   ================================================================================ "<< std::endl << std::endl;

delete histo;

return 0;
}




double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut, Particle *particle,
                          EMModel *emModel, Hist *histo, bool angular) {
  double ekin       = primaryEnergy;
  double dirx       = 0.0;   // direction
  double diry       = 0.0;
  double dirz       = 1.0;
  double gamProdCut = matCut->GetProductionCutsInEnergy()[0]; // gamma production threshold
  double cost       = 1.0; // used in angular sampling
  double z          = 1.0; // used in angular sampling
  int    gvcode     = particle->GetInternalCode();            // internal code of the primary particle


  // Set up a dummy Geant::GeantTaskData and its geantphysics::PhysicsData member: they are needed in the final state
  // sampling
  Geant::GeantTaskData *td = new Geant::GeantTaskData(1,1,1);
  PhysicsData *phd = new PhysicsData();
  td->fPhysicsData = phd;
  // Set up a the primary light track for brem.
  LightTrack primaryLT;
  // And the secondary track container LightTrack secondaryLT;
  std::vector<LightTrack> secLt;  // dummy because we fill secondaries into Geant::GeantTaskData::PhysicsData
  // init time
  clock_t  start_time = clock();
  for (long int i=0; i<numSamples; ++i) {
     // we will use members:
     //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
     //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
     //  fGVcode       <==>  fGVcode   // internal particle code
     //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
     //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
     //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
     //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
     primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
     primaryLT.SetKinE(ekin);
     primaryLT.SetGVcode(gvcode);
//       primaryLT.SetTrackIndex(0); // not important now
     primaryLT.SetDirX(dirx);
     primaryLT.SetDirY(diry);
     primaryLT.SetDirZ(dirz);
//       primaryLT.SetTotalMFP(1.0); // not important now
     //
     // clean the number of secondary tracks used (in PhysicsData)
     td->fPhysicsData->SetNumUsedSecondaries(0);
     //
     // invoke the interaction
     int numSecs = emModel->SampleSecondaries(primaryLT,secLt,td);
     // get the secondary track i.e. the gamma
     if (numSecs>0) {
       LightTrack &secondaryLT = ((td->fPhysicsData->GetListOfSecondaries())[0]);
       double egamma = secondaryLT.GetKinE();
       if (egamma<gamProdCut) {
         std::cerr<<"  ***  gammae = "<<egamma << " < gamProdCut = "<<gamProdCut <<std::endl;
         exit(-1);
       }
       if (!angular) {
         histo->Fill(std::log10(egamma/ekin),1.0);
       }
       cost = secondaryLT.GetDirZ();
       z    = 0.5*(1.0-cost);
       if (angular && z>0.) {
         double val = std::log10(z);
         if (val>-12)
           histo->Fill(val,1.0);
       }
     }
   }
   clock_t end_time = clock();
   //std::cerr<<" --- Time = "<<(end_time-start_time)/(double(CLOCKS_PER_SEC))<<std::endl;
   return (end_time-start_time)/(double(CLOCKS_PER_SEC));
}




void help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  bremTest like GeantV application for testing bremsstrahlung models for e-/e+ in the given user physics-list"
           << std::endl;
  std::cout<<"\n  Usage: TestEm0_GV [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}
