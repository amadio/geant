
/**
 * @brief  TestNudy0: simple GeantV-Low energy neutron test using Nudy.
 * @author Harphool Kumawat
 * @date   June 2018
 *
 * A simple test for Nudy package where low energy neutron (< 20MeV) data
 * are processed from the ufor a given particle, target material,
 * particle kinetic energy and secondary production threshold. EM processes, that has
 * discrete part and assigned to the given particle in the given physics-list, are
 * extracted from the physics-list and the integrated quantities are computed by calling
 * EMProcess/EMModel interface methods. Quantites, that are available in tables (built
 * at the physics initialisation), also extracted (i.e. interpolated) from these tables
 * by using the corresponding interface methods. Note, that atomic cross sections are
 * computed only for materials that has a single elemnt.
 *
 * Run ./TestNudy0_GV --help    for more details!
 */

#include <iostream>
#include <iomanip>
#include <vector>

#include <getopt.h>
#include <err.h>

#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/MaterialCuts.h"
#include "Geant/Isotope.h"
// vecgeom includes
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"

#include "Geant/Region.h"
#include "Geant/PhysicsListManager.h"
#include "Geant/PhysicsList.h"
#include "Geant/PhysicsParameters.h"

// include the user defined physics list
#include "UserPhysicsList.h"

#include "Geant/PhysicsManagerPerParticle.h"
#include "Geant/PhysicsProcess.h"

#include "Geant/Particle.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"
#include "Geant/Gamma.h"
#include "Geant/Proton.h"
#include "Geant/Neutron.h"
#include "Geant/PionPlus.h"
#include "Geant/PionMinus.h"
#include "Geant/PionZero.h"
#include "Geant/KaonPlus.h"
#include "Geant/KaonMinus.h"
#include "Geant/KaonZero.h"
#include "Geant/KaonShort.h"
#include "Geant/KaonLong.h"

#include "Geant/HadronicProcess.h"
#include "Geant/HadronicFinalStateModel.h"
#include "Geant/HadronicFinalStateModelStore.h"

#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyENDF.h"
#include "Geant/TNudyEndfSigma.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/ElasticScatteringProcess.h"
#include "Geant/DiffuseElasticModel.h"
#include "Geant/GlauberGribovElasticXsc.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/NeutronNudyElasticModel.h"
#include "Geant/NeutronNudyInelasticModel.h"
#include "Geant/NeutronNudyFissionModel.h"
#include "Geant/NeutronNudyCaptureModel.h"
#include "Geant/NeutronNudyElasticXsec.h"
#include "Geant/NeutronNudyInelasticXsec.h"
#include "Geant/NeutronNudyFissionXsec.h"
#include "Geant/NeutronNudyCaptureXsec.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/Typedefs.h"
#include "Geant/TaskData.h"
#include "Hist.h"

using geantphysics::Material;
using geantphysics::Element;
using geantphysics::MaterialCuts; // this is just to print the table
using geantphysics::Isotope;
using geantphysics::PhysicsListManager;
using geantphysics::PhysicsList;
using geantphysics::PhysicsParameters;

using geantphysics::PhysicsManagerPerParticle;
using geantphysics::PhysicsProcess;

using geantphysics::Particle;
using geantphysics::Electron;
using geantphysics::Positron;
using geantphysics::Gamma;
using geantphysics::Proton;
using geantphysics::Neutron;

using geantphysics::HadronicProcess;
using geantphysics::HadronicFinalStateModel;

using geantphysics::LightTrack;
using geantphysics::PhysicsData;
using userapplication::Hist;

//
// default values of the input parameters
static std::string particleName("n");           // primary particle is neutron
static std::string materialName("NIST_MAT_Pb"); // material is lead
static int Z                  = 82;             // Charge number of the isotope
static int N                  = 208;            // Total nucleons of the isotope
static int numHistBins1       = 50;            // number of histogram bins between min/max values
static int numHistBins2       = 1000;            // number of histogram bins between min/max values
static int numHistBins3       = 20;            // number of histogram bins between min/max values
static int numSamples         = 1.e+7;          // number of required final state samples
static double primaryEnergy   = 1E-3;           // primary particle energy in [GeV]
static double prodCutValue    = 1E-11;          // by default in length and internal units i.e. [cm]
static bool isProdCutInLength = true;           // is the production cut value given in length ?
double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut, Isotope *isotope,
                          Particle *primParticle, HadronicFinalStateModel *nudyModel, 
                          Hist *h1, Hist *h2, Hist *h3);

static struct option options[] = {
  {"particle-name    (possible particle: n) - default: n", required_argument, 0, 'p'},
  {"material-name    (with a NIST_MAT_ prefix; see more in material doc.)          - default: NIST_MAT_Pb",
  required_argument, 0, 'm'},
  {"number-of-samples (number of required final state samples)                 - default: 1.e+7", required_argument,
      0, 'f'},
  {"primary-energy   (in internal energy units i.e. [GeV])                         - default: 1E-3", required_argument,
        0, 'E'},
  {"sampling target Z    (proton number: z) - default: 82", required_argument, 0, 'z'},
  {"sampling target A    (Mass number (proton + neutron): A) - default: 208", required_argument, 0, 'a'},
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}};
  void help();
        
        //===========================================================================================//
int main(int argc, char *argv[])
{
  //
  // Get input parameters
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "p:m:E:f:c:z:a:", options, &optidx);
    if (c == -1) break;
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
          if (primaryEnergy <= 0) errx(1, "primary particle energy must be positive");
          break;
        case 'f':
          numSamples = (double)strtof(optarg, NULL);
          if (numSamples <= 0) errx(1, "number of final state samples must be positive");
          break;
        case 'c':
          prodCutValue = (double)strtof(optarg, NULL);
          if (prodCutValue <= 0) errx(1, "production cut value must be positive");
          break;
        case 'e':
          isProdCutInLength = false;
          break;
        case 'z':
          Z = (double)strtof(optarg, NULL);
          if (Z <= 0) errx(1, "Proton no. must be positive");
          break;
        case 'a':
          N = (double)strtof(optarg, NULL);
          if (N <= 0) errx(1, "Mass no. must be positive");
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
  if (particleName == "n") {
    particle = Neutron::Definition();
  } else {
    std::cerr << " ***** ERROR: unknown particle name: " << particleName << std::endl;
    help();
    return 0;
  }
  //
  // Set particle kinetic energy
  double kineticEnergy = primaryEnergy;
  double logKinE       = Math::Log(primaryEnergy);
  
  //============= Initialization i.e. building up and init the physics ========================//
  //
  // Create a dummy vecgeom::geometry with only one volume i.e. world; create a region and set production cuts
  //
  // create a vecgeom::LogicalVolume
  vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1., 1., 1.);
  vecgeom::LogicalVolume worldl(&worldParams);
  // create one region and assigne to the logical volume
  vecgeom::Region *aRegion = new vecgeom::Region("ARegion"); //,iscutinlength, gcut, emcut, epcut);
  worldl.SetRegion(aRegion);
  // set the material pointer in the world logical volume
  worldl.SetMaterialPtr((void *)matDetector);
  vecgeom::GeoManager::Instance().SetWorld(worldl.Place());
  vecgeom::GeoManager::Instance().CloseGeometry();
  
  // print the material table
  // std::cerr<< Material::GetTheMaterialTable();
  
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
  // create the ENDF based elastic model for elastic scattering
  PhysicsList *thePhysicsList = new userapplication::UserPhysicsList("UserDefined-PhysicsList");
  // The user can get the PhysicsParameters form the physics list and change default values:
  //  - unlike in case of usual simulation, now we request to compute the CSDA range table
  thePhysicsList->GetPhysicsParameters()->SetIsComputeCSDARange(true);
  // 2. register the physics list:
  PhysicsListManager::Instance().RegisterPhysicsList(thePhysicsList);
  
  // Build all registered physics lists: during normal physics init., it is done
  // automatically in the PhysicsProcessHandler::Initialize()
  PhysicsListManager::Instance().BuildPhysicsLists();
  //===========================================================================================//
  
  //======================== Getting the required information =================================//
  //
  //  Now get the list of EM processes assigned to the particle and use them.
  //
  std::cout << std::endl << std::endl;
  std::cout << "   ================================================================================ \n"
  << "   ================================    RESULTS    ================================= \n"
  << "   ================================================================================ \n";
  
  // first get the MaterialCuts: we have only one
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(aRegion->GetIndex(), matDetector->GetIndex());
  std::cout << "   Material       =  " << matDetector->GetName() << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   Particle       =  " << particle->GetName() << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   Kinetic energy =  " << kineticEnergy / geant::units::MeV << "  [MeV] " << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   Physics list   =  " << thePhysicsList->GetName() << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  //
  // Get all processes, assigned to this partcile and active in the given region, that has discrete part.
  PhysicsManagerPerParticle *thePhysManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // check if the particle has any processes assined to (in the current region): if not than the thePhysManager =
  // nullptr
  size_t numPostStepCandidateProcesses = 0;
  std::vector<PhysicsProcess *> thePostStepCandProcVect;
  if (thePhysManager) {
    thePostStepCandProcVect        = thePhysManager->GetListPostStepCandidateProcesses();
    numPostStepCandidateProcesses  = thePostStepCandProcVect.size();
  }
  std::cout << "   The particle has " << numPostStepCandidateProcesses << " processes with discrete part assinged to."
  << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  // terminate if there is no any post-step candidate processes assigned to the particle in the given region.
  if (!numPostStepCandidateProcesses) {
    std::cout << "   ================================================================================ " << std::endl
    << std::endl;
    return 0;
  }
  
  double compTotalMacXsec    = 0.0; // computed total macroscopic cross section i.e. summed up per-process mac. x-secs
  double getTotalMacXsec     = 0.0; // same but interpolated from the lambda table that is built at initialization
  double compTotalAtomicXsec = 0.0; // computed total atomic cross section (only if material has 1 element)
  std::vector<std::string> processNameVect;      // process names
  std::vector<double> compMacXsecPerProcessVect; // computed macroscopic cross sections per-process
  std::vector<double> getMacXsecPerProcessVect;  // same but interpolated from the lambda table that is built at init.
  std::vector<double> compAtomicXsectionVect;    // computed atomic cross sections (only if material has 1 element)
  
  //TNudyEndfRecoPoint *recopoint = new TNudyEndfRecoPoint(0,rENDF);
  
  bool isSingleElementMaterial = false;
  if (matCut->GetMaterial()->GetNumberOfElements() == 1) {
    isSingleElementMaterial = true;
  }
  
  for (size_t i = 0; i < thePostStepCandProcVect.size(); ++i) {
    PhysicsProcess *proc = thePostStepCandProcVect[i];
    processNameVect.push_back(proc->GetName());
    
    compMacXsecPerProcessVect.push_back(
      proc->ComputeMacroscopicXSection(matCut, kineticEnergy, particle, particle->GetPDGMass()));
    compTotalMacXsec += compMacXsecPerProcessVect[i];
    getMacXsecPerProcessVect.push_back(proc->GetMacroscopicXSection(matCut, kineticEnergy, logKinE, particle->GetPDGMass()));
    getTotalMacXsec += getMacXsecPerProcessVect[i];
    
    if (isSingleElementMaterial) {
      compAtomicXsectionVect.push_back(static_cast<HadronicProcess *>(proc)->GetAtomicCrossSection(
        particle->GetInternalCode(), kineticEnergy, particle->GetPDGMass(),
                                                                                                    (matCut->GetMaterial()->GetElementVector())[0], matDetector));
      
      compTotalAtomicXsec += compAtomicXsectionVect[i];
    }
  }
  
  processNameVect.push_back("total");
  compMacXsecPerProcessVect.push_back(compTotalMacXsec);
  getMacXsecPerProcessVect.push_back(getTotalMacXsec);
  compAtomicXsectionVect.push_back(compTotalAtomicXsec);
  
  std::cout << "   process name                :";
  for (size_t i = 0; i < processNameVect.size(); ++i) {
    std::cout << std::setw(16) << std::right << processNameVect[i] << std::setw(10) << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
  
  if (isSingleElementMaterial) {
    std::cout << "   cross section per atom      :";
    for (size_t i = 0; i < processNameVect.size(); ++i) {
      std::cout << std::setw(14) << std::scientific << std::right << compAtomicXsectionVect[i] / (geant::units::barn)
      << std::setw(14) << std::left << "     [barn]";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  
  std::cout << "   compCrossSectionPerVolume   :";
  for (size_t i = 0; i < processNameVect.size(); ++i) {
    std::cout << std::setw(14) << std::scientific << std::right
    << compMacXsecPerProcessVect[i] / (1. / geant::units::cm) << std::setw(14) << std::left << "     [1/cm]";
  }
  std::cout << std::endl;
  
  std::cout << "   cross section per volume    :";
  for (size_t i = 0; i < processNameVect.size(); ++i) {
    std::cout << std::setw(14) << std::scientific << std::right << getMacXsecPerProcessVect[i] / (1. / geant::units::cm)
    << std::setw(14) << std::left << "     [1/cm]";
  }
  std::cout << std::endl;
  
  double density = matCut->GetMaterial()->GetDensity() / (geant::units::g / geant::units::cm3); // density in [g/cm3]
  std::cout << "   cross section per mass      :";
  for (size_t i = 0; i < processNameVect.size(); ++i) {
    std::cout << std::setw(14) << std::scientific << std::right
    << getMacXsecPerProcessVect[i] / density / (1. / geant::units::cm) << std::setw(14) << std::left
    << "    [cm2/g]";
  }
  std::cout << std::endl;
  std::cout << std::endl;
  
  std::cout << "   mean free path (length)     :";
  for (size_t i = 0; i < processNameVect.size(); ++i) {
    double lambda = PhysicsProcess::GetAVeryLargeValue();
    if (getMacXsecPerProcessVect[i] > 0.) {
      lambda = 1. / getMacXsecPerProcessVect[i] / geant::units::cm;
    }
    std::cout << std::setw(14) << std::scientific << std::right << lambda << std::setw(14) << std::left
    << "       [cm]";
  }
  std::cout << std::endl;
  
  std::cout << "   mean free path (g/cm2)      :";
  for (size_t i = 0; i < processNameVect.size(); ++i) {
    double lambda = PhysicsProcess::GetAVeryLargeValue();
    if (getMacXsecPerProcessVect[i] > 0.) {
      lambda = getMacXsecPerProcessVect[i] / density / (1. / geant::units::cm);
      lambda = 1. / lambda;
    }
    std::cout << std::setw(14) << std::scientific << std::right << lambda << std::setw(14) << std::left
    << "    [g/cm2]";
  }
  std::cout << std::endl;
  std::cout << std::endl;
  
  std::cout << "   ================================================================================ " << std::endl
  << std::endl;
  // this portion is for the secondary parameters (angle and energy distribution)
  double xMin  = 0.0;
  double xMax  = 180.0; //
  Hist *hisAng = new Hist(xMin, xMax, numHistBins1);
  xMin         = 0.0;
  xMax         = 150; //
  Hist *hisEne = new Hist(xMin, xMax, numHistBins2);
  xMin         = 0.0;
  xMax         = 20; //
  Hist *hisSec = new Hist(xMin, xMax, numHistBins3);
  // one can test models one by one (keep only one active model and comment others)
  // geantphysics::HadronicFinalStateModel *nudyModel = new geantphysics::NeutronNudyElasticModel();
  //geantphysics::HadronicFinalStateModel *nudyModel = new geantphysics::NeutronNudyFissionModel();
  geantphysics::HadronicFinalStateModel *nudyModel = new geantphysics::NeutronNudyInelasticModel();
  
  Isotope *isotope = Isotope::GetIsotope(Z, N);
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   Model name     =  " << nudyModel->GetName() << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   Sampling is running for : "<< numSamples <<" neutrons" << std::endl;
  // call sampling method
  double timeInSec = sampleDistribution(numSamples, kineticEnergy, matCut, isotope, particle, 
                                        nudyModel, hisAng, hisEne, hisSec);
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  std::cout << "   Time of sampling =  " << timeInSec << " [s]" << std::endl;
  std::cout << "   -------------------------------------------------------------------------------- " << std::endl;
  // print out histogram to file: fileName
  char fileName[512];
  sprintf(fileName, "nudy_%s_ang", (isotope->GetName()).c_str());
  FILE *f     = fopen(fileName, "w");
  double norm = 1. / numSamples;
  Hist *histo = hisAng;
  for (int i = 0; i < histo->GetNumBins(); ++i) {
    fprintf(f, "%d\t%.8g\t%.8g\n", i, histo->GetX()[i] + 0.5 * histo->GetDelta(), histo->GetY()[i] * norm);
  }
delete histo;
fclose(f);
sprintf(fileName, "nudy_%s_ene", (isotope->GetName()).c_str());
f     = fopen(fileName, "w");
norm  = 1. / numSamples;
histo = hisEne;
for (int i = 0; i < histo->GetNumBins(); ++i) {
  fprintf(f, "%d\t%.8g\t%.8g\n", i, histo->GetX()[i] + 0.5 * histo->GetDelta(), histo->GetY()[i] * norm);
}
delete histo;
fclose(f);
sprintf(fileName, "nudy_%s_Sec", (isotope->GetName()).c_str());
f     = fopen(fileName, "w");
norm  = 1. / numSamples;
histo = hisSec;
for (int i = 0; i < histo->GetNumBins(); ++i) {
  fprintf(f, "%d\t%.8g\t%.8g\n", i, histo->GetX()[i] + 0.5 * histo->GetDelta(), histo->GetY()[i] * norm);
}
delete histo;
fclose(f);
//
std::cout << "   Histogram is written  into files ................................................" << std::endl;
std::cout << "   -------------------------------------------------------------------------------- " << std::endl;

// end
std::cout << "   ================================================================================ " << std::endl
<< std::endl;
delete nudyModel;

  return 0;
}
double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut, Isotope *isotope, 
                          Particle *particle, HadronicFinalStateModel *nudyModel, 
                          Hist *h1, Hist *h2, Hist *h3)
{
  
  double ekin = primaryEnergy;
  double dirx = 0.0; // direction
  double diry = 0.0;
  double dirz = 1.0;
  int gvcode = particle->GetInternalCode(); // internal code of the primary particle i.e. n
  double mass = particle->GetPDGMass();
  // Set up a dummy geant::TaskData and its geantphysics::PhysicsData member: they are needed in the final state
  // sampling
  geant::TaskData *td = new geant::TaskData(1, 1);
  PhysicsData *phd    = new PhysicsData();
  td->fPhysicsData    = phd;
  // Set up a the primary light track for brem.
  LightTrack primaryLT;
  // init time
  clock_t start_time = clock();
  for (long int i = 0; i < numSamples; ++i) {
    // we will use members:
    //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
    //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
    //  fGVcode       <==>  fGVcode   // internal particle code
    primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
    primaryLT.SetKinE(ekin);
    primaryLT.SetGVcode(gvcode);
    primaryLT.SetDirX(dirx);
    primaryLT.SetDirY(diry);
    primaryLT.SetDirZ(dirz);
    primaryLT.SetMass(mass);
    // clean the number of secondary tracks used (in PhysicsData)
    td->fPhysicsData->ClearSecondaries();
    //
    // invoke the interaction
    int numSecs = nudyModel->SampleFinalState(primaryLT, isotope, td);
    if (numSecs > 0) {
      LightTrack *secondaryLT = td->fPhysicsData->GetListOfSecondaries();
      for (int iSec = 0; iSec != numSecs; ++iSec) {
        h1->Fill(std::acos(secondaryLT[iSec].GetDirZ())/geant::units::deg);
        h2->Fill(secondaryLT[iSec].GetKinE()/geant::units::MeV);
      }
    } else {
      h1->Fill(std::acos(primaryLT.GetDirZ())/geant::units::deg);
      h2->Fill(primaryLT.GetKinE()/geant::units::MeV);
    }
    h3->Fill(numSecs);
  }
  clock_t end_time = clock();
  return (end_time - start_time) / (double(CLOCKS_PER_SEC));
}
void help()
{
  std::cout << "\n " << std::setw(120) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
  std::cout << "  TestNudy0 GeantV application for testing integrated physics quantities" << std::endl; 
  std::cout << "  and sampling secondary angle and energy using a given user physics-list" << std::endl;
  std::cout << "\n  Usage: TestNudy0_GV [OPTIONS] \n" << std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout << "\n " << std::setw(120) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
}
