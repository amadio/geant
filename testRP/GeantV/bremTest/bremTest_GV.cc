/**
 * @brief  bremTest_GV: GeantV-Real-physics test for testing e-/e+ models for bremsstrahlung photon emission.
 * @author M Novak
 * @date   April 2017
 *
 * Run ./bremTest_GV --help for more details!
 * The corresponding quantities/distributions can be obtained by using the bremTest_G4 Geant4 test.
 */


#include <iostream>
#include <iomanip>
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
#include "PhysicsParameters.h"

// just to clear them
#include "ELossTableManager.h"
#include "ELossTableRegister.h"

#include "Particle.h"
#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "EMModel.h"
#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"

#include "LightTrack.h"
#include "PhysicsData.h"

// from geantV
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"


// a simple histogram class: can be changed later
#include "Hist.h"

using geantphysics::Material;
using geantphysics::Element;
using geantphysics::MaterialCuts;

using geantphysics::PhysicsParameters;

using geantphysics::Particle;
using geantphysics::Electron;
using geantphysics::Positron;
using geantphysics::Gamma;

using geantphysics::ELossTableManager;
using geantphysics::ELossTableRegister;

// the two brem. model
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
static std::string   bremModelName("bremSB");             // name of the bremsstrahlung model to test
static bool          isUseRejection    = false;           // use rejection sampling instead of sampling tables
static int           numHistBins       = 100;             // number of histogram bins between min/max values
static double        numSamples        = 1.e+7;           // number of required final state samples
static double        primaryEnergy     = 0.1;             // primary particle energy in [GeV]
static double        prodCutValue      = 0.1;             // by default in length and internal units i.e. [cm]
static bool          isProdCutInLength = true;            // is the production cut value given in length ?

static struct option options[] = {
  {"particle-name     (possible particle names: e-, e+)                        - default: e-"                 , required_argument, 0, 'p'},
  {"material-name     (with a NIST_MAT_ prefix; see more in material doc.)     - default: NIST_MAT_Pb"        , required_argument, 0, 'm'},
  {"primary-energy    (in internal energy units i.e. [GeV])                    - default: 0.1"                , required_argument, 0, 'E'},
  {"sampling-type     (flag to switch to rejection sampling)                   - default: no"                 , no_argument      , 0, 'r'},
  {"number-of-samples (number of required final state samples)                 - default: 1.e+7"              , required_argument, 0, 'f'},
  {"number-of-bins    (number of bins in the histogram)                        - default: 100"                , required_argument, 0, 'n'},
  {"model-name        (bremSB, bremRel)                                        - default: bremSB"             , required_argument, 0, 'b'},
  {"cut-vale          (secondary production threshold value for all particles) - default: 0.1"                , required_argument, 0, 'c'},
  {"cut-in-energy     (is the production cut value given in energy ? )         - default: false"              , no_argument      , 0, 'e'},
  {"help"                                                                                                     , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
void help();

//***********************************************************************************************//
//***** THIS WILL BE MODEL SPECIFIC: contains the final state sampling and hist. building  ******//
// method to create photon energy distribution using a SeltzerBergerBremsModel as input argument
double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut,
                          Particle *primParticle, EMModel *model, Hist *histo1, Hist *histo2,
                          Hist *histo3, Hist *histo4);
//***********************************************************************************************//


//===========================================================================================//
//
int main(int argc, char *argv[]) {
  //
  //============================== Get input parameters =====================================//
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "erh:m:E:f:n:c:p:b:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'm':
       materialName = optarg;
       break;
    case 'r':
      isUseRejection = true;
      break;
    case 'E':
      primaryEnergy = strtod(optarg, NULL);
      if (primaryEnergy<=0)
        errx(1, "primary particle energy must be positive");
      break;
    case 'f':
      numSamples = strtod(optarg, NULL);
      if (numSamples<=0)
        errx(1, "number of final state samples must be positive");
      break;
    case 'n':
      numHistBins = (int)strtof(optarg, NULL);
      if (numHistBins<=0)
        errx(1, "number of histogram bins must be positive");
      break;
    case 'c':
      prodCutValue = strtod(optarg, NULL);
      if (prodCutValue<=0)
        errx(1, "production cut value must be positive");
      break;
    case 'p':
       particleName = optarg;
        break;
    case 'b':
      bremModelName = optarg;
      break;
    case 'e':
      isProdCutInLength = false;
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
  //===========================================================================================//

  //============================= Set user defined input data =================================//
  // Create target material: which is supposed to be a NIST Material
  Material *matDetector = Material::NISTMaterial(materialName);
  //
  // Set particle kinetic energy
  double kineticEnergy    = primaryEnergy;
  //
  // Set production cuts if needed
  bool   iscutinlength    = isProdCutInLength;
  double gcut             = prodCutValue;
  double emcut            = prodCutValue;
  double epcut            = prodCutValue;
  //===========================================================================================//

  if (!(bremModelName=="bremSB" || bremModelName=="bremRel")) {
    std::cout << "  *** unknown brem. model name = " << bremModelName << std::endl;
    help();
    return 0;
  }

  // Create primnary particle
  Particle    *particle   = nullptr;
  bool         isElectron = true;
  std::string  pname;
  if (particleName=="e-") {
    particle  = Electron::Definition();
    pname     = "electron";
  } else if (particleName=="e+") {
    particle  = Positron::Definition();
    pname     = "positron";
  } else {
    std::cout<< "  *** unknown particle name = " << particleName << std::endl;
    help();
    return 0;
  }


  //============= Initialization i.e. building up and init the physics ========================//
  // Create a dummy vecgeom::geometry:
  //  - with only one volume i.e. world
  //  - create a region and set production cuts
  //
  // create a vecgeom::LogicalVolume
  vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1.,1.,1.);
  vecgeom::LogicalVolume worldl(&worldParams);
  // create one region and assigne to the logical volume
  vecgeom::Region *aRegion = new vecgeom::Region("ARegion",iscutinlength, gcut, emcut, epcut);
  worldl.SetRegion(aRegion);
  // set the material pointer in the world logical volume
  worldl.SetMaterialPtr((void*)matDetector);
  vecgeom::GeoManager::Instance().SetWorld(worldl.Place());
  vecgeom::GeoManager::Instance().CloseGeometry();
  // Create all(we have only one) MaterialCuts
  MaterialCuts::CreateAll();
  //===========================================================================================//

  // if primary particle energy < gamma production cut => there is no secondary gamma production
  // So get the MaterialCuts of the target: we have only one
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(aRegion->GetIndex(),matDetector->GetIndex());
  // and get the gamma production cut energy
  double gammaCutEnergy = matCut->GetProductionCutsInEnergy()[0];
  if (kineticEnergy<=gammaCutEnergy) {
    std::cout<< " *** Primary energy = " << kineticEnergy/geant::MeV
             << " [MeV] is <= gamma production cut = " << gammaCutEnergy/geant::MeV
             << " [MeV] so there is no secondary gamma production at this energy!"
             << std::endl;
    return 0;
  }



  //*******************************************************************************************//
  //************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
  //
  // Create a SeltzerBergerBremsModel model for e-:
  // - Create a Seltzer-Berger bremsstrahlung model
  EMModel *emModel   = nullptr;
  if (bremModelName=="bremSB") {
    emModel = new SeltzerBergerBremsModel(isElectron);
    // - Set low/high energy usage limits to their min/max possible values
    emModel->SetLowEnergyUsageLimit ( 1.0*geant::keV);
    emModel->SetHighEnergyUsageLimit(10.0*geant::GeV);
  } else {
    emModel = new RelativisticBremsModel();
    // - Set low/high energy usage limits to their min/max possible values
    emModel->SetLowEnergyUsageLimit (  1.0*geant::GeV);
    emModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
  }
  emModel->SetUseSamplingTables(!isUseRejection);
  //
  //*******************************************************************************************//

  // check is primary energy is within the usage limits of the model
  if (kineticEnergy<emModel->GetLowEnergyUsageLimit() || kineticEnergy>emModel->GetHighEnergyUsageLimit()) {
    std::cout<< " *** Primary energy = " << kineticEnergy/geant::GeV
             << " [GeV] should be the min/max energy usage limits of the selected model: \n"
             << "   - model name              = " << emModel->GetName() << " \n"
             << "   - low energy usage limit  = " << emModel->GetLowEnergyUsageLimit()/geant::GeV<< " [GeV]\n"
             << "   - high energy usage limit = " << emModel->GetHighEnergyUsageLimit()/geant::GeV<< " [GeV]\n"
             << "  there is no secondary gamma production otherwise!"
             << std::endl;
    return 0;
  }


  //=========== Set the active regions of the model and one physics-parameter object ==========//
  // - Set the model to be active in region index 0
  (emModel->GetListActiveRegions()).resize(1); // one region
  (emModel->GetListActiveRegions())[0] = true; // make it active there
  // - Create one PhysicsParameters object (with defult values)
  PhysicsParameters *physPars = new PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  // - Initialisation of the model
  emModel->Initialize();
  //===========================================================================================//


  //===========================================================================================//
  //== Use the EMModel interface methods of the model to compute some integrated quantities  ==//
  //
  std::cout<< "  "<< matCut->GetMaterial() << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   MaterialCuts: \n"   << matCut;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Particle       =  " << particle->GetName() << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Kinetic energy =  " << kineticEnergy/geant::MeV << "  [MeV] " << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Model name     =  " << emModel->GetName() << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Rejection ?    =  " << isUseRejection << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;  
  // check if we compute atomic-cross section: only for single elemnt materials
  bool isSingleElementMaterial = false;
  if (matCut->GetMaterial()->GetNumberOfElements()==1) {
    isSingleElementMaterial = true;
  }
  //
  // Note: restrictedDEDX and unrestrictedDEDX will be non-zero value only in case of EMModels that has the ComputeDEDX
  //       method implemented i.e. in case of energy loss intercations (i.e. ioni. and brem.)
  double restrictedDEDX          = 0.0;
  double unRestrictedDEDX        = 0.0;
  // Note: atomicCrossSection is computed only in case of materials that has single element
  double atomicCrossSection      = 0.0;
  double macroscopicCrossSection = 0.0;
  //
  // use the model to compute restricted stopping power
  restrictedDEDX   = emModel->ComputeDEDX(matCut, kineticEnergy, particle);
  // use the model to compute un-restricted stopping power
  unRestrictedDEDX = emModel->ComputeDEDX(matCut, kineticEnergy, particle, true);
  // use the model to compute atomic cross section (only in case of single element materials)
  if (isSingleElementMaterial) {
    const Element *elem = (matCut->GetMaterial()->GetElementVector())[0];
    atomicCrossSection  = emModel->ComputeXSectionPerAtom(elem, matCut, kineticEnergy, particle);
  }
  // use the model to compute macroscopic cross section
  macroscopicCrossSection = emModel->ComputeMacroscopicXSection(matCut, kineticEnergy, particle);
  //
  // print out integrated quantities:
  // -atomic cross section
  if (isSingleElementMaterial) {
    std::cout<< "   cross section per atom      :";
    std::cout<< std::setw(14) << std::scientific << std::right << atomicCrossSection/(geant::barn)
             << std::setw(14) << std::left << "     [barn]";
    std::cout<<std::endl;
  }
  //
  // -macroscopic cross section
  std::cout<< "   cross section per volume    :";
  std::cout<< std::setw(14) << std::scientific << std::right << macroscopicCrossSection/(1./geant::cm)
           << std::setw(14) << std::left << "     [1/cm]";
  std::cout<<std::endl;
  //
  // -restricted stopping power
  std::cout<< "   resricted dE/dx  (MeV/cm)   :";
  std::cout<< std::setw(14) << std::scientific << std::right << restrictedDEDX/(geant::MeV/geant::cm)
           << std::setw(14) << std::left << "   [MeV/cm]";
  std::cout<<std::endl;
  //
  // -unrestricted stopping power
  std::cout<< "   unresricted dE/dx (MeV/cm)  :";
  std::cout<< std::setw(14) << std::scientific << std::right << unRestrictedDEDX/(geant::MeV/geant::cm)
           << std::setw(14) << std::left << "   [MeV/cm]";
  std::cout<<std::endl;
  //===========================================================================================//



  //*******************************************************************************************//
  //************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
  //
  std::string hname = "brem_" + bremModelName + "_GV_";
  // energy distribution(k) is sampled in varibale log10(k/primaryEnergy)
  // angular distribution(theta) is sampled in variable log10(1-cos(theta)*0.5)
  //
  // set up a histogram for the secondary gamma energy(k) : log10(k/primaryEnergy)
  double xMin = std::log10(gammaCutEnergy/primaryEnergy);
  double xMax = 0.1;
  Hist *histo_gamma_energy = new Hist(xMin, xMax, numHistBins);
  //
  // set up histogram for the secondary gamma direction(theta) : log10(1-cos(theta)*0.5)
  xMin     = -12.;
  xMax     = 0.5;
  Hist *histo_gamma_angular = new Hist(xMin, xMax, numHistBins);
  //
  // set up a histogram for the post interaction primary e-/e+ energy(E1) : log10(E1/primaryEnergy)
  xMin     = -12;//std::log10(1.-gammaCutEnergy/primaryEnergy);
  xMax     = 0.1;
  Hist *histo_prim_energy = new Hist(xMin, xMax, numHistBins);
  //
  // set up a histogram for the post interaction primary e-/e+ direction(theta) : log10(1-cos(theta)*0.5)
  xMin     = -16.;
  xMax     = 0.5;
  Hist *histo_prim_angular = new Hist(xMin, xMax, numHistBins);

  // start sampling
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Sampling is running : .....................................................      " << std::endl;
  // call sampling method
  double timeInSec = sampleDistribution(numSamples, kineticEnergy, matCut, particle, emModel, histo_gamma_energy,
                                        histo_gamma_angular, histo_prim_energy, histo_prim_angular);
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Time of sampling =  " << timeInSec << " [s]" << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

  std::cout<< "   Writing histograms into files. " << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

  // print out histogram to file: fileName
  char fileName[512];
  sprintf(fileName,"brem_%s_GV_gamma_energy_%s",bremModelName.c_str(),(matCut->GetMaterial()->GetName()).c_str());
  FILE *f     = fopen(fileName,"w");
  Hist *histo = histo_gamma_energy;
  double norm = 0.25/numSamples;
  for (int i=0; i<histo->GetNumBins(); ++i) {
   fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
  }
  fclose(f);
  delete histo;
  //
  sprintf(fileName,"brem_%s_GV_gamma_angular_%s",bremModelName.c_str(),(matCut->GetMaterial()->GetName()).c_str());
  f     = fopen(fileName,"w");
  histo = histo_gamma_angular;
  norm  = 1./numSamples;
  for (int i=0; i<histo->GetNumBins(); ++i) {
   fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
  }
  fclose(f);
  delete histo;
  //
  sprintf(fileName,"brem_%s_GV_%s_energy_%s",bremModelName.c_str(),pname.c_str(),(matCut->GetMaterial()->GetName()).c_str());
  f     = fopen(fileName,"w");
  histo = histo_prim_energy;
  norm  = 0.25/numSamples;
  for (int i=0; i<histo->GetNumBins(); ++i) {
   fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
  }
  fclose(f);
  delete histo;
  //
  sprintf(fileName,"brem_%s_GV_%s_angular_%s",bremModelName.c_str(),pname.c_str(),(matCut->GetMaterial()->GetName()).c_str());
  f     = fopen(fileName,"w");
  histo = histo_prim_angular;
  norm  = 1./numSamples;
  for (int i=0; i<histo->GetNumBins(); ++i) {
   fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
  }
  fclose(f);
  delete histo;
  //*******************************************************************************************//

  // end
  std::cout << "   ================================================================================ "
            << std::endl << std::endl;

  // delete some objects
  delete emModel;

  PhysicsParameters::Clear();
  // clear the ELossTableManager(will alos delete all ELossTable-s) and ELossTableRegister
  ELossTableManager::Instance().Clear();
  ELossTableRegister::Instance().Clear();
  MaterialCuts::ClearAll();
  Material::ClearAllMaterials(); // will delete all Elements and Isotoes as well

return 0;
}


void help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  Model-level GeantV test for testing GeantV e-/e+ models for bremsstrahlung photon emission."
           << std::endl;
  std::cout<<"\n  Usage: bremTest_GV [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}



//*******************************************************************************************//
//************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
//
// implementation of the final state distribution sampling
double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut, Particle *primParticle,
                          EMModel *emModel, Hist *histo1, Hist *histo2, Hist *histo3, Hist *histo4) {
  double ekin       = primaryEnergy;
  double dirx       = 0.0;   // direction
  double diry       = 0.0;
  double dirz       = 1.0;
//  double gamProdCut = matCut->GetProductionCutsInEnergy()[0]; // gamma production threshold
  int    gvcode     = primParticle->GetInternalCode();        // internal code of the primary particle i.e. e-

  // Set up a dummy Geant::GeantTaskData and its geantphysics::PhysicsData member: they are needed in the final state
  // sampling
  Geant::GeantTaskData *td = new Geant::GeantTaskData(1,1);
  PhysicsData *phd = new PhysicsData();
  td->fPhysicsData = phd;
  // Set up a the primary light track for brem.
  LightTrack primaryLT;
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
     int numSecs = emModel->SampleSecondaries(primaryLT, td);
     // get the secondary track i.e. the gamma
     if (numSecs>0) {
       std::vector<LightTrack> &secondaryLT = td->fPhysicsData->GetListOfSecondaries();
       // reduced gamma energy
       double eGamma = secondaryLT[0].GetKinE()/ekin;
       if (eGamma>0.0) {
         histo1->Fill(std::log10(eGamma),1.0);
       }
       double costGamma = secondaryLT[0].GetDirZ();
       costGamma = 0.5*(1.0-costGamma);
       if (costGamma>0.0) {
         costGamma = std::log10(costGamma);
         if (costGamma>-12.) {
           histo2->Fill(costGamma,1.0);
         }
       }
       // go for the post interaction primary
       double ePrim  = primaryLT.GetKinE()/ekin;
       if (ePrim>0.0) {
         ePrim = std::log10(ePrim);
         if (ePrim>-12.0) {
           histo3->Fill(ePrim,1.0);
         }
       }
       double costPrim = primaryLT.GetDirZ();
       costPrim = 0.5*(1.0-costPrim);
       if (costPrim>0.0) {
         costPrim = std::log10(costPrim);
         if (costPrim>-16.) {
           histo4->Fill(costPrim,1.0);
         }
       }
     }
   }
   clock_t end_time = clock();
   return (end_time-start_time)/(double(CLOCKS_PER_SEC));
}
//*******************************************************************************************//
