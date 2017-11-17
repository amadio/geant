/**
 * @brief  elasticTest_GV: GeantV-Real-physics test for testing models for gamma photons elastic into e-/e+ pairs.
 *
 * @author W Pokorski
 * @date   June 2017
 *
 * Run ./elasticTest_GV --help for more details!
 * The corresponding quantities/distributions can be obtained by using the elasticTest_G4 Geant4 test.
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
#include "Isotope.h"

// vecgeom includes
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"

#include "Region.h"
#include "PhysicsParameters.h"

#include "ELossTableManager.h"
#include "ELossTableRegister.h"

#include "Particle.h"
#include "Proton.h"
#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "DiffuseElasticModel.h"

#include "LightTrack.h"
#include "PhysicsData.h"

// a simple histogram class: can be changed later
#include "Hist.h"

using geantphysics::Material;
using geantphysics::Element;
using geantphysics::Isotope;

using geantphysics::PhysicsParameters;

using geantphysics::Particle;
using geantphysics::Proton;
using geantphysics::Electron;
using geantphysics::Positron;
using geantphysics::Gamma;

// the hadronic model
using geantphysics::HadronicFinalStateModel;
using geantphysics::DiffuseElasticModel;

using geantphysics::LightTrack;
using geantphysics::PhysicsData;

using userapplication::Hist;

//
// default values of the input parameters
static std::string   materialName("NIST_MAT_C");         // material is lead
static int           Z = 6;
static int           N = 12;
static int           numHistBins       = 100;             // number of histogram bins between min/max values
static double        numSamples        = 1.e+7;           // number of required final state samples
static double        primaryEnergy     = 1.0;             // primary particle energy in [GeV]

static struct option options[] = {
  {"material-name     (with a NIST_MAT_ prefix; see more in material doc.)     - default: NIST_MAT_Pb"        , required_argument, 0, 'm'},
  {"primary-energy    (in internal energy units i.e. [GeV])                    - default: 0.1"                , required_argument, 0, 'E'},
  {"number-of-samples (number of required final state samples)                 - default: 1.e+7"              , required_argument, 0, 'f'},
  {"number-of-bins    (number of bins in the histogram)                        - default: 100"                , required_argument, 0, 'n'},
  {"help"                                                                                                     , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
void help();

//***********************************************************************************************//
//***** THIS WILL BE MODEL SPECIFIC: contains the final state sampling and hist. building  ******//
// method to create photon energy distribution using a SeltzerBergerBremsModel as input argument
double sampleDistribution(double numSamples, double primaryEnergy, Isotope *isotope,
                          Particle *primParticle, HadronicFinalStateModel *elModel, Hist *h1);
//***********************************************************************************************//


//===========================================================================================//
//
int main(int argc, char *argv[]) {
  //
  //============================== Get input parameters =====================================//
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "h:m:E:f:n:", options, &optidx);
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
  //===========================================================================================//

  //============================= Set user defined input data =================================//


  Isotope* isotope = Isotope::GetIsotope( Z, N );

  std::cout << "mass " << isotope->GetIsoMass() << std::endl;

  // Create target material: which is supposed to be a NIST Material
  Material *matDetector = Material::NISTMaterial(materialName);
  //
  // Create primnary gamma particle
  Particle *particle = Proton::Definition();
  //
  // Set particle kinetic energy
  double kineticEnergy    = primaryEnergy;
  //

  //============= Initialization i.e. building up and init the physics ========================//
  // Create a dummy vecgeom::geometry:
  //  - with only one volume i.e. world
  //
  // create a vecgeom::LogicalVolume
  vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1.,1.,1.);
  vecgeom::LogicalVolume worldl(&worldParams);

  // create a dummy vecgeom::media for this volume:
  // Only the vecgeom::material index (that will be 0) is used later to identify our material during the init.
  worldl.SetMaterialPtr((void*)matDetector);
  //  worldl.SetTrackingMediumPtr(new vecgeom::Medium((matDetector->GetName()).c_str(),new vecgeom::Material(),new double()));
  vecgeom::GeoManager::Instance().SetWorld(worldl.Place());
  vecgeom::GeoManager::Instance().CloseGeometry();

  //===========================================================================================//


  //*******************************************************************************************//
  //************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
  //
  // Create DiffuseElasticModel model
  HadronicFinalStateModel *elModel = new DiffuseElasticModel();
  // - Set low/high energy usage limits
  elModel->SetLowEnergyUsageLimit (100.0*geant::eV);
  elModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
  //
  //*******************************************************************************************//

  // - Initialisation of the model
  elModel->Initialize();
  //===========================================================================================//


  //===========================================================================================//
  //
  // Get the MaterialCuts of the target: we have only one
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Isotope: \n"   << *isotope  << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Particle       =  " << particle->GetName() << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Kinetic energy =  " << kineticEnergy/geant::MeV << "  [MeV] " << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Model name     =  " << elModel->GetName() << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  //

  //*******************************************************************************************//
  //************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
  //
  // Create a histogram to store the t
  //
  double xMin =  0.0;
  double xMax =  1.0; //
  Hist *histo_t    = new Hist(xMin, xMax, numHistBins);

  //
  //
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Sampling is running : .....................................................      " << std::endl;
  // call sampling method
  double timeInSec = sampleDistribution(numSamples, kineticEnergy, isotope, particle, elModel, histo_t);
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  std::cout<< "   Time of sampling =  " << timeInSec << " [s]" << std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

  // print out histogram to file: fileName
  char fileName[512];
  sprintf(fileName,"t_%s",(isotope->GetName()).c_str());
  FILE *f     = fopen(fileName,"w");
  double norm = 1./numSamples;
  Hist *histo = histo_t;
  for (int i=0; i<histo->GetNumBins(); ++i) {
   fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
  }
  delete histo;
  fclose(f);
  //
  std::cout<< "   Histogram is written  into files ................................................"<< std::endl;
  std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
  //*******************************************************************************************//

  // end
  std::cout << "   ================================================================================ "
            << std::endl << std::endl;

  // delete the EMModel objects
  delete elModel;

  PhysicsParameters::Clear();
  Material::ClearAllMaterials(); // will delete all Elements and Isotoes as well

return 0;
}


void help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  Model-level GeantV test for testing GeantV Bethe-Heitler model for e-/e+"
           <<"  pair production by photons."
           << std::endl;
  std::cout<<"\n  Usage: elasticTest_GV [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}



//*******************************************************************************************//
//************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
//
// implementation of the final state distribution sampling
double sampleDistribution(double numSamples, double primaryEnergy, Isotope *isotope,
                          Particle *primParticle, HadronicFinalStateModel *elModel, Hist *h1) {
  double ekin       = primaryEnergy;
  double dirx       = 0.0;   // direction
  double diry       = 0.0;
  double dirz       = 1.0;
  int    gvcode     = primParticle->GetInternalCode();  // internal code of the primary particle
  double mass       = primParticle->GetPDGMass();

  // Set up a dummy Geant::GeantTaskData and its geantphysics::PhysicsData member: they are needed in the final state
  // sampling
  Geant::GeantTaskData *td = new Geant::GeantTaskData(1,1);
  PhysicsData         *phd = new PhysicsData();
  td->fPhysicsData         = phd;
  // Set up a the primary light track for brem.
  LightTrack primaryLT;
  // init time
  clock_t  start_time = clock();
  for (long int i=0; i<numSamples; ++i) {
     // we will use members:
     //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
     //  fGVcode       <==>  fGVcode   // internal particle code
     //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
     //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
     //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
     //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
     primaryLT.SetKinE(ekin);
     primaryLT.SetGVcode(gvcode);
//       primaryLT.SetTrackIndex(0); // not important now
     primaryLT.SetDirX(dirx);
     primaryLT.SetDirY(diry);
     primaryLT.SetDirZ(dirz);
     primaryLT.SetMass(mass);
//       primaryLT.SetTotalMFP(1.0); // not important now
     //
     // clean the number of secondary tracks used (in PhysicsData)
     td->fPhysicsData->SetNumUsedSecondaries(0);
     //
     // invoke the interaction
     elModel->SampleFinalState(primaryLT, isotope, td);

     // show updated primaryLT
     h1->Fill(std::acos(primaryLT.GetDirZ()));

   }
   clock_t end_time = clock();
   return (end_time-start_time)/(double(CLOCKS_PER_SEC));
}
//*******************************************************************************************//
