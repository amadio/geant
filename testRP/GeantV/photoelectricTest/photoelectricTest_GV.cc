/**
 * @brief  photoelectricTest_GV: GeantV-Real-physics application to test models for gamma's photoelectric.
 * @author M Bandieramonte from M Novak example
 * @date   April 2017
 *
 * A GeantV-Real-physics application to test models for gamma's photoelectric available in the
 * given physics list: sauterGavrilaPhotoElectric (for the moment) The models are extracted from
 * the given physics list. Note, that the corresponding Geant4 application does not use a
 * physics-list to build-up and initialise the models but instead it creates and initialises
 * the models directly.
 * Both the energy and angular distributions of the emitted photoelectron can be tested
 * with the application (see the available options below by running with --help option)
 *
 * Run ./photoelectricTest_GV --help    for more details!
 * The corresponding Geant4 test application is photoelectricTest_G4.
 */



#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <ctime>
#include <sstream>

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
#include "SauterGavrilaPhotoElectricModel.h"

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

// the photoelectric model
using geantphysics::EMModel;
using geantphysics::SauterGavrilaPhotoElectricModel;

using geantphysics::LightTrack;
using geantphysics::PhysicsData;

using userapplication::Hist;

//
// default values of the input parameters
static std::string   particleName("gamma");                  // primary particle is gamma
static std::string   materialName("NIST_MAT_Pb");         // material is lead
static std::string   photoElectricModelName("SauterGavrilaPhotoElectric");      // name of the sauterGavrilaPhotoElectric model to test
static int           numHistBins       = 100;             // number of histogram bins between min/max values
static double        numSamples        = 1.e+07;           // number of required final state samples - 1.e+07;
static double        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
static double        prodCutValue      = 0.1;             // by default in length and internal units i.e. [cm]
static bool          isProdCutInLength = true;            // is the production cut value given in length ?
static bool          isAlias           = false;            // is the Alias sampling active ?

static struct option options[] = {

    {"particle-name     (possible particle names: gamma)                         - default: gamma"               , required_argument, 0, 'p'},
    {"material-name     (with a NIST_MAT_ prefix; see more in material doc.)     - default: NIST_MAT_Pb"         , required_argument, 0, 'm'},
    {"primary-energy    (in internal energy units i.e. [GeV])                    - default: 0.01"                , required_argument, 0, 'E'},
    {"number-of-samples (number of required final state samples)                 - default: 1.e+7"               , required_argument, 0, 'f'},
    {"number-of-bins    (number of bins in the histogram)                        - default: 100"                 , required_argument, 0, 'n'},
    {"model-name        (sauterGavrilaPhotoElectric)                             - default: sauterGavrilaPhotoElectric"       , required_argument, 0, 'b'},
    {"isAlias           (is the Alias sampling active ?)                         - default: false"               , required_argument, 0, 's'},
    {"cut-value         (secondary production threshold value for all particles) - default: 0.1"                 , required_argument, 0, 'c'},
    {"cut-in-energy     (is the production cut value given in energy ? )         - default: false"               , no_argument      , 0, 'e'},
    {"help"                                                                                                      , no_argument      , 0, 'h'},

    {0, 0, 0, 0}
};
void help();

//***********************************************************************************************//
//***** THIS WILL BE MODEL SPECIFIC: contains the final state sampling and hist. building  ******//
// method to create photon energy distribution using a SautergavrilaPhotoElectricModel as input argument
double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut,
                          Particle *primParticle, EMModel *model, Hist *histo1, Hist *histo2, Hist *histo3);

//***********************************************************************************************//

double CalculateDiffCrossSection(double tau, double cosTheta);


//===========================================================================================//
//
int main(int argc, char *argv[]) {
    //
    //============================== Get input parameters =====================================//
    while (true) {
        int c, optidx = 0;
        c = getopt_long(argc, argv, "eh:m:E:f:n:c:p:b:s:", options, &optidx);
        if (c == -1)
            break;
        switch (c) {
            case 0:
                c = options[optidx].val;
                /* fall through */
            case 's':
                isAlias = true;
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
            case 'p':
                particleName = optarg;
                break;
            case 'b':
                photoElectricModelName = optarg;
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

    //if (!(photoElectricModelName=="photoElectricSB" || photoElectricModelName=="photoElectricRel")) {
    //  std::cout << "  *** unknown photoElectric. model name = " << photoElectricModelName << std::endl;
    //  help();
    //  return 0;
    //}

    // Create primary particle
    Particle *particle = Gamma::Definition();
    //bool         isElectron = true;
    std::string  pname;



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
    //double gammaCutEnergy = matCut->GetProductionCutsInEnergy()[0];
    /*if (kineticEnergy<=gammaCutEnergy) {
     std::cout<< " *** Primary energy = " << kineticEnergy/geant::MeV
     << " [MeV] is <= gamma production cut = " << gammaCutEnergy/geant::MeV
     << " [MeV] so there is no secondary gamma production at this energy!"
     << std::endl;
     return 0;
     }*/



    //*******************************************************************************************//
    //************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
    //
    // Create a SauterGavrilaPhotoElectricModel model for gammas:
    // - Create a SauterGavrilaPhotoElectricModel model
    std::cout<<"Creating the model SauterGavrilaPhotoElectricModel\n";
    EMModel *emModel = new SauterGavrilaPhotoElectricModel(photoElectricModelName, true); //true to use Alias Sampling method
    EMModel *emModel_rej = new SauterGavrilaPhotoElectricModel(photoElectricModelName, false); //true to use Alias Sampling method
    // - Set low/high energy usage limits to their min/max possible values
    emModel->SetLowEnergyUsageLimit ( 0.01*geant::keV);

    emModel->SetHighEnergyUsageLimit(100.0*geant::GeV);


    emModel_rej->SetLowEnergyUsageLimit ( 0.01*geant::keV);

    emModel_rej->SetHighEnergyUsageLimit(100.0*geant::GeV);
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

    //=========== Set the active regions of the model and one physics-parameter object ==========//
    // - Set the model to be active in region index 0
    (emModel_rej->GetListActiveRegions()).resize(1); // one region
    (emModel_rej->GetListActiveRegions())[0] = true; // make it active there
    // - Create one PhysicsParameters object (with defult values)
    PhysicsParameters *physPars_rej = new PhysicsParameters();
    // - Set it to be active in region index 0
    (physPars_rej->GetListActiveRegions()).resize(1);
    (physPars_rej->GetListActiveRegions())[0] = true;
    // - Initialisation of the model
    emModel_rej->Initialize();
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
    std::cout<< "   Alias sampling =  " << isAlias << std::endl;
    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

    // check if we compute atomic-cross section: only for single elemnt materials
    bool isSingleElementMaterial = false;
    if (matCut->GetMaterial()->GetNumberOfElements()==1) {
        isSingleElementMaterial = true;
    }
    //


    // Note: atomicCrossSection is computed only in case of materials that has single element
    double atomicCrossSection      = 0.0;
    double macroscopicCrossSection = 0.0;

    //
    // use the model to compute atomic cross section (only in case of single element materials)
    std::cout<< std::setw(14) << std::scientific <<std::endl;
    if (isSingleElementMaterial) {

        const Element *elem = (matCut->GetMaterial()->GetElementVector())[0];
        std::cout<<std::setw(14)<<std::scientific<<std::endl;
        atomicCrossSection  = emModel->ComputeXSectionPerAtom(elem, matCut, kineticEnergy, particle);

    }
    
    //UNCOMMENT TO TEST CrossSectionPerVolume method
    //clock_t  start = clock();
    //for (int i= 0; i<stat; i++)
    // use the model to compute macroscopic cross section
    macroscopicCrossSection = emModel->ComputeMacroscopicXSection(matCut, kineticEnergy, particle);
    //clock_t  end = clock();
    //std::cout<<"ComputeMacroscopicXSection ex-time: "<<(end-start)/(double(CLOCKS_PER_SEC))<<std::endl;
    
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
    //===========================================================================================//



    //*******************************************************************************************//
    //************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
    //


    //std::string hname = photoElectricModelName + "_GV_" + str;
    // energy distribution(k) is sampled in varibale log10(k/primaryEnergy)
    // angular distribution(theta) is sampled in variable log10(1-cos(theta)*0.5)
    //
    // set up a histogram for the secondary photoelectron energy(k) : log10(k/primaryEnergy)
    double xMin = -3.0; //std::log10(primaryEnergy);//std::log10(gammaCutEnergy/primaryEnergy);
    double xMax = 3.0;
    Hist *histo_photoelectron_energy = new Hist(xMin, xMax, numHistBins);
    Hist *histo_photoelectron_energy_rej = new Hist(xMin, xMax, numHistBins);
    //
    // set up histogram for the secondary photoelectron direction(theta) : log10(1-cos(theta)*0.5)
    xMin     = -12.;
    xMax     = 0.5;
    Hist *histo_photoelectron_angular = new Hist(xMin, xMax, numHistBins);
    Hist *histo_photoelectron_angular_rej = new Hist(xMin, xMax, numHistBins);

    xMin     = -1.0001;
    xMax     = 1.0001;
    Hist *histo_angle = new Hist(xMin, xMax, numHistBins);
    Hist *histo_angle_rej = new Hist(xMin, xMax, numHistBins);



    //THE PRIMARY DOES NOT SURVIVE TO THE P.E. PROCESS SO WE DON'T NEED OTHER HISTOGRAMS - THEY MUST BE EMPTY
    ////The primary does not survive to photoelectric effect, so these other histo are not needed
    // set up a histogram for the post interaction primary gamma energy(E1) : log10(E1/primaryEnergy)
    //xMin     = -12;//std::log10(1.-gammaCutEnergy/primaryEnergy);
    //xMax     = 0.1;
    //Hist *histo_prim_energy = new Hist(xMin, xMax, numHistBins);
    //
    // set up a histogram for the post interaction primary e-/e+ direction(theta) : log10(1-cos(theta)*0.5)
    //xMin     = -16.;
    //xMax     = 0.5;
    //Hist *histo_prim_angular = new Hist(xMin, xMax, numHistBins);

    // start sampling
    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
    std::cout<< "   Sampling is running : .....................................................      " << std::endl;


    // call sampling method
    double timeInSec = sampleDistribution(numSamples, kineticEnergy, matCut, particle, emModel, histo_photoelectron_energy, histo_photoelectron_angular, histo_angle);
    double timeInSec_rej = sampleDistribution(numSamples, kineticEnergy, matCut, particle, emModel_rej, histo_photoelectron_energy_rej, histo_photoelectron_angular_rej, histo_angle_rej);
    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
    std::cout<< "   Time of sampling Alias =  " << timeInSec << " [s]" << std::endl;
    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
    std::cout<< "   Time of sampling Rejection =  " << timeInSec_rej << " [s]" << std::endl;
    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

    std::cout<< "   Writing histograms into files. " << std::endl;
    std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;

    // print out histogram to file: fileName
    char fileName[512];
    std::ostringstream strs;
    strs << primaryEnergy*1000;
    std::string str = strs.str();

    sprintf(fileName,"photoElectric_%s_GV_photoelectron_energy_%s_%sMeV.ascii",photoElectricModelName.c_str(),(matCut->GetMaterial()->GetName()).c_str(), str.c_str());
    FILE *f     = fopen(fileName,"w");
    Hist *histo = histo_photoelectron_energy;
    double norm = 0.25/numSamples;
    for (int i=0; i<histo->GetNumBins(); ++i) {
        fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
        //printf("%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
    }
    fclose(f);
    delete histo;

    //
    sprintf(fileName,"photoElectric_%s_GV_photoelectron_angular_%s_%sMeV.ascii",photoElectricModelName.c_str(),(matCut->GetMaterial()->GetName()).c_str(), str.c_str());
    f     = fopen(fileName,"w");
    histo = histo_photoelectron_angular;
    norm  = 1./numSamples;
    for (int i=0; i<histo->GetNumBins(); ++i) {
        fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
    }
    fclose(f);
    delete histo;


    sprintf(fileName,"photoElectric_%s_GV_photoelectron_angular_rejection_%s_%sMeV.ascii",photoElectricModelName.c_str(),(matCut->GetMaterial()->GetName()).c_str(), str.c_str());
    f     = fopen(fileName,"w");
    histo = histo_photoelectron_angular_rej;
    norm  = 1./numSamples;
    for (int i=0; i<histo->GetNumBins(); ++i) {
        fprintf(f,"%d\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm);
    }
    fclose(f);
    delete histo;

    double xsec[numHistBins];
    double cosTheta;

    //sprintf(fileName,"GV_%s_cosTheta_%s_%sMeV.ascii",photoElectricModelName.c_str(),(matCut->GetMaterial()->GetName()).c_str(), str.c_str());
    sprintf(fileName,"cosTheta_%sMeV.ascii", str.c_str());
    f     = fopen(fileName,"w");
    histo = histo_angle;
    norm  = 1./numSamples;
    double sum=0;
    for (int i=0; i<histo->GetNumBins(); ++i) {
        cosTheta=histo->GetX()[i]+0.5*histo->GetDelta();//+0.5*deltaTheta;
        xsec[i]= CalculateDiffCrossSection(kineticEnergy/geant::kElectronMassC2, cosTheta);
        sum+=xsec[i];
    }

    for (int i=0; i<histo->GetNumBins(); ++i) {
        //fprintf(f,"%d\t%.8g\t%.8g\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i]*norm, xsec[i]/sum,(histo->GetY()[i]*norm)/(xsec[i]/sum) );
        fprintf(f,"%d\t%.8g\t%.8g\t%.8g\t%.8g\n",i,histo->GetX()[i]+0.5*histo->GetDelta(),histo->GetY()[i], histo_angle_rej->GetY()[i], ( (histo->GetY()[i])/(histo_angle_rej->GetY()[i]) ) );

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
    std::cout<<"  Model-level GeantV test for testing GeantV e-/e+ models for photoElectric e- emission."
    << std::endl;
    std::cout<<"\n  Usage: photoElectricTest_GV [OPTIONS] \n"<<std::endl;
    for (int i = 0; options[i].name != NULL; i++) {
        printf("\t-%c  --%s\n", options[i].val, options[i].name);
    }
    std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}

double CalculateDiffCrossSection(double tau, double cosTheta)
{

    // Based on Geant4 : G4SauterGavrilaAngularDistribution
    // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
    // input  : energy0  (incoming photon energy)
    // input  : cosTheta (cons(theta) of photo-electron)
    // output : dsigma   (differential cross section, K-shell only)

    //double tau = energy0 / geant::kElectronMassC2;

    //gamma and beta: Lorentz factors of the photoelectron
    double gamma = tau + 1.0;
    double beta = std::sqrt(tau * (tau + 2.0)) / gamma;

    double z = 1 - beta * cosTheta;
    double z2 = z * z;
    double z4 = z2 * z2;
    double y = 1 - cosTheta * cosTheta; //sen^2(theta)

    double dsigma = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
    return dsigma;

}


//*******************************************************************************************//
//************                 THIS CONTAINS MODEL SPECIFIC PARTS                 ***********//
//
// implementation of the final state distribution sampling

double sampleDistribution(double numSamples, double primaryEnergy, const MaterialCuts *matCut, Particle *primParticle,
                          EMModel *emModel, Hist *histo1, Hist *histo2, Hist *histo3) {


    double ekin       = primaryEnergy;
    double dirx       = 0.0;   // direction
    double diry       = 0.0;
    double dirz       = 1.0;
    int    gvcode     = primParticle->GetInternalCode();        // internal code of the primary particle i.e. e-

    // Set up a dummy Geant::GeantTaskData and its geantphysics::PhysicsData member: they are needed in the final state
    // sampling
    Geant::GeantTaskData *td = new Geant::GeantTaskData(1,1);
    PhysicsData *phd = new PhysicsData();
    td->fPhysicsData = phd;
    // Set up a the primary light track for photoElectric.
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
            double ePhotoElectron = secondaryLT[0].GetKinE();
            if (ePhotoElectron>0.0) {
                ePhotoElectron = log10((ekin - ePhotoElectron)*1000000);
                histo1->Fill(ePhotoElectron,1.0);

            }
            double costPhotoElectron = secondaryLT[0].GetDirZ();
            //if(costPhotoElectron>1) costPhotoElectron=1;
            //else if (costPhotoElectron<-1) costPhotoElectron=-11;
            //histo3->Fill(costPhotoElectron, 1.0);
            costPhotoElectron = 0.5*(1.0-costPhotoElectron);
            if (costPhotoElectron>0.0) {
                costPhotoElectron = std::log10(costPhotoElectron);
                if (costPhotoElectron>-12.) {
                    histo2->Fill(costPhotoElectron,1.0);
                }
            }
        }
    }


    clock_t end_time = clock();
    return (end_time-start_time)/(double(CLOCKS_PER_SEC));
}
//*******************************************************************************************//
