#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

// GeantV include
#include "Geant/RunManager.h"

// real-physics includes
#include "Geant/PhysicsProcessHandler.h"
#include "Geant/PhysicsListManager.h"
#include "Geant/MSCModel.h"

// application includes
#include "TestEm5.h"
#include "TestEm5DetectorConstruction.h"
#include "TestEm5PrimaryGenerator.h"
#include "TestEm5PhysicsList.h"

#include "Geant/UserFieldConstruction.h"

// some helper methods to get the possible input arguments and configure the user defined components of the application,
// set up the run manager and run the simulation.
void GetInputArguments(int argc, char *argv[]);
void SetupUserPhysicsList(userapplication::TestEm5PhysicsList *physlist);
void SetupUserDetector(userapplication::TestEm5DetectorConstruction *detector);
void SetupUserField(geant::RunManager *runMgr);
void SetupUserPrimaryGenerator(userapplication::TestEm5PrimaryGenerator *primarygun);
void SetupUserApplication(userapplication::TestEm5 *app);
geant::RunManager *RunManager();

// The main application: gets the possible input arguments, sets up the run-manager, physics-list, detector, primary
//                       generator, application and starts the simulation.
int main(int argc, char *argv[])
{
  //
  // Get optional input arguments
  GetInputArguments(argc, argv);
  //
  // Create run manager and set configuration parameters
  geant::RunManager *runManager = RunManager();
  //
  // Create a user defined pysics list object, set its configurable parameters and register it in the global
  // PhysicsListManager
  // NOTE: if the user not register its own physics list then the default physics list will be used
  userapplication::TestEm5PhysicsList *userPhysList = new userapplication::TestEm5PhysicsList("testEm5-user-physics");
  SetupUserPhysicsList(userPhysList);
  geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(userPhysList);
  //
  // create the user detector construction object, set its configurable parameters and register in the RunManager
  userapplication::TestEm5DetectorConstruction *detector = new userapplication::TestEm5DetectorConstruction(runManager);
  SetupUserDetector(detector);
  runManager->SetDetectorConstruction(detector);
  //
  // create the user primary generator object, set its configurable parameters and register in the RunManager
  userapplication::TestEm5PrimaryGenerator *primaryGenerator = new userapplication::TestEm5PrimaryGenerator(detector);
  SetupUserPrimaryGenerator(primaryGenerator);
  runManager->SetPrimaryGenerator(primaryGenerator);
  //
  // set up mag. field propagation
  SetupUserField(runManager);
  //
  // create the testEm5 user application object, set its configurable parameters and register in the RunManager
  userapplication::TestEm5 *testEm5Application = new userapplication::TestEm5(runManager, detector, primaryGenerator);
  SetupUserApplication(testEm5Application);
  runManager->SetUserApplication(testEm5Application);
  //
  // run the simulation
  runManager->RunSimulation();
  //
  // delete the run manager at the end of the simulation
  delete runManager;
  return 0;
}

//
// Optional input arguments that make possible the configuration of detector(parDet), primary generator(parGun), the
// application(parApp), run configuration(parConfig) and some physics processes(parProcess):
//
// detector parameters
std::string parDetTargetMaterial = ""; // i.e. default application value
double parDetTargetThickness     = 0.; // i.e. default application value
double parDetProductionCuts      = 0.; // i.e. default application value
double parDetGammaProdCut        = 0.; // i.e. default application value
double parDetElectronProdCut     = 0.; // i.e. default application value
double parDetPositronProdCut     = 0.; // i.e. default application value
//
// primary generator parameters (primary particle gun)
std::string parGunPrimaryParticleName = ""; // i.e. default application value
double parGunPrimaryKinEnergy         = 0.; // i.e. default application value
//
// TestEm5 application parameters
std::string parAppHist1FileName = "";  // i.e. default application value
double parAppHist1NumBins       = 0.;  // i.e. default application value
double parAppHist1MinVal        = -1.; // i.e. default application value
double parAppHist1MaxVal        = -1.; // i.e. default application value
//
// run configuration parameters
int parConfigNumBufferedEvt   = 4;    // number of events taken to be transported on the same time (buffered)
int parConfigNumRunEvt        = 4000; // total number of events to be transported during the run
int parConfigNumPrimaryPerEvt = 1000; // number of primary particles per event
int parConfigNumThreads       = 4;    // number of working threads
int parConfigNumPropagators   = 1;    // number of propagators per working threads
bool parConfigVectorizedGeom  = 0;    // activate geometry basketizing
int parConfigNumPerBasket     = 16;   // default number of particles per basket
//
// physics process configuration parameters:
std::string parProcessMSCStepLimit = ""; // i.e. default application value
double parProcessStepMaxValue      = 0.; // i.e. default application value

//
// field configuration parameters
int parFieldActive      = 0;         // activate magnetic field
int parFieldUseRK       = 0;         // use Runge-Kutta instead of helix
double parFieldEpsRK    = 0.0003;    // Revised / reduced accuracy - vs. 0.0003 default
int parFieldBasketized  = 1;         // basketize magnetic field
float parFieldVector[3] = {0, 0, 2}; // default constant field value {0,0,2} [kilogauss]

static struct option options[] = {{"det-Target-Material-Name", required_argument, 0, 'a'},
                                  {"det-Target-Thickness", required_argument, 0, 'b'},
                                  {"det-Production-Cuts", required_argument, 0, 'c'},
                                  {"det-Gamma-Production-Cut", required_argument, 0, 'd'},
                                  {"det-Electron-Production-Cut", required_argument, 0, 'e'},
                                  {"det-Positron-Production-Cut", required_argument, 0, 'f'},

                                  {"gun-Primary-Kinetic-Energy", required_argument, 0, 'g'},
                                  {"gun-Primary-Particle-Name", required_argument, 0, 'h'},

                                  {"app-Hist1-File-Name", required_argument, 0, 'i'},
                                  {"app-Hist1-Bin-Number", required_argument, 0, 'j'},
                                  {"app-Hist1-Minimum-Value", required_argument, 0, 'k'},
                                  {"app-Hist1-Maximum-Value", required_argument, 0, 'l'},

                                  {"field-active", required_argument, 0, 'E'},
                                  {"field-vector", required_argument, 0, 'F'},
                                  {"field-use-RK", required_argument, 0, 'G'},
                                  {"field-eps-RK", required_argument, 0, 'H'},
                                  {"field-basketized", required_argument, 0, 'I'},

                                  {"config-Number-Of-Buffered-Events", required_argument, 0, 'M'},
                                  {"config-number-of-buffered-events", required_argument, 0, 'm'},
                                  {"config-Total-Number-Of-Events", required_argument, 0, 'N'},
                                  {"config-total-number-of-events", required_argument, 0, 'n'},
                                  {"config-number-of-primary-per-events", required_argument, 0, 'o'},
                                  {"config-Number-Of-Primary-Per-Events", required_argument, 0, 'O'},
                                  {"config-Number-Of-Threads", required_argument, 0, 'p'},
                                  {"config-Number-Of-Propagators", required_argument, 0, 'q'},
                                  {"config-Vectorized-Geom", required_argument, 0, 'r'},
                                  {"config-Tracks-Per-Basket", required_argument, 0, 's'},

                                  {"process-MSC-Step-Limit", required_argument, 0, 't'},
                                  {"process-Step-Max-Value", required_argument, 0, 'u'},

                                  {"help", required_argument, 0, 'x'},
                                  {0, 0, 0, 0}};

enum DIR_OPTIONS { DIR_X_OPT = 0, DIR_Y_OPT, DIR_Z_OPT };
char *const dir_token[] = {[DIR_OPTIONS::DIR_X_OPT] = (char *const) "x", [DIR_OPTIONS::DIR_Y_OPT] = (char *const) "y",
                           [DIR_OPTIONS::DIR_Z_OPT] = (char *const) "z", NULL};

void help()
{
  printf("\nUsage: runTestEm5 [OPTIONS] INPUT_FILE\n\n");
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

void GetInputArguments(int argc, char *argv[])
{
  int errfnd = 0;
  char *subopts;
  char *value;

  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "", options, &optidx);
    if (c == -1) break;
    switch (c) {
    case 0:
      c = options[optidx].val;
    case 'a':
      parDetTargetMaterial = optarg;
      break;
    case 'b':
      parDetTargetThickness = (double)strtof(optarg, NULL);
      break;
    case 'c':
      parDetProductionCuts = (double)strtof(optarg, NULL);
      break;
    case 'd':
      parDetGammaProdCut = (double)strtof(optarg, NULL);
      break;
    case 'e':
      parDetElectronProdCut = (double)strtof(optarg, NULL);
      break;
    case 'f':
      parDetPositronProdCut = (double)strtof(optarg, NULL);
      break;
    case 'g':
      parGunPrimaryKinEnergy = (double)strtof(optarg, NULL);
      if (parGunPrimaryKinEnergy <= 0) errx(1, "primary particle energy must be positive");
      break;
    case 'h':
      parGunPrimaryParticleName = optarg;
      break;
    case 'i':
      parAppHist1FileName = optarg;
      break;
    case 'j':
      parAppHist1NumBins = (int)strtol(optarg, NULL, 10);
      break;
    case 'k':
      parAppHist1MinVal = (double)strtof(optarg, NULL);
      break;
    case 'l':
      parAppHist1MaxVal = (double)strtof(optarg, NULL);
      break;
    case 'm':
    case 'M':
      parConfigNumBufferedEvt = (double)strtof(optarg, NULL);
      break;
    case 'n':
    case 'N':
      parConfigNumRunEvt = (double)strtof(optarg, NULL);
      break;
    case 'o':
    case 'O':
      parConfigNumPrimaryPerEvt = (double)strtof(optarg, NULL);
      break;
    case 'p':
      parConfigNumThreads = (double)strtof(optarg, NULL);
      break;
    case 'q':
      parConfigNumPropagators = (double)strtof(optarg, NULL);
      break;
    case 'r':
      parConfigVectorizedGeom = (bool)strtof(optarg, NULL);
      break;
    case 's':
      parConfigNumPerBasket = (int)strtol(optarg, NULL, 10);
      break;
    case 't':
      parProcessMSCStepLimit = optarg;
      break;
    case 'u':
      parProcessStepMaxValue = (double)strtof(optarg, NULL);
      break;
    //---- Field
    case 'E':
      parFieldActive = (int)strtol(optarg, NULL, 10);
      break;
    case 'F': // field direction sub-optarg
      subopts = optarg;
      while (*subopts != '\0' && !errfnd) {
        switch (getsubopt(&subopts, dir_token, &value)) {
        case DIR_OPTIONS::DIR_X_OPT:
          parFieldVector[0] = strtod(value, NULL);
          break;
        case DIR_OPTIONS::DIR_Y_OPT:
          parFieldVector[1] = strtod(value, NULL);
          break;
        case DIR_OPTIONS::DIR_Z_OPT:
          parFieldVector[2] = strtod(value, NULL);
          break;
        default:
          fprintf(stderr, "No match found for token: [%s] among DIR_OPTIONS", value);
          errfnd = 1;
          exit(0);
          break;
        }
      }
      break;
    case 'G':
      parFieldUseRK = (int)strtol(optarg, NULL, 10);
      break;
    case 'H':
      parFieldEpsRK = strtod(optarg, NULL);
      break;
    case 'I':
      parFieldBasketized = (int)strtol(optarg, NULL, 10);
      break;
    // -- help
    case 'x':
      help();
      exit(0);
    default:
      help();
      errx(1, "unknown option %c", c);
      exit(0);
    }
  }
}

geant::RunManager *RunManager()
{
  // create the GeantConfiguration object and the RunManager object
  geant::GeantConfig *runConfig = new geant::GeantConfig();
  geant::RunManager *runManager = new geant::RunManager(parConfigNumPropagators, parConfigNumThreads, runConfig);
  // create the real physics main manager/interface object and set it in the RunManager
  runManager->SetPhysicsInterface(new geantphysics::PhysicsProcessHandler());
  //
  // Set parameters of the GeantConfig object:
  runConfig->fNtotal            = parConfigNumRunEvt;
  runConfig->fNbuff             = parConfigNumBufferedEvt;
  runConfig->fNaverage          = parConfigNumPrimaryPerEvt;
  runConfig->fMaxPerBasket      = parConfigNumPerBasket;
  runConfig->fUseVectorizedGeom = parConfigVectorizedGeom;
  //
  // Some additional parameters that have values in this application different than their default
  //
  // this should be true by default from now on since we use only V3
  runConfig->fUseV3         = true;
  runConfig->fNminThreshold = 5 * parConfigNumThreads;
  // Set threshold for tracks to be reused in the same volume
  runConfig->fNminReuse = 100000;
  // Number of steps after the particle is killed:
  // for msc if we run in single scattering setings high number of steps are perfectly possible
  runConfig->fNstepsKillThr = 100000000;
  // Activate standard scoring
  runConfig->fUseStdScoring = true;

  return runManager;
}

void SetupUserDetector(userapplication::TestEm5DetectorConstruction *detector)
{
  if (parDetTargetMaterial != "") detector->SetTargetMaterialName(parDetTargetMaterial);
  if (parDetTargetThickness > 0.) detector->SetTargetThickness(parDetTargetThickness);
  if (parDetProductionCuts > 0.) {
    detector->SetProductionCuts(parDetProductionCuts);
  } else {
    if (parDetGammaProdCut > 0.) detector->SetGammaProductionCut(parDetGammaProdCut);
    if (parDetElectronProdCut > 0.) detector->SetElectronProductionCut(parDetElectronProdCut);
    if (parDetPositronProdCut > 0.) detector->SetPositronProductionCut(parDetPositronProdCut);
  }
}

void SetupUserField(geant::RunManager *runMgr)
{
  auto config = runMgr->GetConfig();
  if (parFieldActive) {
    // Create magnetic field and needed classes for trajectory integration
    auto fieldConstructor = new geant::UserFieldConstruction();
    fieldConstructor->UseConstantMagField(parFieldVector, "kilogauss");

    config->fUseRungeKutta      = parFieldUseRK;
    config->fEpsilonRK          = parFieldEpsRK;
    config->fUseVectorizedField = parFieldBasketized;

    runMgr->SetUserFieldConstruction(fieldConstructor);
    printf("main: Created uniform field and set up field-propagation.\n");
  } else {
    printf("main: no magnetic field configured.\n");
    config->fUseRungeKutta      = false;
    config->fUseVectorizedField = false;
  }
}

void SetupUserPrimaryGenerator(userapplication::TestEm5PrimaryGenerator *primarygun)
{
  // it needs to be consistent with GeantConfig::fNaverage i.e. number of primary particles per event !!!
  if (parConfigNumPrimaryPerEvt > 0) primarygun->SetNumberOfPrimaryParticlePerEvent(parConfigNumPrimaryPerEvt);
  if (parGunPrimaryParticleName != "") primarygun->SetPrimaryParticleName(parGunPrimaryParticleName);
  if (parGunPrimaryKinEnergy > 0.) primarygun->SetPrimaryParticleEnergy(parGunPrimaryKinEnergy);
}

void SetupUserApplication(userapplication::TestEm5 *app)
{
  if (parAppHist1FileName != "") app->SetHist1FileName(parAppHist1FileName);
  if (parAppHist1MinVal > 0.) app->SetHist1Min(parAppHist1MinVal);
  if (parAppHist1MaxVal >= 0.) app->SetHist1Max(parAppHist1MaxVal);
  if (parAppHist1NumBins > 0) app->SetHist1NumBins(parAppHist1NumBins);
}

void SetupUserPhysicsList(userapplication::TestEm5PhysicsList *physlist)
{
  if (parProcessMSCStepLimit != "") {
    if (parProcessMSCStepLimit == "UseSafety") {
      physlist->SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm::kUseSaftey);
    } else if (parProcessMSCStepLimit == "ErrorFree") {
      physlist->SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm::kErrorFree);
    } else if (parProcessMSCStepLimit == "UseDistanceToBoundary") {
      physlist->SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm::kUseDistanceToBoundary);
    } else {
      std::cerr << " **** ERROR runTestEm5::SetupUserPhysicsList() \n"
                << "   unknown MSC stepping algorithm = " << parProcessMSCStepLimit << std::endl;
      exit(-1);
    }
  }
  if (parProcessStepMaxValue > 0.) {
    physlist->SetStepMaxValue(parProcessStepMaxValue);
  }
}
