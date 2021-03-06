#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "Geant/RunManager.h"
#include "Geant/ExternalFramework.h"

#include "Geant/PhysicsProcessHandler.h"
#include "Geant/PhysicsListManager.h"

// Class for constant B-field
#include "Geant/UserFieldConstruction.h"

// FULL-CMS application
#include "CMSFullApp.h"
#include "CMSDetectorConstruction.h"
#include "CMSFieldConstruction.h"
#include "CMSParticleGun.h"
#include "CMSPhysicsList.h"

// some helper methods to get the possible input arguments and configure the user defined components of the application,
// set up the run manager and run the simulation.
void GetArguments(int argc, char *argv[]);
void SetupDetectorConstruction(cmsapp::CMSDetectorConstruction *det);
void SetupFieldConfig(geant::RunManager *runMgr);
void SetupPrimaryGenerator(cmsapp::CMSParticleGun *gun);
void SetupApplication(cmsapp::CMSFullApp *app);
geant::RunManager *RunManager();

//
// Optional input arguments that make possible the configuration of detector(parDet), primary generator(parGun) and
// run configuration(parConfig) :
//
// detector parameters
std::string parDetGDMFile = ""; // i.e. default application values
std::string parFieldFile  = "";
//
// primary generator parameters (primary particle gun)
std::string parGunPrimaryParticleName = "";           // i.e. default application value
int parGunPrimaryPerEvent             = 0;            // i.e. default application value
double parGunPrimaryKinEnergy         = 0.;           // i.e. default application value
double parGunPrimaryDir[3]            = {0., 0., 0.}; // i.e. default application value
//
// run configuration parameters
int parConfigNumBufferedEvt     = 4;  // number of events taken to be transported on the same time (buffered)
int parConfigNumRunEvt          = 10; // total number of events to be transported during the run
int parConfigNumThreads         = 1;  // number of working threads
int parConfigNumPropagators     = 1;  // number of propagators per working threads
int parConfigNumTracksPerBasket = 16; // default number of tracks per basket
int parConfigBsizeFLD           = 0;  // default basket size for field propagator
int parConfigBsizeMSC           = 0;  // default basket size for MSC
int parConfigBsizePHY           = 0;  // default basket size for physics
int parConfigIsPerformance      = 0;  // run without any user actions
int parConfigVectorizedGeom     = 0;  // activate geometry basketizing
int parConfigVectorizedPhysics  = 0;  // activate physics basketizing
int parConfigVectorizedMSC      = 0;  // activate MSC basketizing
int parConfigExternalLoop       = 0;  // activate external loop mode
int parVerboseTracking          = 0;  // verbosity for *every* track

int parConfigMonitoring      = 0; // activate some monitoring
int parConfigSingleTrackMode = 0; // activate single track mode
//
// field configuration parameters
int parFieldType        = 0;      // field type: 0-no field 1-constant field 2-field map
double parFieldEpsRK    = 0.0003; // Revised / reduced accuracy - vs. 0.0003 default
int parFieldBasketized  = 0;      // basketize magnetic field
int parUseRungeKutta    = 0;
float parFieldVector[3] = {0., 0., 38.}; // Constant field value [kiloGauss]

//
//
// The main application: gets the possible input arguments, sets up the run-manager, detector, primary generator,
//                       application and starts the simulation.
int main(int argc, char *argv[])
{
  //
  // Read in user arguments
  GetArguments(argc, argv);
  //
  // Create and configure run manager
  geant::RunManager *runMgr = RunManager();
  //
  // Defining two different physics lists:
  // - ~ 85 - 90 % of the steps are done in region #15 (EcalRegion) out of the 22 regions
  // - a physics list with models using rejection will be used everywhere but region #15
  // - a physics list with models using sampling tables will be used only in region #15
  std::vector<bool> physListActiveRegionList1(22, 1);
  physListActiveRegionList1[15] = 0;
  std::vector<bool> physListActiveRegionList2(22, 0);
  physListActiveRegionList2[15] = 1;
  //
  // Register user defined physics lists for the full CMS application
  // Activating them in different regions - using sampling-table based model only in the "most active" regions
  bool useEMModelsWithSamplingTables = true;
  auto physListRejection =
      new cmsapp::CMSPhysicsList(*runMgr->GetConfig(), "with-rejection", !useEMModelsWithSamplingTables);
  // Switch off physics basketizing (processes except MSC) for the list of models using rejection tables since
  physListRejection->SetBasketizing(false);
  geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(physListRejection, physListActiveRegionList1);

  auto physListTables =
      new cmsapp::CMSPhysicsList(*runMgr->GetConfig(), "with-sampling-tables", useEMModelsWithSamplingTables);
  // This list will use the global parConfigVectorizedPhysics value
  geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(physListTables, physListActiveRegionList2);
  //
  // Create detector construction
  cmsapp::CMSDetectorConstruction *det = new cmsapp::CMSDetectorConstruction(runMgr);
  SetupDetectorConstruction(det);
  runMgr->SetDetectorConstruction(det);
  //
  // Create field  construction  & Get field flags
  //  and activate integration of tracks in field
  SetupFieldConfig(runMgr);

  //
  // Create primary generator
  cmsapp::CMSParticleGun *gun = new cmsapp::CMSParticleGun();
  SetupPrimaryGenerator(gun);
  runMgr->SetPrimaryGenerator(gun);
  //
  // Create application
  cmsapp::CMSFullApp *app = new cmsapp::CMSFullApp(runMgr, gun);
  SetupApplication(app);
  runMgr->SetUserApplication(app);
  cmsapp::CMSParticleGun::Print();
  //
  // Run the simulation
  if (parConfigExternalLoop) {
    userfw::Framework fw(parConfigNumPropagators * parConfigNumThreads, parConfigNumRunEvt, runMgr,
                         runMgr->GetPrimaryGenerator());
    fw.Run();
  } else {
    runMgr->RunSimulation();
    delete runMgr;
  }

  return 0;
}

static struct option options[] = {{"gun-set-primary-energy", required_argument, 0, 'a'},
                                  {"gun-set-primary-type", required_argument, 0, 'b'},
                                  {"gun-set-primary-per-event", required_argument, 0, 'c'},
                                  {"gun-set-primary-direction", required_argument, 0, 'd'},

                                  {"det-set-gdml", required_argument, 0, 'e'},
                                  {"det-set-field", required_argument, 0, 'f'},

                                  {"config-number-of-buffered-events", required_argument, 0, 'm'},
                                  {"config-total-number-of-events", required_argument, 0, 'n'},
                                  {"config-number-of-threads", required_argument, 0, 'p'},
                                  {"config-number-of-propagators", required_argument, 0, 'q'},
                                  {"config-tracks-per-basket", required_argument, 0, 'r'},
                                  {"config-basket-fld", required_argument, 0, 'g'},
                                  {"config-basket-msc", required_argument, 0, 'i'},
                                  {"config-basket-phy", required_argument, 0, 'j'},
                                  {"config-run-performance", required_argument, 0, 's'},
                                  {"config-vectorized-geom", required_argument, 0, 't'},
                                  {"config-external-loop", required_argument, 0, 'u'},
                                  {"config-vectorized-physics", required_argument, 0, 'v'},
                                  {"config-vectorized-MSC", required_argument, 0, 'V'},
                                  {"config-monitoring", required_argument, 0, 'x'},
                                  {"config-single-track", required_argument, 0, 'y'},
                                  {"verbose-tracking", required_argument, 0, 'z'},

                                  {"field-type", required_argument, 0, 'E'},
                                  {"field-vector", required_argument, 0, 'F'},
                                  {"field-use-RK", required_argument, 0, 'G'},
                                  {"field-eps-RK", required_argument, 0, 'I'},
                                  {"field-basketized", required_argument, 0, 'J'},

                                  {"help", no_argument, 0, 'h'},
                                  {0, 0, 0, 0}};

enum PRIMDIR_OPTIONS { PRIMDIR_X_OPT = 0, PRIMDIR_Y_OPT, PRIMDIR_Z_OPT };
char *const primdir_token[] = {[PRIMDIR_OPTIONS::PRIMDIR_X_OPT] = (char *const) "x",
                               [PRIMDIR_OPTIONS::PRIMDIR_Y_OPT] = (char *const) "y",
                               [PRIMDIR_OPTIONS::PRIMDIR_Z_OPT] = (char *const) "z",
                               NULL};

enum MAGFIELD_DIR_OPTIONS { DIR_X_OPT = 0, DIR_Y_OPT, DIR_Z_OPT };
char *const magfield_dir_token[] = {[MAGFIELD_DIR_OPTIONS::DIR_X_OPT] = (char *const) "x",
                                    [MAGFIELD_DIR_OPTIONS::DIR_Y_OPT] = (char *const) "y",
                                    [MAGFIELD_DIR_OPTIONS::DIR_Z_OPT] = (char *const) "z",
                                    NULL};

void help()
{
  printf("\nUsage: fullCMSApp [OPTIONS] INPUT_FILE\n\n");
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

void GetArguments(int argc, char *argv[])
{
  // std::cout << "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" << std::endl;
  // for the set-absorber sub-options
  int errfnd = 0;
  char *subopts;
  char *value;

  while (true) {
    int c, optidx = 0;
    //
    c = getopt_long(argc, argv, "", options, &optidx);
    //
    if (c == -1) break;
    //
    switch (c) {
    case 0:
      c = options[optidx].val; /* fall through */
    //---- Primary generator
    case 'a':
      parGunPrimaryKinEnergy = strtod(optarg, NULL);
      if (parGunPrimaryKinEnergy <= 0) errx(1, "primary particle energy must be positive");
      break;
    case 'b':
      parGunPrimaryParticleName = optarg;
      break;
    case 'c':
      parGunPrimaryPerEvent = (int)strtol(optarg, NULL, 10);
      break;
    case 'd': // primary direction sub-optarg
      subopts = optarg;
      while (*subopts != '\0' && !errfnd) {
        switch (getsubopt(&subopts, primdir_token, &value)) {
        case PRIMDIR_OPTIONS::PRIMDIR_X_OPT:
          parGunPrimaryDir[0] = strtod(value, NULL);
          break;
        case PRIMDIR_OPTIONS::PRIMDIR_Y_OPT:
          parGunPrimaryDir[1] = strtod(value, NULL);
          break;
        case PRIMDIR_OPTIONS::PRIMDIR_Z_OPT:
          parGunPrimaryDir[2] = strtod(value, NULL);
          break;
        default:
          fprintf(stderr, "No match found for token: [%s] among PRIMDIR_OPTIONS", value);
          errfnd = 1;
          exit(0);
          break;
        }
      }
      // std::cout<< " primary dir = (" << parGunPrimaryDir[0] <<", "<<parGunPrimaryDir[1]<<",
      // "<<parGunPrimaryDir[2]<<")" << std::endl;
      break;
    //---- Detector
    case 'e':
      parDetGDMFile = optarg;
      break;
    case 'f':
      parFieldFile = optarg;
      break;

    //---- Run configuration
    case 'g':
      parConfigBsizeFLD = (int)strtol(optarg, NULL, 10);
      break;
    case 'i':
      parConfigBsizeMSC = (int)strtol(optarg, NULL, 10);
      break;
    case 'j':
      parConfigBsizePHY = (int)strtol(optarg, NULL, 10);
      break;
    case 'm':
      parConfigNumBufferedEvt = (int)strtol(optarg, NULL, 10);
      break;
    case 'n':
      parConfigNumRunEvt = (int)strtol(optarg, NULL, 10);
      break;
    case 'p':
      parConfigNumThreads = (int)strtol(optarg, NULL, 10);
      std::cout << "  Argument read: ConfigNumThreads = " << parConfigNumThreads << std::endl;
      break;
    case 'q':
      parConfigNumPropagators = (int)strtol(optarg, NULL, 10);
      break;
    case 'r':
      parConfigNumTracksPerBasket = (int)strtol(optarg, NULL, 10);
      break;
    case 's':
      parConfigIsPerformance = (int)strtol(optarg, NULL, 10);
      break;
    case 't':
      parConfigVectorizedGeom = (int)strtol(optarg, NULL, 10);
      break;
    case 'v':
      parConfigVectorizedPhysics = (int)strtol(optarg, NULL, 10);
      break;
    case 'x':
      parConfigMonitoring = (int)strtol(optarg, NULL, 10);
      break;
    case 'y':
      parConfigSingleTrackMode = (int)strtol(optarg, NULL, 10);
      break;
    case 'V':
      parConfigVectorizedMSC = (int)strtol(optarg, NULL, 10);
      break;
    case 'u':
      parConfigExternalLoop = (int)strtol(optarg, NULL, 10);
      break;
    case 'z':
      parVerboseTracking = (int)strtol(optarg, NULL, 10);
      // std::cout << "  Argument read: VerboseTracking = " << parVerboseTracking << std::endl;
      break;
    //---- Field
    case 'E':
      parFieldType = (int)strtol(optarg, NULL, 10);
      break;
    case 'F': // field direction sub-optarg
      subopts = optarg;
      while (*subopts != '\0' && !errfnd) {
        switch (getsubopt(&subopts, magfield_dir_token, &value)) {
        case MAGFIELD_DIR_OPTIONS::DIR_X_OPT:
          parFieldVector[0] = strtof(value, NULL);
          break;
        case MAGFIELD_DIR_OPTIONS::DIR_Y_OPT:
          parFieldVector[1] = strtof(value, NULL);
          break;
        case MAGFIELD_DIR_OPTIONS::DIR_Z_OPT:
          parFieldVector[2] = strtof(value, NULL);
          break;
        default:
          fprintf(stderr, "No match found for token: [%s] among MAGFIELD_DIR_OPTIONS", value);
          errfnd = 1;
          exit(0);
          break;
        }
      }
      break;
    case 'G':
      parUseRungeKutta = (int)strtol(optarg, NULL, 10);
      break;
    case 'I':
      parFieldEpsRK = strtod(optarg, NULL);
      break;
    case 'J':
      parFieldBasketized = (int)strtol(optarg, NULL, 10);
      break;
    //---- Help
    case 'h':
      help();
      exit(0);
      break;
    default:
      help();
      errx(1, "unknown option %c", c);
    }
  }
}

geant::RunManager *RunManager()
{
  // create the GeantConfiguration object and the RunManager object
  geant::GeantConfig *runConfig = new geant::GeantConfig();
  std::cout << " Instantiation RunManager with : \n"
            << "    # Threads     = " << parConfigNumThreads << std::endl
            << "    # Propagators = " << parConfigNumPropagators << std::endl;
  geant::RunManager *runManager = new geant::RunManager(parConfigNumPropagators, parConfigNumThreads, runConfig);
  //
  // Set parameters of the GeantConfig object:
  runConfig->fNtotal = parConfigNumRunEvt;
  runConfig->fNbuff  = parConfigNumBufferedEvt;
  //  runConfig->fNaverage      = parConfigNumPrimaryPerEvt;
  //
  // Some additional parameters that have values in this application different than their default
  //
  // this should be true by default from now on since we use only V3
  runConfig->fUseV3         = true;
  runConfig->fNminThreshold = 5 * parConfigNumThreads;
  // Set threshold for tracks to be reused in the same volume
  runConfig->fNminReuse  = 100000;
  runConfig->fNperBasket = parConfigNumTracksPerBasket;
  // Activate monitoring
  if (parConfigMonitoring > 0) runConfig->fMonHandlers = true;
  // Activate vectorized geometry
  runConfig->fUseVectorizedGeom = parConfigVectorizedGeom;
  if (parConfigVectorizedGeom == 2) runConfig->fUseSDGeom = true;

  runConfig->fUseVectorizedPhysics = parConfigVectorizedPhysics;
  runConfig->fNvecPHY              = parConfigBsizePHY;
  if (parConfigVectorizedPhysics == 2) runConfig->fUseSDPhysics = true;

  runConfig->fUseVectorizedMSC = parConfigVectorizedMSC;
  runConfig->fNvecMSC          = parConfigBsizeMSC;
  if (parConfigVectorizedMSC == 2) runConfig->fUseSDMSC = true;

  // create the real physics main manager/interface object and set it in the RunManager
  runManager->SetPhysicsInterface(new geantphysics::PhysicsProcessHandler(*runConfig));

  runConfig->fUseStdScoring   = false;
  runConfig->fSingleTrackMode = parConfigSingleTrackMode;

  return runManager;
}

void SetupDetectorConstruction(cmsapp::CMSDetectorConstruction *det)
{
  if (parDetGDMFile != "") {
    det->SetGDMLFile(parDetGDMFile);
  }
}

void SetupFieldConfig(geant::RunManager *runMgr)
{
  geant::cxx::UserFieldConstruction *fieldCtion = nullptr;
  auto config                                   = runMgr->GetConfig();
  bool useField                                 = false;
  std::cout << "=== Field type: ";
  switch (parFieldType) {
  case 0:
    std::cout << "no magnetic field\n";
    break;
  case 1:
    std::cout << "uniform magnetic field of (" << parFieldVector[0] << ", " << parFieldVector[1] << ", "
              << parFieldVector[2] << ") kilogauss\n";
    useField = true;
    {
      // fieldCtion = new CMSFieldConstruction(useUniformField);  // Default value
      auto uniformFldCtion = new geant::UserFieldConstruction();
      uniformFldCtion->UseConstantMagField(parFieldVector, "kilogauss");
      fieldCtion = uniformFldCtion;
    }
    break;
  case 2:
    std::cout << "field from sampled CMS field map\n";
    useField = true;
    {
      auto cmsFieldCtion = new CMSFieldConstruction();
      if (parFieldFile != "") cmsFieldCtion->SetFileForField(parFieldFile);
      fieldCtion = cmsFieldCtion;
    }
    break;
  default:
    std::cout << "Unknown field option. Defaulting to no field.\n";
    parFieldType = 0;
  }

  config->fUseRungeKutta = parUseRungeKutta;

  runMgr->SetUserFieldConstruction(fieldCtion);

  if (useField) {
    // Create magnetic field and needed classes for trajectory integration
    config->fEpsilonRK          = parFieldEpsRK;
    config->fUseVectorizedField = parFieldBasketized;
    config->fNvecFLD            = parConfigBsizeFLD;
    if (parFieldBasketized == 2) config->fUseSDField = true;
    std::cout << "=== Created magnetic field and set up field-propagation.\n";
  } else {
    config->fUseVectorizedField = false;
    std::cout << "=== No magnetic field configured.\n";
    runMgr->SetUserFieldConstruction(nullptr);
  }
}

void SetupPrimaryGenerator(cmsapp::CMSParticleGun *gun)
{
  if (parGunPrimaryPerEvent > 0) gun->SetNumPrimaryPerEvt(parGunPrimaryPerEvent);
  if (parGunPrimaryParticleName != "") gun->SetPrimaryName(parGunPrimaryParticleName);
  if (parGunPrimaryKinEnergy > 0.) gun->SetPrimaryEnergy(parGunPrimaryKinEnergy);
  if ((parGunPrimaryDir[0] || parGunPrimaryDir[1] || parGunPrimaryDir[2])) gun->SetPrimaryDirection(parGunPrimaryDir);
}

void SetupApplication(cmsapp::CMSFullApp *app)
{
  if (parConfigIsPerformance) {
    app->SetPerformanceMode(true);
  }
}
