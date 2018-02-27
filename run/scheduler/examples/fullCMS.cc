#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "Geant/RunManager.h"
#include "Geant/ExternalFramework.h"

#include "PhysicsProcessHandler.h"
#include "PhysicsListManager.h"

//FULL-CMS application
#include "CMSFullApp.h"
#include "CMSDetectorConstruction.h"
#include "CMSParticleGun.h"
#include "CMSPhysicsList.h"



// some helper methods to get the possible input arguments and configure the user defined components of the application,
// set up the run manager and run the simulation.
void GetArguments(int argc, char *argv[]);
void SetupDetectorConstruction (cmsapp::CMSDetectorConstruction* det);
void SetupPrimaryGenerator     (cmsapp::CMSParticleGun* gun);
void SetupApplication          (cmsapp::CMSFullApp* app);
geant::RunManager* RunManager();

//
// Optional input arguments that make possible the configuration of detector(parDet), primary generator(parGun) and
// run configuration(parConfig) :
//
// detector parameters
std::string parDetGDMFile             = "";   // i.e. default application values
//
// primary generator parameters (primary particle gun)
std::string parGunPrimaryParticleName =  "";  // i.e. default application value
int         parGunPrimaryPerEvent     =   0;  // i.e. default application value
double      parGunPrimaryKinEnergy    =  0.;  // i.e. default application value
double      parGunPrimaryDir[3]       = {0.,0.,0.}; // i.e. default application value
//
// run configuration parameters
int         parConfigNumBufferedEvt   =     4;  // number of events taken to be transported on the same time (buffered)
int         parConfigNumRunEvt        =    10;  // total number of events to be transported during the run
int         parConfigNumThreads       =     4;  // number of working threads
int         parConfigNumPropagators   =     1;  // number of propagators per working threads
int         parConfigNumTracksPerBasket =  16;  // default number of tracks per basket
int         parConfigIsPerformance    =     0;  // run without any user actions
int         parConfigVectorizedGeom   =     0;  // activate geometry basketizing
int         parConfigExternalLoop     =     0;  // activate external loop mode
//
//
// The main application: gets the possible input arguments, sets up the run-manager, detector, primary generator,
//                       application and starts the simulation.
int main(int argc, char *argv[]) {
  //
  // Read in user arguments
  GetArguments(argc, argv);
  //
  // Create and configure run manager
  geant::RunManager* runMgr = RunManager();
  //
  // Create user defined physics list for the full CMS application
  geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(new cmsapp::CMSPhysicsList());
  //
  // Create detector construction
  cmsapp::CMSDetectorConstruction *det = new cmsapp::CMSDetectorConstruction(runMgr);
  SetupDetectorConstruction(det);
  runMgr->SetDetectorConstruction(det);
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
    userfw::Framework fw(parConfigNumPropagators*parConfigNumThreads, parConfigNumRunEvt, runMgr, runMgr->GetPrimaryGenerator());
    fw.Run();
  } else {
    runMgr->RunSimulation();
    delete runMgr;
  }

  return 0;
}

static struct option options[] = {
         {"gun-set-primary-energy"             , required_argument, 0, 'a'},
				 {"gun-set-primary-type"               , required_argument, 0, 'b'},
         {"gun-set-primary-per-event"          , required_argument, 0, 'c'},
         {"gun-set-primary-direction"          , required_argument, 0, 'd'},

         {"det-set-gdml"                       , required_argument, 0, 'e'},

         {"config-number-of-buffered-events"   , required_argument, 0, 'm'},
         {"config-total-number-of-events"      , required_argument, 0, 'n'},
         {"config-number-of-threads"           , required_argument, 0, 'p'},
         {"config-number-of-propagators"       , required_argument, 0, 'q'},
         {"config-tracks-per-basket"           , required_argument, 0, 'r'},
         {"config-run-performance"             , required_argument, 0, 's'},
         {"config-vectorized-geom"             , required_argument, 0, 't'},
         {"config-external-loop"               , required_argument, 0, 'u'},

         {"help", no_argument, 0, 'h'},
         {0, 0, 0, 0}
};

enum PRIMDIR_OPTIONS { PRIMDIR_X_OPT = 0, PRIMDIR_Y_OPT, PRIMDIR_Z_OPT};
char *const primdir_token[] = {
   [PRIMDIR_OPTIONS::PRIMDIR_X_OPT]   = (char* const)"x",
   [PRIMDIR_OPTIONS::PRIMDIR_Y_OPT]   = (char* const)"y",
   [PRIMDIR_OPTIONS::PRIMDIR_Z_OPT]   = (char* const)"z",
   NULL
};

void help() {
  printf("\nUsage: fullCMSApp [OPTIONS] INPUT_FILE\n\n");
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}


void GetArguments(int argc, char *argv[]) {
  //std::cout << "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" << std::endl;
  // for the set-absorber sub-options
  int    errfnd  =  0;
  char  *subopts;
  char  *value;

  while (true) {
    int c, optidx=0;
    //
    c = getopt_long(argc, argv, "", options, &optidx);
    //
    if (c==-1)
      break;
    //
    switch (c) {
      case 0:
        c = options[optidx].val; /* fall through */
      //---- Primary generator
      case 'a':
        parGunPrimaryKinEnergy = strtod(optarg, NULL);
        if (parGunPrimaryKinEnergy<=0)
          errx(1, "primary particle energy must be positive");
        break;
      case 'b':
        parGunPrimaryParticleName = optarg;
        break;
      case 'c':
        parGunPrimaryPerEvent     = (int)strtol(optarg, NULL, 10);
        break;
      case 'd':  // primary direction sub-optarg
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
       //std::cout<< " primary dir = (" << parGunPrimaryDir[0] <<", "<<parGunPrimaryDir[1]<<", "<<parGunPrimaryDir[2]<<")" << std::endl;
       break;
      //---- Detector
      case 'e':
        parDetGDMFile = optarg;
    	  break;
      //---- Run configuration
      case 'm':
        parConfigNumBufferedEvt   = (int)strtol(optarg, NULL, 10);
        break;
      case 'n':
        parConfigNumRunEvt        = (int)strtol(optarg, NULL, 10);
        break;
      case 'p':
        parConfigNumThreads       = (int)strtol(optarg, NULL, 10);
        break;
      case 'q':
        parConfigNumPropagators   = (int)strtol(optarg, NULL, 10);
        break;
      case 'r':
        parConfigNumTracksPerBasket   = (int)strtol(optarg, NULL, 10);
        break;
      case 's':
        parConfigIsPerformance    = (int)strtol(optarg, NULL, 10);
        break;
      case 't':
        parConfigVectorizedGeom   = (int)strtol(optarg, NULL, 10);
        break;
      case 'u':
        parConfigExternalLoop     = (int)strtol(optarg, NULL, 10);
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



geant::RunManager* RunManager() {
  // create the GeantConfiguration object and the RunManager object
  geant::GeantConfig*     runConfig  = new geant::GeantConfig();
  geant::RunManager* runManager = new geant::RunManager(parConfigNumPropagators, parConfigNumThreads, runConfig);
  // create the real physics main manager/interface object and set it in the RunManager
  runManager->SetPhysicsInterface(new geantphysics::PhysicsProcessHandler());
  //
  // Set parameters of the GeantConfig object:
  runConfig->fNtotal        = parConfigNumRunEvt;
  runConfig->fNbuff         = parConfigNumBufferedEvt;
//  runConfig->fNaverage      = parConfigNumPrimaryPerEvt;
  //
  // Some additional parameters that have values in this application different than their default
  //
  // this should be true by default from now on since we use only V3
  runConfig->fUseV3         = true;
  runConfig->fNminThreshold = 5*parConfigNumThreads;
  // Set threshold for tracks to be reused in the same volume
  runConfig->fNminReuse     = 100000;
  runConfig->fNperBasket  = parConfigNumTracksPerBasket;
  //
  // Activate vectorized geometry
  runConfig->fUseVectorizedGeom = parConfigVectorizedGeom;

  return runManager;
}


void SetupDetectorConstruction(cmsapp::CMSDetectorConstruction* det) {
  if (parDetGDMFile!="") {
    det->SetGDMLFile(parDetGDMFile);
  }
}


void SetupPrimaryGenerator(cmsapp::CMSParticleGun* gun) {
  if (parGunPrimaryPerEvent>0)
    gun->SetNumPrimaryPerEvt(parGunPrimaryPerEvent);
  if (parGunPrimaryParticleName!="")
    gun->SetPrimaryName(parGunPrimaryParticleName);
  if (parGunPrimaryKinEnergy>0.)
    gun->SetPrimaryEnergy(parGunPrimaryKinEnergy);
  if ((parGunPrimaryDir[0] || parGunPrimaryDir[1] || parGunPrimaryDir[2]))
    gun->SetPrimaryDirection(parGunPrimaryDir);
}


void SetupApplication(cmsapp::CMSFullApp *app) {
  if (parConfigIsPerformance) {
    app->SetPerformanceMode(true);
  }
}
