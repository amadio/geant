#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "GeantRunManager.h"

#include "PhysicsProcessHandler.h"
#include "PhysicsListManager.h"

//FULL-LHCb application
#include "LHCbFullApp.h"
#include "LHCbDetectorConstruction.h"
#include "LHCbParticleGun.h"
#include "LHCbPhysicsList.h"



// some helper methods to get the possible input arguments and configure the user defined components of the application,
// set up the run manager and run the simulation.
void GetArguments(int argc, char *argv[]);
void SetupDetectorConstruction (lhcbapp::LHCbDetectorConstruction* det);
void SetupPrimaryGenerator     (lhcbapp::LHCbParticleGun* gun);
void SetupApplication          (lhcbapp::LHCbFullApp* app);
Geant::GeantRunManager* RunManager();


// The main application: gets the possible input arguments, sets up the run-manager, detector, primary generator,
//                       application and starts the simulation.
int main(int argc, char *argv[]) {
   ROOT::EnableThreadSafety();
  //
  // Read in user arguments
  GetArguments(argc, argv);
  //
  // Create and configure run manager
  Geant::GeantRunManager* runMgr = RunManager();
  
  // Create user defined physics list for the full LHCb application
  geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(new lhcbapp::LHCbPhysicsList());
  
  //
  // Create detector construction
  lhcbapp::LHCbDetectorConstruction *det = new lhcbapp::LHCbDetectorConstruction(runMgr);
  SetupDetectorConstruction(det);
  runMgr->SetDetectorConstruction(det);
  //
  // Create primary generator
  lhcbapp::LHCbParticleGun *gun = new lhcbapp::LHCbParticleGun();
  SetupPrimaryGenerator(gun);
  runMgr->SetPrimaryGenerator(gun);
  //
  // Create application
  lhcbapp::LHCbFullApp *app = new lhcbapp::LHCbFullApp(runMgr, gun);
  SetupApplication(app);
  runMgr->SetUserApplication(app);
  lhcbapp::LHCbParticleGun::Print();
  //
  // Run the simulation
  runMgr->RunSimulation();
  //
  // Delete the run manager
  delete runMgr;

  return 0;
}


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
bool        parConfigIsPerformance    = false;  // run without any user actions
//
//
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
         {"config-run-performance"             , no_argument      , 0, 'r'},

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
  printf("\nUsage: fullLHCbApp [OPTIONS] INPUT_FILE\n\n");
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
        parConfigIsPerformance    = true;
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



Geant::GeantRunManager* RunManager() {
  // create the GeantConfiguration object and the GeantRunManager object
  Geant::GeantConfig*     runConfig  = new Geant::GeantConfig();
  Geant::GeantRunManager* runManager = new Geant::GeantRunManager(parConfigNumPropagators, parConfigNumThreads, runConfig);
  // create the real physics main manager/interface object and set it in the GeantRunManager
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
  //
  // Activate standard scoring
  runConfig->fUseStdScoring = true;

  return runManager;
}


void SetupDetectorConstruction(lhcbapp::LHCbDetectorConstruction* det) {
  if (parDetGDMFile!="") {
    det->SetGDMLFile(parDetGDMFile);
  }
}


void SetupPrimaryGenerator(lhcbapp::LHCbParticleGun* gun) {
  if (parGunPrimaryPerEvent>0)
    gun->SetNumPrimaryPerEvt(parGunPrimaryPerEvent);
  if (parGunPrimaryParticleName!="")
    gun->SetPrimaryName(parGunPrimaryParticleName);
  if (parGunPrimaryKinEnergy>0.)
    gun->SetPrimaryEnergy(parGunPrimaryKinEnergy);
  if ((parGunPrimaryDir[0] || parGunPrimaryDir[1] || parGunPrimaryDir[2]))
    gun->SetPrimaryDirection(parGunPrimaryDir);
}


void SetupApplication(lhcbapp::LHCbFullApp *app) {
  if (parConfigIsPerformance) {
    app->SetPerformanceMode(true);
  }
}
