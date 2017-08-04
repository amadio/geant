///////////////////////////////////////////////////////////////////////////////////////////
////
////                      caloAppRP.cc
////                      Created: 1 August 2017
////                      Author: Ryan Schmitz
////
//// Description: A (linear) calorimeter implemented using VecGeom libraries. This executable handles input
//			arguments provided by a user (see list below), and passes arguments between
//			classes. These arguments can be provided either through a macro (see caloApp.mac)
//			or directly from terminal commands.
//
//			Note: In order to add more than 3 absorbers per layer, one only needs to copy the
//			associated input arguments provided here and increment the number associated with the
//			new absorber (e.g. copy argument det-absorber-3-material to det-absorber-4-material
//			and change every 3 to a 4). Other parts of the application have been written to
//			automatically account for this potential change.
////
////////////////////////////////////////////////////////////////////////////////////////////

#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "GeantRunManager.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

// realphysics
#include "PhysicsProcessHandler.h"
#include "PhysicsListManager.h"
#include "MSCModel.h"

// application
#include "CaloApp.h"
#include "CaloDetectorConstruction.h"
#include "CaloPrimaryGenerator.h"
#include "UserPhysicsList.h"
#ifdef GEANT_TBB
#include "TaskMgrTBB.h"
#endif
//using namespace Geant;

void GetArguments(int argc, char *argv[]);
void SetupUserPhysicsList(userapplication::UserPhysicsList *physlist);
void SetupUserDetector        (userapplication::CaloDetectorConstruction*   detector);
void SetupUserPrimaryGenerator(userapplication::CaloPrimaryGenerator*     primarygun);
void SetupUserApplication     (userapplication::CaloApp*                         app);
void PrintRunInfo();
Geant::GeantRunManager* RunManager();

int main(int argc, char *argv[]) {

  // Read in user arguments
  GetArguments(argc, argv);

  // Create and configure run manager
  Geant::GeantRunManager* runMgr = RunManager();

  // Create user defined physics list
  userapplication::UserPhysicsList *userPhysList = new userapplication::UserPhysicsList("CaloPhysicsList");
  SetupUserPhysicsList(userPhysList);
  geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(userPhysList);

  // Create Calo detector construction
  userapplication::CaloDetectorConstruction *calo = new userapplication::CaloDetectorConstruction(runMgr);
  SetupUserDetector(calo);
  runMgr->SetDetectorConstruction(calo);

  // Create Calo primary generator
  userapplication::CaloPrimaryGenerator *caloGun = new userapplication::CaloPrimaryGenerator(calo);
  SetupUserPrimaryGenerator(caloGun);
  runMgr->SetPrimaryGenerator(caloGun);

  // Create the real physics Calo application
  userapplication::CaloApp *caloApplication = new userapplication::CaloApp(runMgr,calo,caloGun);
  SetupUserApplication(caloApplication);
  runMgr->SetUserApplication(caloApplication);
  
  // Print basic parameters for the simulation
  PrintRunInfo();

  // Run the simulation
  runMgr->RunSimulation();

  // Delete the run manager
  delete runMgr;
  return 0; 
}

int 		n_buffered    = 4;
int 		n_events      = 10;
double 		n_primaries_per_event = 100.; // average number of tracks per event
int 		n_threads     = 4;
int 		n_propagators = 1;
int 		n_track_max   = 500;
int 		n_learn_steps = 0;

int 		n_layers      = 10;
int 		n_absorbers   = 2;
const int 	maxAbs        = 5;
double 		caloYZ 	     = 10; //10cm
std::string 	absorberMat[maxAbs+1];
int 		absorberThickness[maxAbs+1];
bool 		userMat[maxAbs+1];
bool 		userThickness[maxAbs+1];

std::string 	histName="CaloHist";
double 		histMin	      = 0;
double 		histMax       = 1;
int 		histNumBins   = 100;

double 		prodCutLength = 0;
double 		prodCutEnergy = 0;
double 		prodCutGamma  = 0;
double 		prodCutElectron = 0;
double 		prodCutPositron = 0;

std::string 	particleProcessMSCStepLimit;
int 		particleProcessStepMaxValue = 0;

bool 		performance = true;
bool 		monitor = false, debug = false, coprocessor = false, tbbmode = false, userLayers=false, userAbsorbers=false, userYZ=false;

double 		primary_energy = 0.1; // [GeV]
std::string 	primary_type  = "e-";

Geant::TaskBroker *broker = nullptr;

static struct option options[] = {{"gun-primary-energy", required_argument, 0, 'a'},
				  {"gun-primary-type", required_argument, 0, 'b'},

				  {"det-yzLength", required_argument, 0, 'c'},
                                  {"det-numLayers", required_argument, 0, 'd'},
                                  {"det-numAbsorbers", required_argument, 0, 'e'},
                                  {"det-absorber-1-material", required_argument, 0, 'f'},
                                  {"det-absorber-2-material", required_argument, 0, 'g'},
                                  {"det-absorber-3-material", required_argument, 0, 'i'},
                                  {"det-absorber-1-thickness", required_argument, 0, 'j'},
                                  {"det-absorber-2-thickness", required_argument, 0, 'k'},
                                  {"det-absorber-3-thickness", required_argument, 0, 'l'},
                                  
                                  {"config-number-events", required_argument, 0, 'm'},
                                  {"config-number-primaries-per-event", required_argument, 0, 'n'},
                                  {"config-number-max-tracks-per-basket", required_argument, 0, 'o'},
                                  {"config-number-threads", required_argument, 0, 'p'},
                                  {"config-number-propagators", required_argument, 0, 'q'},
                                  {"config-number-buffered-events", required_argument, 0, 'r'},
				  {"config-flag-monitor", no_argument, 0, 's'},
                                  {"config-flag-debug", no_argument, 0, 't'},
                                  {"config-flag-coprocessor", no_argument, 0, 'u'},
                                  {"config-flag-tbbmode", no_argument, 0, 'v'},

				  {"hist-name", required_argument, 0, 'w'},
				  {"hist-bin-min", required_argument, 0, 'x'},
				  {"hist-bin-max", required_argument, 0, 'y'},
				  {"hist-bin-number", required_argument, 0, 'z'},

				  {"det-prod-cut-length", required_argument, 0, 'A'},
				  {"det-prod-cut-energy", required_argument, 0, 'B'},
				  {"det-prod-cut-gamma", required_argument, 0, 'C'},
				  {"det-prod-cut-electron", required_argument, 0, 'D'},
				  {"det-prod-cut-positron", required_argument, 0, 'E'},

				  {"particle-process-MSC-step-limit", required_argument, 0, 'F'},
				  {"particle-process-step-max-value", required_argument, 0, 'G'},

                                  {"help", no_argument, 0, 'h'},
                                  {0, 0, 0, 0}};

 
void help() {
  printf("\nUsage: runAppRP [OPTIONS] INPUT_FILE\n\n");

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

void PrintRunInfo(){

  // Print run information
    std::cout<< "\n\n"
             << " =============================================================================== \n"
             << "  primary GV code      : " << primary_type << "  [22 => e-; 23 => e+; 42 => gamma]\n"
             << "  primary energy       : " << primary_energy  << " [GeV] \n"
             << "  #events              : " << n_events << " \n"
             << "  #primaries per event : " << n_primaries_per_event << " \n"
             << "  total # primaries    : " << n_events*n_primaries_per_event << " \n"
             << " ===============================================================================\n\n";
}

void GetArguments(int argc, char *argv[]) {
  std::cout << "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" << std::endl;

  while (true) {
    int c, optidx=0;

    c = getopt_long(argc, argv, "a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:A:B:C:D:E:F:G:h", options, &optidx);

    if (c == -1)
      break;

    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'a':
      primary_energy = (double)strtof(optarg, NULL);

      if (primary_energy <= 0)
        errx(1, "primary particle energy must be positive");
      break;

    case 'b':
      primary_type = optarg;
      break;
  
    case 'c':
	userYZ=true;
	caloYZ=(double)strtof(optarg,NULL);
	break;

    case 'd':
      userLayers=true;
      n_layers = (int)strtol(optarg, NULL, 10);

	if (n_layers <= 0)
	  errx(1,"number of layers must be positive");
	break;

    case 'e':
      userAbsorbers=true;
      n_absorbers = (int)strtol(optarg, NULL, 10);
	if (n_absorbers <= 0)
	  errx(1,"number of absorbers must be positive");
	else if (n_absorbers>3)
	  errx(1,"greater than 3 absorbers not supported under current settings. Add more arguments to macro and try again.");
	break;

    case 'f':
	userMat[1]=true;	
	absorberMat[1]=optarg;
        break;

    case 'g':
	userMat[2]=true;
	absorberMat[2]=optarg;
        break;

    case 'i':
	if (n_absorbers<3)
		errx(1,"ERROR: Can't edit absorber 3 since there are less than 3 absorbers");
	userMat[3]=true;
	absorberMat[3]=optarg;	
        break;


    case 'j':
	userThickness[1]=true;
	absorberThickness[1]=(double)strtof(optarg,NULL);
	if (absorberThickness[1]<=0)
		errx(1,"absorber thickness must be positive");
	break;
    
    case 'k':
	userThickness[2]=true;
	absorberThickness[2]=(double)strtof(optarg,NULL);
	if (absorberThickness[2]<=0)
		errx(1,"absorber thickness must be positive");
	break;

    case 'l':
	userThickness[3]=true;
	if (n_absorbers<3)
		errx(1,"ERROR: Can't edit absorber 3 since there are less than 3 absorbers");
	absorberThickness[3]=(double)strtof(optarg,NULL);
	if (absorberThickness[3]<=0)
		errx(1,"absorber thickness must be positive");
	break;

    case 'm':
      n_events = (int)strtol(optarg, NULL, 10);

      if (n_events <= 0)
        errx(1, "number of events must be positive");
      break;

    case 'n':
      n_primaries_per_event = (int)strtol(optarg, NULL,10);

      if (n_primaries_per_event <= 0)
        errx(1, "number of primary tracks per event must be positive");
      break;

    case 'o':
      n_track_max = (int)strtol(optarg, NULL, 10);

      if (n_track_max < 1)
        errx(1, "max number of tracks per basket must be positive");
      break;

    case 'p':
      n_threads = (int)strtol(optarg, NULL, 10);

      if (n_threads < 1)
        errx(1, "number of threads must be positive");
      break;

    case 'q':
      n_propagators = (int)strtol(optarg, NULL, 10);
      break;

    case 'r':
      n_buffered = (int)strtol(optarg, NULL, 10);

      if (n_buffered < 1)
        errx(1, "number of buffered events must be positive");
      break;

    case 's':
      monitor = true;
      break;

    case 't':
      debug = true;
      break;

    case 'u':
      coprocessor = true;
      break;

    case 'v':
      tbbmode = true;
      break;

    case 'w':
      histName = optarg;
      break;

    case 'x':
      histMin=(double)strtof(optarg,NULL);
      break;
    
    case 'y':
      histMax=(double)strtof(optarg,NULL);
      break;

    case 'z':
      histNumBins = (int)strtol(optarg, NULL, 10);

      if (histNumBins < 0)
        errx(1, "number of bins must be positive");
      break;
    
    case 'A':
      prodCutLength=(double)strtof(optarg,NULL);
      break;
    
    case 'B':
      prodCutEnergy=(double)strtof(optarg,NULL);
      break;
    
    case 'C':
      prodCutGamma=(double)strtof(optarg,NULL);
      break;
    
    case 'D':
      prodCutElectron=(double)strtof(optarg,NULL);
      break;
    
    case 'E':
      prodCutPositron=(double)strtof(optarg,NULL);
      break;
    
    case 'F':
      particleProcessMSCStepLimit = optarg;
      break;

    case 'G':
      particleProcessStepMaxValue = (int)strtol(optarg, NULL, 10);

      if (histNumBins < 0)
        errx(1, "number of bins must be positive");
      break;
    
    case 'h':
      help();
      exit(0);
      break;

    default:
      help();
      errx(1, "unknown option %c", c);
    }
  }
  if (n_absorbers>=3){
	if(!userMat[3] || !userThickness[3]){
	std::cerr << "Error: defined a calorimeter with 3 or more absorbers, but didn't define the absorbers fully! Set the material and thickness and try again\n";
  	}	
  }

  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32,128,1);
    broker = gpuBroker;
    nthreads += gpuBroker->GetNstream()+1;
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }
}


Geant::GeantRunManager* RunManager() {
  Geant::GeantConfig* config=new Geant::GeantConfig();
  config->fNtotal        = n_events;
  config->fNbuff         = n_buffered;
  config->fUseMonitoring = monitor;
  config->fNminThreshold = 5*n_threads;
  config->fNmaxBuffSpill = 128;   // New configuration parameter!!!
  config->fUseV3         = true;
  config->SetMonitored(Geant::GeantConfig::kMonQueue, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonMemory, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonBasketsPerVol, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonVectors, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonConcurrency, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonTracksPerEvent, monitor);
  config->fNaverage = n_primaries_per_event;   // Average number of primaries per event

  //
  // for msc if we run in single scattering setings:
  config->fNstepsKillThr = 100.*config->fNstepsKillThr;

  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  config->fPriorityThr = 0.05;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  config->fNperBasket = 16;   // Initial vector size (tunable)

  // This is now the most important parameter for memory considerations
  config->fMaxPerBasket = n_track_max;   // Maximum vector size (tunable)
//  config->fEmin = 3.E-6; // [3 KeV] energy cut : not used in real physics
//  config->fEmax = primaryEnergy;  // [100 GeV] NOTE: probably it is not needed anymore with the real physics

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
  config->fLearnSteps = n_learn_steps;
  if (performance) config->fLearnSteps = 0;
   // Activate I/O
  config->fFillTree = false;
   // Activate debugging using -DBUG_HUNT=ON in your cmake build
  if (debug) {
    config->fUseDebug = true;
    config->fDebugTrk = 1;
  }
// Activate standard scoring
  config->fUseStdScoring = true;
  if (performance) config->fUseStdScoring = false;
  // Monitor the application
  config->fUseAppMonitoring = false;

  // Set threshold for tracks to be reused in the same volume
  config->fNminReuse = 100000;

  // Create run manager
  Geant::GeantRunManager *runManager = new Geant::GeantRunManager(n_propagators, n_threads, config);
  if (broker) runManager->SetCoprocessorBroker(broker);
  // Create the (new) real physics main manager object and set it in the GeantPropagator
  runManager->SetPhysicsInterface( new geantphysics::PhysicsProcessHandler() );
#ifdef GEANT_TBB
  if (tbbmode)
    runManager->SetTaskMgr( new TaskMgrTBB() );
#endif

  return runManager;
}

void SetupUserDetector(userapplication::CaloDetectorConstruction* det) {
  if (userLayers) det->SetNumLayers(n_layers);
  if (userAbsorbers) det->SetNumAbsorbers(n_absorbers);
  if (userYZ) det->SetDetectorYZ(caloYZ);
  if (prodCutLength>0) det->SetProductionCutsByLength(prodCutLength);
  if (prodCutEnergy>0) det->SetProductionCutsByEnergy(prodCutEnergy);
  if (prodCutGamma>0) det->SetProductionCutsByLength(prodCutGamma);
  if (prodCutElectron>0) det->SetProductionCutsByLength(prodCutElectron);
  if (prodCutPositron>0) det->SetProductionCutsByLength(prodCutPositron);

  for (int k=1; k<=n_absorbers; k++) {
        if (userMat[k])
                det->SetAbsorberMaterialName(k,absorberMat[k]);
        if (userThickness[k]){
                det->SetAbsorberThickness(k,absorberThickness[k]*geant::mm);
        }
  }
}

void SetupUserPrimaryGenerator(userapplication::CaloPrimaryGenerator* primarygun) {
  primarygun->SetNumberOfPrimaryParticlePerEvent(n_primaries_per_event);
  primarygun->SetPrimaryParticleName(primary_type);
  primarygun->SetPrimaryParticleEnergy(primary_energy);
}

void SetupUserApplication(userapplication::CaloApp* caloApp) {
  caloApp->SetHist1FileName(histName);
  if(histMin!=0) caloApp->SetHist1Min(histMin);
  if(histMax!=0) caloApp->SetHist1Max(histMax);
  if(histNumBins!=0) caloApp->SetHist1NumBins(histNumBins);
}

void SetupUserPhysicsList(userapplication::UserPhysicsList *userPhysList){
  if (particleProcessMSCStepLimit!="") {
    if (particleProcessMSCStepLimit=="UseSafety") {
      userPhysList->SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm::kUseSaftey);
    } else if (particleProcessMSCStepLimit=="ErrorFree") {
      userPhysList->SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm::kErrorFree);
    } else if (particleProcessMSCStepLimit=="UseDistanceToBoundary") {
      userPhysList->SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm::kUseDistanceToBoundary);
    } else {
      std::cerr<< " **** ERROR caloAppRP::SetupUserPhysicsList() \n"
               << "   unknown MSC stepping algorithm = " << particleProcessMSCStepLimit
               << std::endl;
      exit(-1);
    }
  }
    if (particleProcessStepMaxValue>0.) {
      userPhysList->SetStepMaxValue(particleProcessStepMaxValue);
    }
}

