#ifdef USE_ROOT
#include "Rtypes.h"
#include "TGeoManager.h"
#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "GeantRunManager.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

// realphysics
#include "ParticleGun.h"
#include "PhysicsProcessHandler.h"

#include "ExN03ApplicationRP.h"

#ifdef GEANT_TBB
#include "TaskMgrTBB.h"
#endif

//using namespace Geant;

static int n_events      = 20;
static int n_buffered    = 4;
static int n_threads     = 4;
static int n_track_max   = 500;
static int n_learn_steps = 0;
static int n_reuse       = 100000;
static int n_propagators = 1;
static bool monitor = false, score = false, debug = false, coprocessor = false, tbbmode = false, usev3 = true;
static double n_avrg_tracks_per_evt = 500.; // average number of tracks per event
static double primary_energy        = 100; // [GeV]

static struct option options[] = {{"primary-energy", required_argument, 0, 'E'},
                                  {"events", required_argument, 0, 'e'},
                                  {"avrg-tracks-per-evt", required_argument, 0, 'a'},
                                  {"geometry", required_argument, 0, 'g'},
                                  {"learn-steps", required_argument, 0, 'l'},
                                  {"max-tracks-per-basket", required_argument, 0, 'B'},
                                  {"monitor", no_argument, 0, 'm'},
                                  {"debug", no_argument, 0, 'd'},
                                  {"nbuffer", required_argument, 0, 'b'},
                                  {"score", no_argument, 0, 's'},
                                  {"threads", required_argument, 0, 't'},
                                  {"coprocessor", required_argument, 0, 'r'},
                                  {"tbbmode", required_argument, 0, 'i'},
                                  {"reuse", required_argument, 0, 'u'},
                                  {"propagators", required_argument, 0, 'p'},
                                  {"v2", no_argument, 0, 'v'},
                                  {"help", no_argument, 0, 'h'},
                                  {0, 0, 0, 0}};

void help() {
  printf("\nUsage: runAppRP [OPTIONS] INPUT_FILE\n\n");

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

int main(int argc, char *argv[]) {
  std::cout << "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" << std::endl;
  std::string exn03_geometry_filename("ExN03.root");

  if (argc == 1) {
    help();
    exit(0);
  }

  while (true) {
    int c, optidx = 0;

    c = getopt_long(argc, argv, "hvE:e:a:g:l:B:b:t:r:i:u:p:", options, &optidx);

    if (c == -1)
      break;

    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'E':
      primary_energy = (double)strtof(optarg, NULL);

      if (primary_energy <= 0)
        errx(1, "primary particle energy must be positive");
      break;

    case 'e':
      n_events = (int)strtol(optarg, NULL, 10);

      if (n_events <= 0)
        errx(1, "number of events must be positive");
      break;

    case 'a':
      n_avrg_tracks_per_evt = (double)strtof(optarg, NULL);

      if (n_avrg_tracks_per_evt <= 0)
        errx(1, "average number of primary tracks per event must be positive");
      break;

    case 'g':
      exn03_geometry_filename = optarg;
      break;

    case 'l':
      n_learn_steps = (int)strtol(optarg, NULL, 10);

      if (n_learn_steps <= 0)
        errx(1, "number of learning steps must be positive");

      break;

    case 'B':
      n_track_max = (int)strtol(optarg, NULL, 10);

      if (n_track_max < 1)
        errx(1, "max number of tracks per basket must be positive");

      break;

    case 'm':
      monitor = true;
      break;

    case 'b':
      n_buffered = (int)strtol(optarg, NULL, 10);

      if (n_buffered < 1)
        errx(1, "number of buffered events must be positive");
      break;

    case 't':
      n_threads = (int)strtol(optarg, NULL, 10);

      if (n_threads < 1)
        errx(1, "number of threads must be positive");
      break;

    case 's':
      score = true;
      break;

    case 'r':
      coprocessor = optarg;
      break;

    case 'i':
      tbbmode = true;
      break;

    case 'u':
      n_reuse = (int)strtol(optarg, NULL, 10);
      break;

    case 'p':
      n_propagators = (int)strtol(optarg, NULL, 10);
      break;

    case 'v':
      usev3 = false;
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
  bool performance = true;

  Geant::TaskBroker *broker = nullptr;
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

  //
  // Parameters for the primary generator
  int     gvParticleCode =  22;    // internal code of the primary particle:  22 ==> e-
  double  primaryEnergy  = primary_energy;          // kinetic energy of the primary particles in [GeV]
  double  avNPrimPerEvt  = n_avrg_tracks_per_evt;   // Average number of tracks per event
  double  xPos           =  -8.;   // x-position of the particle gun
  double  yPos           =   0.;   // y-position of the particle gun
  double  zPos           =   0.;   // z-position of the particle gun
  double  xDir           =   1.;   // x-direction of the particle gun
  double  yDir           =   0.;   // x-direction of the particle gun
  double  zDir           =   0.;   // x-direction of the particle gun


  Geant::GeantConfig* config=new Geant::GeantConfig();

  config->fGeomFileName  = exn03_geometry_filename;
  config->fNtotal        = n_events;
  config->fNbuff         = n_buffered;
  config->fUseMonitoring = monitor;
  config->fNminThreshold = 5*n_threads;
  config->fNmaxBuffSpill = 128;   // New configuration parameter!!!
  config->fUseV3         = usev3;
  config->SetMonitored(Geant::GeantConfig::kMonQueue, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonMemory, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonBasketsPerVol, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonVectors, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonConcurrency, monitor);
  config->SetMonitored(Geant::GeantConfig::kMonTracksPerEvent, monitor);
  config->fNaverage = avNPrimPerEvt;   // Average number of tracks per event

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
  config->fNminReuse = n_reuse;

  // Create run manager
  Geant::GeantRunManager *runMgr = new Geant::GeantRunManager(n_propagators, n_threads, config);
  if (broker) runMgr->SetCoprocessorBroker(broker);
  // Create the (new) real physics main manager object and set it in the GeantPropagator
  runMgr->SetPhysicsInterface( new geantphysics::PhysicsProcessHandler() );

  // for vector physics -OFF now
  // runMgr->SetVectorPhysicsProcess(new GVectorPhysicsProcess(config->fEmin, nthreads));

  // Create the real physics ParticleGun object and set it in the run manager
  runMgr->SetPrimaryGenerator( new geantphysics::ParticleGun(avNPrimPerEvt, gvParticleCode, primaryEnergy, xPos, yPos,
                                                             zPos, xDir, yDir, zDir) );
  // Create the real physics ExN03 application
  runMgr->SetUserApplication ( new ExN03ApplicationRP(runMgr) );
#ifdef GEANT_TBB
  if (tbbmode)
    runMgr->SetTaskMgr( new TaskMgrTBB() );
#endif
// print run information
    std::cout<< "\n\n"
             << " =============================================================================== \n"
             << "  primary GV code      : " << gvParticleCode << "  [22 => e-; 23 => e+; 42 => gamma]\n"
             << "  primary energy       : " << primaryEnergy  << " [GeV] \n"
             << "  primary position     : " << "[ " << xPos << ", " << yPos << ", " << zPos << " ] \n"
             << "  primary direction    : " << "[ " << xDir << ", " << yDir << ", " << zDir << " ] \n"
             << "  #events              : " << n_events << " \n"
             << "  #primaries per event : " << n_avrg_tracks_per_evt << " \n"
             << "  total # primaries    : " << n_events*n_avrg_tracks_per_evt << " \n"
             << " ===============================================================================\n\n";
  // run the simulation
  runMgr->RunSimulation();
  // delete the run manager
  delete runMgr;
  return 0;
}
#endif
