#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

#include <err.h>
#include <getopt.h>
#include <unistd.h>
#ifdef USE_ROOT
#include "Rtypes.h"
#include "TGeoManager.h"
#endif
#include "GeantRunManager.h"
#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TTabPhysProcess.h"
#include "CMSApplication.h"
#ifdef GEANT_TBB
#include "TaskMgrTBB.h"
#endif
#include "CMSFieldConstruction.h"

using namespace Geant;

static int n_events = 10;
static int n_buffered = 16;
static int n_threads = 4;
static int n_track_max = 64;
static int n_learn_steps = 100000;
static int n_reuse = 100000;
static int n_propagators = 1;
static int max_memory = 4000; /* MB */
static bool monitor = false, score = false, debug = false, coprocessor = false;
static bool tbbmode = false, usev3 = true, usenuma = false;

static struct option options[] = {{"events", required_argument, 0, 'e'},
                                  {"hepmc-event-file", required_argument, 0, 'E'},
                                  {"fstate", required_argument, 0, 'f'},
                                  {"field-file", required_argument, 0, 'F'},
                                  {"geometry", required_argument, 0, 'g'},
                                  {"learn-steps", required_argument, 0, 'l'},
                                  {"max-tracks-per-basket", required_argument, 0, 'B'},
                                  {"use-Runge-Kutta", no_argument, 0, 'K'},
                                  {"use-CMS-field", no_argument, 0, 'C'},
                                  {"use-Uniform-field", no_argument, 0, 'M'},
                                  {"monitor", no_argument, 0, 'm'},
                                  {"debug", no_argument, 0, 'd'},
                                  {"max-memory", required_argument, 0, 'M'},
                                  {"nbuffer", required_argument, 0, 'b'},
                                  {"score", no_argument, 0, 's'},
                                  {"threads", required_argument, 0, 't'},
                                  {"xsec", required_argument, 0, 'x'},
                                  {"coprocessor", required_argument, 0, 'r'},
                                  {"tbbmode", required_argument, 0, 'i'},
                                  {"reuse", required_argument, 0, 'u'},
                                  {"propagators", required_argument, 0, 'p'},
                                  {"v3", required_argument, 0, 'v'},
                                  {"numa", required_argument, 0, 'n'},
                                  {0, 0, 0, 0}};

void help() {
  printf("\nUsage: cmsapp [OPTIONS] INPUT_FILE\n\n");

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

int main(int argc, char *argv[]) {
  std::cout << "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" << std::endl;
  std::string cms_geometry_filename("cms2015.root");
  std::string xsec_filename("xsec_FTFP_BERT_G496p02_1mev.root");
  std::string fstate_filename("fstate_FTFP_BERT_G496p02_1mev.root");
  std::string hepmc_event_filename("pp14TeVminbias.root");
  std::string field_filename("CMSmagneticField.txt");
  bool useRungeKutta= false;
  bool useCMSfield=   false;   //  If false, use a constant magnetic field

  if (argc == 1) {
    help();
    exit(0);
  }

  while (true) {
    int c, optidx = 0;

    // c = getopt_long(argc, argv, "E:e:f:g:l:B:mM:b:t:x:r:i:u:p:v:n:", options, &optidx);
    //   = getopt_long(argc, argv, "cE:e:f:F:g:l:B:mM:b:t:x:r:Ku", options, &optidx);    
    c = getopt_long(argc, argv, "cE:e:f:F:g:l:B:mM:b:t:x:r:i:Ku:p:v:n:", options, &optidx);
    
    if (c == -1)
      break;

    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */

    case 'b':
      n_buffered = (int)strtol(optarg, NULL, 10);

      if (n_buffered < 1)
        errx(1, "number of buffered events must be positive");
      break;

    case 'B':
      n_track_max = (int)strtol(optarg, NULL, 10);

      if (n_track_max < 1)
        errx(1, "max number of tracks per basket must be positive");
      break;

    case 'c':
    case 'C':
      useCMSfield = true;
      break;

    case 'U':  // Uniform field
      useCMSfield = false;
      break;

    case 'e':
      n_events = (int)strtol(optarg, NULL, 10);

      if (n_events <= 0)
        errx(1, "number of events must be positive");
      break;

    case 'E':
      hepmc_event_filename = optarg;
      break;

    case 'f':
      fstate_filename = optarg;
      break;

    case 'F':
      field_filename = optarg;
      break;

    case 'g':
      cms_geometry_filename = optarg;
      break;

    case 'K':
      useRungeKutta = true;
      break;

    case 'l':
      n_learn_steps = (int)strtol(optarg, NULL, 10);

      if (n_learn_steps <= 0)
        errx(1, "number of learning steps must be positive");
      break;

    case 'm':
      monitor = true;
      break;

    case 'M':
      max_memory = (int)strtol(optarg, NULL, 10);

      if (max_memory < 128)
        errx(1, "max memory is too low");
      break;

    case 'r':
      coprocessor = optarg;
      break;

    case 's':
      score = true;
      break;

    case 't':
      n_threads = (int)strtol(optarg, NULL, 10);

      if (n_threads < 1)
        errx(1, "number of threads must be positive");

      break;

    case 'x':
      xsec_filename = optarg;
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
      usev3 = bool(strtol(optarg, NULL, 10));
      if (!usev3) n_buffered = 4;
      break;

    case 'n':
      usenuma = bool(strtol(optarg, NULL, 10));
      break;

    default:
      errx(1, "unknown option %c", c);
    }
  }

  bool performance = true;   

  TaskBroker *broker = nullptr;
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

  GeantConfig* config=new GeantConfig();
  
  config->fGeomFileName = cms_geometry_filename;
  config->fNtotal = n_events;
  config->fNbuff = n_buffered;
  // Default value is 1. (0.1 Tesla)
  // config->fBmag = 40.; // 4 Tesla

  // V3 options
  config->fNmaxBuffSpill = 8;  // New configuration parameter!!!
  config->fNstackLanes = 10;
  config->fUseV3 = usev3;
  config->fUseNuma = usenuma;

  // Enable use of RK integration in field for charged particles
  config->fUseRungeKutta = false;
  // prop->fEpsilonRK = 0.001;  // Revised / reduced accuracy - vs. 0.0003 default

  config->fNminThreshold=5 * n_threads;
  config->fUseMonitoring = monitor;
  config->fNaverage = 500;

  config->SetMonitored(GeantConfig::kMonQueue, false);
  config->SetMonitored(GeantConfig::kMonMemory, monitor);
  config->SetMonitored(GeantConfig::kMonBasketsPerVol, false);
  config->SetMonitored(GeantConfig::kMonVectors, false);
  config->SetMonitored(GeantConfig::kMonConcurrency, false);
  config->SetMonitored(GeantConfig::kMonTracksPerEvent, false);
  config->SetMonitored(GeantConfig::kMonTracks, false);
  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  config->fPriorityThr = 0.1;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  config->fNperBasket = 16; // Initial vector size

  // This is now the most important parameter for memory considerations
  config->fMaxPerBasket = n_track_max;

  // Maximum user memory limit [MB]
  config->fMaxRes = max_memory;
  if (config) config->fMaxRes = 0;
  // CRITICAL: the energy cut must correspond to xsec/final state files !  
  config->fEmin = 0.001; // [1 MeV] energy cut
  config->fEmax = 0.01;  // 10 MeV
  printf("CMSApp> Production threshold set to %7.2g GeV\n", config->fEmin);  
  if (debug) {
    config->fUseDebug = true;
    config->fDebugTrk = 1;
    //propagator->fDebugEvt = 0;
    //propagator->fDebugStp = 0;
    //propagator->fDebugRep = 10;
  }
  config->fUseMonitoring = monitor;

  // Set threshold for tracks to be reused in the same volume
  config->fNminReuse = n_reuse;

  // Activate standard scoring   
  config->fUseStdScoring = true;
  if (performance) config->fUseStdScoring = false;
  config->fLearnSteps = n_learn_steps;
  if (performance) config->fLearnSteps = 0;

  // Activate I/O
  config->fFillTree = false;
  config->fTreeSizeWriteThreshold = 100000;
  // Activate old version of single thread serialization/reading
  //   config->fConcurrentWrite = false;

  // Create run manager
  GeantRunManager *runMgr = new GeantRunManager(n_propagators, n_threads, config);
  if (broker) runMgr->SetCoprocessorBroker(broker);
  // Create the tab. phys process.
  runMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str()));

  // CMS magnetic field
  // propagator->fBmag = 40.; // 4 Tesla

  //  Enable use of RK integration in field for charged particles
  // config->fUseRungeKutta = false;
  config->fUseRungeKutta = useRungeKutta;

  config->fEpsilonRK = 0.0003;  // Revised / reduced accuracy - vs. 0.0003 default 

  if( useCMSfield ) {
     CMSFieldConstruction* CmsFieldCtr= new CMSFieldConstruction();
     CmsFieldCtr->SetFileForField(field_filename);
     printf("CMSApp: Setting CMS-detector-construction to GeantPropagator\n");
     runMgr // ->GetDetectorConstruction()
        ->SetUserFieldConstruction(CmsFieldCtr);
     // printf("Calling CreateFieldAndSolver from runCMS_new.C");
     // CmsFieldCtr->CreateFieldAndSolver(propagator->fUseRungeKutta);
  } else {
     UserFieldConstruction* fieldCtr= new UserFieldConstruction();
     float fieldVec[3] = { 0.0f, 0.0f, 38.0f };
     fieldCtr->UseConstantMagField( fieldVec, "kilogauss" );
     printf("CMSApp: Setting generic detector-construction to GeantPropagator - created field= %f %f %f.\n",
            fieldVec[0], fieldVec[1], fieldVec[2] );
     runMgr // ->GetDetectorConstruction()
        ->SetUserFieldConstruction(fieldCtr);
  }

#ifdef USE_VECGEOM_NAVIGATOR
#ifdef USE_ROOT
//  runMgr->LoadVecGeomGeometry();
#else
//  runMgr->LoadGeometry(cms_geometry_filename.c_str());
#endif
#endif

  if (hepmc_event_filename.empty()) {
    runMgr->SetPrimaryGenerator( new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0) );
  } else {
    // propagator->fPrimaryGenerator->SetEtaRange(-2.,2.);
    // propagator->fPrimaryGenerator->SetMomRange(0.,0.5);
    // propagator->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");
    runMgr->SetPrimaryGenerator( new HepMCGenerator(hepmc_event_filename) );
  }

  CMSApplication *CMSApp = new CMSApplication(runMgr);
  runMgr->SetUserApplication( CMSApp );
  if (score) {
    CMSApp->SetScoreType(CMSApplication::kScore);
  } else {
    CMSApp->SetScoreType(CMSApplication::kNoScore);
  }
#ifdef GEANT_TBB
  if (tbbmode)
    runMgr->SetTaskMgr( new TaskMgrTBB() );
#endif

  runMgr->RunSimulation();
//  propagator->PropagatorGeom(cms_geometry_filename.c_str(), n_threads, monitor);
  return 0;
}
