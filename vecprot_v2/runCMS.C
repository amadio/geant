// The following in ROOT v6 equivalent to gSystem->Load("../lib/libGeant_v");
// R__LOAD_LIBRARY(libGeant_v)

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

// Autoload the library early so that GeantPropagator is defined when applicable.
namespace Geant {
inline namespace cxx {
class TaskBroker;
class GeantPropagator;
}
}

void runCMS(const int ncputhreads=4,
            const bool performance=true,
            const char *geomfile="../cmstrack/cms2015.root",
            const char *xsec="xsec_FTFP_BERT_G496p02_1mev.root",
            const char *fstate="fstate_FTFP_BERT_G496p02_1mev.root",
            bool coprocessor = COPROCESSOR_REQUEST,
	    const char *eventfile="pp14TeVminbias.root",
	    const float magfield=40.,
	    const int ntotal=10                                    // Number of events to be transported
)
{
   // gSystem->Load("libPhysics");
   // gSystem->Load("libHist");
   // gSystem->Load("libThread");
   // gSystem->Load("libGeom");
   // gSystem->Load("libVMC");
   // gSystem->Load("../lib/libGeant_v");
   // gSystem->Load("../lib/libXsec");
   // gSystem->Load("../lib/libGeantExamples");

   printf("\n===================================== Input parameters ========================================\n");
   printf("Number of threads:                 %d\n",ncputhreads);
   printf("Number of event requested:         %d\n",ntotal);
   printf("Geometry file:                     %s\n",geomfile);
   printf("Cross sections:                    %s\n",xsec);
   printf("Final state:                       %s\n",fstate);
   printf("Input event file:                  %s\n",eventfile);
   printf("Mag field (kGauss)                 %f\n",magfield);
   printf("===================================== Input parameters ========================================\n\n");


//=============================================================================
// PERFORMANCE MODE SWITCH: no scoring, no memory cleanup thread, no monitoring
//=============================================================================
//   bool performance = true;

   int nthreads = ncputhreads;
   int nbuffered  = 5;   // Number of buffered events (tunable [1,ntotal])
   int npropagators = 1;

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

   Geant::GeantConfig* config = new Geant::GeantConfig();
   config->fGeomFileName = geomfile;
   config->fNtotal = ntotal;
   config->fNbuff = nbuffered;
   config->fBmag = magfield; // 4 Tesla

   //  Enable use of RK integration in field for charged particles
   config->fUseRungeKutta = false;
   // config->fEpsilonRK = 0.001;  // Revised / reduced accuracy - vs. 0.0003 default

   GeantRunManager *runMgr = new GeantRunManager(npropagators, nthreads, config);
   if (broker) runMgr->SetCoprocessorBroker(broker);

   // Monitor different features
   config->fNminThreshold = 5*nthreads;
   config->SetMonitored(GeantConfig::kMonQueue,          true & (!performance));
   config->SetMonitored(GeantConfig::kMonMemory,         false & (!performance));
   config->SetMonitored(GeantConfig::kMonBasketsPerVol,  false & (!performance));
   config->SetMonitored(GeantConfig::kMonVectors,        false & (!performance));
   config->SetMonitored(GeantConfig::kMonConcurrency,    false & (!performance));
   config->SetMonitored(GeantConfig::kMonTracksPerEvent, false & (!performance));
   config->SetMonitored(GeantConfig::kMonTracks,         false & (!performance));
   bool graphics = (config->GetMonFeatures()) ? true : false;
   config->fUseMonitoring = graphics;

   // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
   // If set to 0 takes the default value of 0.01
   config->fPriorityThr = 0.1;

   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
   config->fNperBasket = 16;   // Initial vector size (tunable)

   // This is now the most important parameter for memory considerations
   config->fMaxPerBasket = 64;   // Maximum vector size (tunable)

   // Minimum number of tracks in a basket that stay in the same volume and
   // therefore can be re-used with the same thread without re-basketizing.
   config->fNminReuse = 10000;   // (default in confi is 4)

   // Kill threshold - number of steps allowed before killing a track 
   //                  (includes internal geometry steps)
   config->fNstepsKillThr = 100000;

   // Maximum user memory limit [MB]
   config->fMaxRes = 4000;
   if (performance) config->fMaxRes = 0;

   config->fEmin = 0.001; // [1 MeV] energy cut

   config->fEmax = 0.01; // 10 MeV
   // Create the tab. phys process.
   runMgr->SetPhysicsProcess( new Geant::TTabPhysProcess("tab_phys", xsec, fstate));
//   config->fPrimaryGenerator = new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0);
   //   config->fPrimaryGenerator = new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0);
//   config->fPrimaryGenerator = new GunGenerator(1, 0, 1., 0, 0, 0, 0.362783697740757, 0.259450124768640, 0.882633622956438);
   std::string s(eventfile);
   runMgr->SetPrimaryGenerator(new Geant::HepMCGenerator(s));
//   config->fPrimaryGenerator->SetEtaRange(-2.4,2.4);
//   config->fPrimaryGenerator->SetMomRange(0.,0.5);
   //   config->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
   config->fLearnSteps = 100000;
   if (performance) config->fLearnSteps = 0;

   CMSApplication *app = new CMSApplication(runMgr);
   // Activate I/O
   config->fFillTree = false;
   config->fTreeSizeWriteThreshold = 100000;
   // Activate old version of single thread serialization/reading
//   config->fConcurrentWrite = false;
   app->SetScoreType(CMSApplication::kScore);
//   if (performance) app->SetScoreType(CMSApplication::kNoScore);
   runMgr->SetUserApplication( app );

//   gROOT->ProcessLine(".x factory.C+");   
// Activate debugging using -DBUG_HUNT=ON in your cmake build
   config->fDebugEvt = 0;
   config->fDebugTrk = 0;
   config->fDebugStp = 0;
   config->fDebugRep = 10;
   
// Activate standard scoring   
   config->fUseStdScoring = false; // true;
   if (performance) config->fUseStdScoring = false;
   config->fUseMonitoring = graphics;
   runMgr->RunSimulation();
   // config->PropagatorGeom(nthreads);
   delete config;
}
