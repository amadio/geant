// The following in ROOT v6 equivalent to gSystem->Load("../lib/libGeant_v");
// R__LOAD_LIBRARY(libGeant_v)

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

// Autoload the library early so that GeantPropagator is defined when applicable.
class TaskBroker;
class GeantPropagator;

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
   TGeoManager::Import(geomfile);
   
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

   GeantPropagator *prop = GeantPropagator::NewInstance(ntotal, nbuffered, nthreads);
   prop->fBmag = magfield; // 4 Tesla

   //  Enable use of RK integration in field for charged particles
   prop->fUseRungeKutta = false;
   // prop->fEpsilonRK = 0.001;  // Revised / reduced accuracy - vs. 0.0003 default 

   if (broker) prop->SetTaskBroker(broker);

   // Monitor different features
   prop->SetNminThreshold(5*nthreads);
   prop->SetMonitored(GeantPropagator::kMonQueue,          true & (!performance));
   prop->SetMonitored(GeantPropagator::kMonMemory,         false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonBasketsPerVol,  false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonVectors,        false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonConcurrency,    false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonTracksPerEvent, false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonTracks,         false & (!performance));
   bool graphics = (prop->GetMonFeatures()) ? true : false;
   prop->fUseMonitoring = graphics;   

   // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
   // If set to 0 takes the default value of 0.01
   prop->fPriorityThr = 0.1;

   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
   prop->fNperBasket = 16;   // Initial vector size (tunable)

   // This is now the most important parameter for memory considerations
   prop->fMaxPerBasket = 64;   // Maximum vector size (tunable)

   // Minimum number of tracks in a basket that stay in the same volume and
   // therefore can be re-used with the same thread without re-basketizing.
   prop->fNminReuse = 10000;   // (default in propagator is 4)

   // Kill threshold - number of steps allowed before killing a track 
   //                  (includes internal geometry steps)
   prop->fNstepsKillThr = 100000;

   // Maximum user memory limit [MB]
   prop->fMaxRes = 4000;
   if (performance) prop->fMaxRes = 0;

   prop->fEmin = 0.001; // [1 MeV] energy cut

   prop->fEmax = 0.01; // 10 MeV
   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
//   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   //   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
//   prop->fPrimaryGenerator = new GunGenerator(1, 0, 1., 0, 0, 0, 0.362783697740757, 0.259450124768640, 0.882633622956438);
   std::string s(eventfile);
   prop->fPrimaryGenerator = new HepMCGenerator(s);
//   prop->fPrimaryGenerator->SetEtaRange(-2.4,2.4);
//   prop->fPrimaryGenerator->SetMomRange(0.,0.5);
   //   prop->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
   prop->fLearnSteps = 100000;
   if (performance) prop->fLearnSteps = 0;

   CMSApplication *app = new CMSApplication(prop);
   // Activate I/O
   prop->fFillTree = false;
   prop->fTreeSizeWriteThreshold = 100000;
   // Activate old version of single thread serialization/reading
//   prop->fConcurrentWrite = false;
   app->SetScoreType(CMSApplication::kScore);
//   if (performance) app->SetScoreType(CMSApplication::kNoScore);
   prop->fApplication = app;

//   gROOT->ProcessLine(".x factory.C+");   
// Activate debugging using -DBUG_HUNT=ON in your cmake build
   prop->fDebugEvt = 0;
   prop->fDebugTrk = 0;
   prop->fDebugStp = 0;
   prop->fDebugRep = 10;
   
// Activate standard scoring   
   prop->fUseStdScoring = false; // true;
   if (performance) prop->fUseStdScoring = false;
   prop->fUseMonitoring = graphics;
   prop->PropagatorGeom(nthreads);
   delete prop;
}
