#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

// Autoload the library early so that GeantPropagator is defined when applicable.
class TaskBroker;
class GeantPropagator;

void runFastSim( int ncputhreads=4,
                 bool performance=true,
	         const char *geomfile="Par02FullDetector.root",
                 bool coprocessor = COPROCESSOR_REQUEST
               )
{

//=============================================================================
// PERFORMANCE MODE SWITCH: no scoring, no memory cleanup thread, no monitoring
//=============================================================================
//   bool performance = true;

   int nthreads   = ncputhreads;
   int ntotal     = 10000;  // Number of events to be transported
   int nbuffered  = 10;  // Number of buffered events (tunable [1,ntotal])
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
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered, nthreads);

   if (broker) prop->SetTaskBroker(broker);

   // Monitor different features
   prop->SetNminThreshold(5*nthreads);
   prop->SetMonitored(GeantPropagator::kMonQueue,          true & (!performance));
   prop->SetMonitored(GeantPropagator::kMonMemory,         false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonBasketsPerVol,  false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonVectors,        false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonConcurrency,    false & (!performance));
   prop->SetMonitored(GeantPropagator::kMonTracksPerEvent, false & (!performance));
   bool graphics = (prop->GetMonFeatures()) ? true : false;
   prop->fUseMonitoring = graphics;
   prop->fNaverage = 10;   // Average number of tracks per event
  
   // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
   // If set to 0 takes the default value of 0.01
   prop->fPriorityThr = 0.05;

   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
   prop->fNperBasket = 16;   // Initial vector size (tunable)

   // This is now the most important parameter for memory considerations
   prop->fMaxPerBasket = 256;   // Maximum vector size (tunable)

   prop->fEmin = 1.E-3; // [1 MeV] energy cut
   prop->fEmax = 50.0;  // Particle gun energy in GeV

   prop->fBmag = 1.0;  // Magnetic field in Tesla units

   // Create the fast-sim process.
   prop->fProcess = new FastSimProcess;

   // for vector physics -OFF now
   // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin, nthreads);

   prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 11, prop->fEmax, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );  //e-
   //prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 22, prop->fEmax, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );  //gamma
   //prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 13, prop->fEmax, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );  //mu-
   //prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 13, prop->fEmax, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 );  //mu- along x
   //prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 13, prop->fEmax, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 );  //mu- along z
   //prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 211, prop->fEmax, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );  //pi+
   //prop->fPrimaryGenerator = new GunGenerator( prop->fNaverage, 111, prop->fEmax, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );  //pi0

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
   prop->fLearnSteps = 0;
   if (performance) prop->fLearnSteps = 0;
   
   prop->fApplication = new FastSimApplication;  // It does the analysis

   // Activate I/O
   prop->fFillTree = false;
   // Activate old version of single thread serialization/reading
   prop->fConcurrentWrite = false;

   // Activate debugging using -DBUG_HUNT=ON in your cmake build
   prop->fDebugEvt = 0;
   prop->fDebugTrk = 0;
   prop->fDebugStp = 0;
   prop->fDebugRep = 10;

// Activate standard scoring   
   prop->fUseStdScoring = false;
   if (performance) prop->fUseStdScoring = false;
   // Monitor the application
   prop->fUseAppMonitoring = false;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
   delete prop;
}
