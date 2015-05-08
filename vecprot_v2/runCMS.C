void runCMS(Int_t nthreads=4,
	 const char *geomfile="../cmstrack/cms2015.root",
	 const char *xsec="xsec_FTFP_BERT_G496p02_1mev.root",
	 const char *fstate="fstate_FTFP_BERT_G496p02_1mev.root")
{
   gSystem->Load("libPhysics");
   gSystem->Load("libHist");
   gSystem->Load("libThread");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("../lib/libGeant_v");
   gSystem->Load("../lib/libXsec");
   gSystem->Load("../lib/libGeantExamples");

//=============================================================================
// PERFORMANCE MODE SWITCH: no scoring, no memory cleanup thread, no monitoring
//=============================================================================
   Bool_t performance = false;

   Int_t ntotal   = 10;  // Number of events to be transported
   Int_t nbuffered  = 5;   // Number of buffered events (tunable [1,ntotal])
   TGeoManager::Import(geomfile);
   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   // Monitor different features
   wmgr->SetNminThreshold(5*nthreads);
   wmgr->SetMonitored(WorkloadManager::kMonQueue,          true & (!performance));
   wmgr->SetMonitored(WorkloadManager::kMonMemory,         false & (!performance));
   wmgr->SetMonitored(WorkloadManager::kMonBasketsPerVol,  false & (!performance));
   wmgr->SetMonitored(WorkloadManager::kMonVectors,        false & (!performance));
   wmgr->SetMonitored(WorkloadManager::kMonConcurrency,    false & (!performance));
   wmgr->SetMonitored(WorkloadManager::kMonTracksPerEvent, false & (!performance));
   wmgr->SetMonitored(WorkloadManager::kMonTracks,         false & (!performance));
   Bool_t graphics = (wmgr->GetMonFeatures()) ? true : false;
   prop->fUseMonitoring = graphics;   

   // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
   // If set to 0 takes the default value of 0.01
   prop->fPriorityThr = 0.1;

   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
   prop->fNperBasket = 16;   // Initial vector size (tunable)

   // This is now the most important parameter for memory considerations
   prop->fMaxPerBasket = 64;   // Maximum vector size (tunable)

   // Maximum user memory limit [MB]
   prop->fMaxRes = 4000;
   if (performance) prop->fMaxRes = 0;

   prop->fEmin = 0.001; // [1 MeV] energy cut

   prop->fEmax = 0.01; // 10 MeV
   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
//   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   //   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   std::string s = "pp14TeVminbias.root";
//   prop->fPrimaryGenerator = new GunGenerator(1, 0, 1., 0, 0, 0, 0.362783697740757, 0.259450124768640, 0.882633622956438);
   prop->fPrimaryGenerator = new HepMCGenerator(s);
//   prop->fPrimaryGenerator->SetEtaRange(-2.4,2.4);
//   prop->fPrimaryGenerator->SetMomRange(0.,0.5);
   //   prop->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
   prop->fLearnSteps = 100000;
   if (performance) prop->fLearnSteps = 0;

   CMSApplication *app = new CMSApplication();
   app->SetScoreType(CMSApplication::kScore);
   if (performance) app->SetScoreType(CMSApplication::kNoScore);
   prop->fApplication = app;

//   gROOT->ProcessLine(".x factory.C+");   
// Activate debugging using -DBUG_HUNT=ON in your cmake build
   prop->fDebugEvt = 0;
   prop->fDebugTrk = 0;
   prop->fDebugStp = 0;
   prop->fDebugRep = 10;
   
// Activate standard scoring   
   prop->fUseStdScoring = true;
   if (performance) prop->fUseStdScoring = false;
   prop->fUseMonitoring = graphics;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
   delete prop;
}
