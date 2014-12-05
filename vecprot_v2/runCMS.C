void runCMS(Int_t nthreads=4,
         Bool_t graphics=kTRUE,
//         const char *geomfile="simple_ecal.root")
//         const char *geomfile="http://root.cern.ch/files/cms.root")
	 const char *geomfile="../cmstrack/cms2015.root",
	 const char *xsec="xsec_FTFP_BERT_G496p02.root",
	 const char *fstate="fstate_FTFP_BERT_G496p02.root")
{
   gSystem->Load("libPhysics");
   gSystem->Load("libHist");
   gSystem->Load("libThread");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   //   gSystem->Load("../buildTGeo/lib/libGeant_v");
   // gSystem->Load("../buildTGeo/lib/libXsec");
   gSystem->Load("../lib/libGeant_v");
   gSystem->Load("../lib/libXsec");

   Int_t ntotal   = 5;  // Number of events to be transported
   Int_t nbuffered  = 2;   // Number of buffered events
   TGeoManager::Import(geomfile);
   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
//   prop->fNaverage = 400;   // Average number of tracks per event
   
   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
   prop->fNperBasket = 16;   // Initial vector size

   // This is now the most important parameter for memory considerations
   prop->fMaxPerBasket = 64;   // Maximum vector size

   // Maximum user memory limit [MB]
   prop->fMaxRes = 4000;

//   prop->fEmin = 3.E-6; // [3 KeV] energy cut
   prop->fEmin = 1.E-3; // [1 MeV] energy cut
//   prop->fEmax = 0.03.; // [30MeV] used for now to select particle gun energy
   prop->fEmax = 0.05;
   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
//   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   //   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   std::string s = "pp14TeVminbias.root";
   prop->fPrimaryGenerator = new HepMCGenerator(s);
   //   prop->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");

   prop->fApplication = new MyApplication();

//   gROOT->ProcessLine(".x factory.C+");   
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   prop->fUseMonitoring = graphics;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
