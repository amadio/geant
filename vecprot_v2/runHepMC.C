void runHepMC(Int_t nthreads=4,
         Bool_t graphics=kFALSE,
//         const char *geomfile="simple_ecal.root")
//         const char *geomfile="http://root.cern.ch/files/cms.root")
	 const char *geomfile="ExN03.root",
	 const char *xsec="xsec_FTFP_BERT.root",
	 const char *fstate="fstate_FTFP_BERT.root")
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

   Int_t ntotal   = 20;  // Number of events to be transported
   Int_t nbuffered  = 10;   // Number of buffered events
   TGeoManager::Import(geomfile);
   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   prop->fNaverage = 400;   // Average number of tracks per event
   prop->fNperBasket = 128;   // Vector size
   prop->fEmin = 3.E-6; // [3 KeV] energy cut
//   prop->fEmax = 0.03.; // [30MeV] used for now to select particle gun energy
   prop->fEmax = 0.03;
   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);

   //   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   std::string fname = "pp14TeVminbias.root";
   prop->fPrimaryGenerator = new HepMCGenerator(fname);
   //   prop->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");

   prop->fApplication = new MyApplication();

//   gROOT->ProcessLine(".x factory.C+");   
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   prop->fUseMonitoring = graphics;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
