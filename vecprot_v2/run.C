void run(Int_t nthreads=2, Bool_t graphics=kFALSE, 
//         const char *geomfile="gexam1.root")
//         const char *geomfile="http://root.cern.ch/files/cms.root")
	 const char *geomfile="Ex03.root")
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


   Int_t ntotal   = 100;  // Number of events to be transported
   Int_t nbuffered  = 10;   // Number of buffered events
   TGeoManager::Import(geomfile);


   //   gSystem->Load("../lib/libGeant_v_vecGeom.so");

   Int_t ntotal   = 1000;  // Number of events to be transported
   Int_t nbuffered  = 20;   // Number of buffered events

   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   prop->fNaverage = 10;   // Average number of tracks per event

   prop->fNperBasket = 8;   // Vector size
   prop->fEmin = 1.E-5; // [10KeV] energy cut
   prop->fEmax = 0.03.; // [30MeV] used for now to select particle gun energy

   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", "xsec_FTFP_BERT.root", "fstate_FTFP_BERT.root");

   prop->fApplication = new MyApplication();

//   gROOT->ProcessLine(".x factory.C+");   
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
