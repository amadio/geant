void run(Int_t nthreads=1, Bool_t graphics=kFALSE, 
         const char *geomfile="ExN03.root")
//         const char *geomfile="http://root.cern.ch/files/cms.root")
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeant_v.so");
   gSystem->Load("libXsec.so");

   Int_t ntotal   = 20;  // Number of events to be transported
   Int_t nbuffered  = 10;   // Number of buffered events
   TGeoManager::Import(geomfile);
   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   prop->fNaverage = 100;   // Average number of tracks per event
   prop->fNperBasket = 8;
   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", "xsec_FTFP_BERT.root", "fstate_FTFP_BERT.root");

   prop->fApplication = new MyApplication();

//   gROOT->ProcessLine(".x factory.C+");   
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
