void run(Int_t nthreads=10, Bool_t graphics=kTRUE, const char *geomfile="http://root.cern.ch/files/cms.root")
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeant.so");
   
   GeantPropagator *prop = GeantPropagator::Instance();
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   prop->fNtotal   = 1500;  // Number of events to be transported
   prop->fNevents  = 100;   // Number of buffered events
   prop->fNaverage = 500;   // Average number of tracks per event
   prop->fNperBasket = 10;
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
