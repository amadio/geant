void run(Int_t nthreads=1, Bool_t graphics=kFALSE, 
         const char *geomfile="gexam1.root")
//         "http://root.cern.ch/files/cms.root")
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeant.so");
   gSystem->Load("libUser.so");   

   Int_t ntotal   = 20;  // Number of events to be transported
   Int_t nbuffered  = 10;   // Number of buffered events
   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   prop->fNaverage = 500;   // Average number of tracks per event
   prop->fNperBasket = 8;

   prop->fApplication = new MyApplication();

//   gROOT->ProcessLine(".x factory.C+");   
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
