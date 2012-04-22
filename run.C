void run(Int_t nthreads=3, Bool_t graphics=kTRUE, const char *geomfile="http://root.cern.ch/files/cms.root")
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
   prop->fNevents  = 100;     // Number of events to be transported
   prop->fNaverage = 100;   // Average number of tracks per event
   prop->PropagatorGeom(geomfile, nthreads, graphics);
}   
   
   
