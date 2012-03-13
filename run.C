void run(Int_t nthreads=4, const char *geomfile="http://root.cern.ch/files/cms.root")
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeant.so");
   
   GeantPropagator *prop = GeantPropagator::Instance();
   prop->fNevents  = 10;     // Number of events to be transported
   prop->fNaverage = 1000;   // Average number of tracks per event
   prop->PropagatorGeom(geomfile, nthreads);
}   
   
   
