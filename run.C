void run(const char *geomfile="geometry.root", Int_t nthreads=4)
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeant.so");
   
   GeantPropagator *prop = GeantPropagator::Instance();
   prop->fNevents  = 1;     // Number of events to be transported
   prop->fNaverage = 100;   // Average number of tracks per event
   prop->PropagatorGeom(geomfile, nthreads);
}   
   
   
