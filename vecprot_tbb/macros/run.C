void run(Int_t nthreads=16, Bool_t graphics=kTRUE, const char *geomfile="../../geometry/cms.root",
         Int_t evtot=500,
         Int_t evbuf=100,
         Double_t trav=400.,
         Int_t trperbask=10,
         Int_t min_feeder_arg=50,
         Int_t nprior=5,
         Int_t dispthr=100)
{
   gSystem->Load("libtbb2.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeant.so");

   GeantMainPropagator *mainprop = GeantMainPropagator::Instance();

   // Note that min_feeder is currently hard-coded to
   // TMath::Max(min_feeder_arg, 2*nthreads)

   // n_threads, events_total, events_buffered, tracks_average, max_per_basket
   // min_feeder, n_events_to_prioritize, threshold_to_start_DispTask
   mainprop->SetParams (nthreads, evtot, evbuf, trav, trperbask, min_feeder_arg, nprior, dispthr);

   mainprop->Start(geomfile, graphics);
}

