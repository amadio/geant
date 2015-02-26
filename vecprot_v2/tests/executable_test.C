/*
 * executable_test.C
 *
 */
#include "TThread.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "WorkloadManager.h"
#include "PhysicsProcess.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"
#include "TObject.h"
#include "GeantVApplication.h"
#include "MyApplication.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "TROOT.h"
#include <iostream>
#include "GunGenerator.h"
#include "TTabPhysProcess.h"

void run(Int_t nthreads=4,
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
   // gSystem->Load("../buildTGeo/lib/libGeant_v");
   // gSystem->Load("../buildTGeo/lib/libXsec");
   gSystem->Load("../lib/libGeant_v");
   gSystem->Load("../lib/libXsec");
   // for vector physics - OFF now
   // gSystem->Load("../lib/libVphysproc");

   Int_t ntotal   = 20;  // Number of events to be transported
   Int_t nbuffered  = 10;   // Number of buffered events
   TGeoManager::Import(geomfile);
   
   GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   prop->fNaverage = 500;   // Average number of tracks per event
   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
   prop->fNperBasket = 16;   // Initial vector size
   // This is now the most important parameter for memory considerations
   prop->fMaxPerBasket = 256;   // Initial vector size
   prop->fEmin = 3.E-6; // [3 KeV] energy cut
   // prop->fEmax = 0.03.; // [30MeV] used for now to select particle gun energy
   prop->fEmax = 0.03;
   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
   // for vector physics -OFF now
   // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin, nthreads);
   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   prop->fApplication = new MyApplication();
   // gROOT->ProcessLine(".x factory.C+");   
   // prop->fUseDebug = kTRUE;
   // prop->fDebugTrk = 1;
   prop->fUseMonitoring = graphics;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
   delete prop;
}   


int main(int argc, char* argv[]) {
   if (argc < 5) { 
      std::cerr << "Wrond usage - please add all possible values & paths to files";
      return 1;
   }
   run(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5]);
   return 0;
}
