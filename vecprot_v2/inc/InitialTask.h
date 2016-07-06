#ifndef INITIALTASK
#define INITIALTASK

#include "Geant/Error.h"

#include "WorkloadManager.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantTaskData.h"
#include "PhysicsProcess.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"
#include "GeantVApplication.h"
#include "GeantFactoryStore.h"
#include "MyHit.h"
#include "TThread.h"
#include "TThreadMergingServer.h"
#include "TThreadMergingFile.h"
#include "FlowControllerTask.h"

#if USE_VECGEOM_NAVIGATOR
#include "base/TLS.h"
#include "management/GeoManager.h"
#include "materials/Medium.h"
#else
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#endif
#include "TaskBroker.h"

#include "tbb/task.h"

class InitialTask : public tbb::task
{
private:

public:
  InitialTask ();
  ~InitialTask ();

  tbb::task* execute ();

};

#endif //INITIALTASK
