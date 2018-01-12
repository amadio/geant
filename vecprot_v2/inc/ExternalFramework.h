#ifndef USER_EXTERNAL_FRAMEWORK_H
#define USER_EXTERNAL_FRAMEWORK_H

#include <thread>
#include <atomic>
#include <vector>
#include "Geant/Config.h"
#include "GeantRunManager.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "PrimaryGenerator.h"
#include "EventSet.h"
#include "GeantEvent.h"
#include "GeantEventServer.h"
#include "GeantConfig.h"

namespace userfw {

class Framework
{
  using GeantRunManager = Geant::GeantRunManager;
  using GeantVApplication = Geant::GeantVApplication;
  using PrimaryGenerator = Geant::PrimaryGenerator;
  using GeantConfig = Geant::GeantConfig;
  using GeantTaskData = Geant::GeantTaskData;
  using EventSet = Geant::EventSet;
  using GeantEvent = Geant::GeantEvent;
  using GeantEventInfo = Geant::GeantEventInfo;
  using GeantTrack = Geant::GeantTrack;

private:
  bool fInitialized = false;                 /** Initialization flag */
  size_t fNthreads = 1;                      /** Number of threads */
  size_t fNevents = 0;                       /** Number of events to be transported */
  std::atomic<size_t> fIevent;               /** Current generated event */
  GeantRunManager *fGeantRunMgr = nullptr;   /** Geant run manager */
  PrimaryGenerator *fGenerator = nullptr;    /** Generator used in external loop mode */
  
public:
  Framework(size_t nthreads, size_t nevents, GeantRunManager *mgr, PrimaryGenerator *gen)
    :fNthreads(nthreads), fNevents(nevents), fGeantRunMgr(mgr), fGenerator(gen) {}
  Framework(const Framework &) = delete;

  ~Framework() {
    delete fGeantRunMgr;
    if (fGenerator != fGeantRunMgr->GetPrimaryGenerator())
      delete fGenerator;
  }
  
  Framework &operator=(const Framework &) = delete;

  GeantRunManager *GetRunMgr() const { return fGeantRunMgr; }

  void Initialize() {
    if (fInitialized) return;
    fIevent.store(0);
    fGeantRunMgr->GetConfig()->fRunMode = GeantConfig::kExternalLoop;
    if (!fGeantRunMgr->IsInitialized()) fGeantRunMgr->Initialize();
    fGenerator->InitPrimaryGenerator();
    fInitialized = true;
  }
  
  EventSet *GenerateEventSet(size_t nevents, GeantTaskData *td)
  {  
    GeantEvent **events = new GeantEvent*[nevents];
    size_t nstored = 0;
    for (size_t i=0 ; i< nevents; ++i) {
      // Book an event number
      const size_t evt = fIevent.fetch_add(1);
      if (evt >= fNevents) break;
      GeantEvent *event = new GeantEvent();
      GeantEventInfo event_info = fGenerator->NextEvent(td);
      while (event_info.ntracks == 0) {
        printf("Discarding empty event\n");
        event_info = fGenerator->NextEvent(td);
      }
      event->SetEvent(evt);
      event->SetNprimaries(event_info.ntracks);
      event->SetVertex(event_info.xvert, event_info.yvert, event_info.zvert);
      for (int itr = 0; itr < event_info.ntracks; ++itr) {
        GeantTrack &track = td->GetNewTrack();
        track.SetParticle(event->AddPrimary(&track));
        track.SetPrimaryParticleIndex(itr);
        fGenerator->GetTrack(itr, track, td);
      }
      events[nstored++] = event;
    }
    if (nstored == 0) return nullptr;
    EventSet *evset = new EventSet(nstored);
    for (size_t i=0 ; i< nstored; ++i)
      evset->AddEvent(events[i]);
    delete [] events;
    return evset;
  }


  static
  bool RunTransportTask(size_t ntotransport, Framework *fw)
  {
    // This is the entry point for the user code to transport as a task a set of
    // events.
  
    GeantRunManager *runMgr = fw->GetRunMgr();
    // First book a transport task from GeantV run manager
    while (1) {
      GeantTaskData *td = runMgr->BookTransportTask();
      if (!td) return false;
  
      // ... then create the event set
      EventSet *evset = fw->GenerateEventSet(ntotransport, td);
      if (!evset) return true;  // work completed
  
      // ... finally invoke the GeantV transport task
      bool transported = runMgr->RunSimulationTask(evset, td);
      // Now we could run some post-transport task
      
      if (!transported) return false;
    }
    return true;
  }

  void Run()
  {
    // Run the external event loop
    Initialize();
    size_t ntotransport = fNevents/fNthreads;
    if (ntotransport < 1) ntotransport = 1;
    printf("=== RUNNING SIMULATION WITH EXTERNAL LOOP\n");
    std::vector<std::thread> workers;
    for (size_t n = 0; n < fNthreads; ++n) {
      workers.emplace_back(RunTransportTask, ntotransport, this);
    }
    for (auto& worker : workers) {
      worker.join();
    }
    printf("=== SIMULATION WITH EXTERNAL LOOP COMPLETED\n");
  }

};

} // userfw

#endif // USER_EXTERNAL_FRAMEWORK_H
