//
//  GeantVProducer.cpp
//  DispatchProcessingDemo
//
//  Created by Chris Jones on 2016/04/25.
//

#include <vector>
#include <memory>
#include <atomic>
#include "ConfiguredProducer.h"
#include "Event.h"
#include "Waiter.h"

#include "tbb/task.h"

#include <base/Vector3D.h>
#include "Geant/Config.h"
#include "GeantRunManager.h"
#include "EventSet.h"

//#include "HepMCGenerator.h"
#include "PhysicsProcessHandler.h"
#include "PhysicsListManager.h"

#include "CMSApplicationTBB.h"
#include "CMSPhysicsList.h"
#include "CMSDetectorConstruction.h"
#include "CMSParticleGun.h"

using namespace Geant;
using namespace cmsapp;

namespace {

  /* Stand in for work done by GeantV via TBB based tasks*/

  class SimulateTask : public tbb::task {
    tbb::task* m_waitingTask;

  public:
    SimulateTask(tbb::task* iWaiting) :
      m_waitingTask(iWaiting) {}

    tbb::task* execute() override {
      printf("SimTask: iWaiting=%p\n", m_waitingTask);

      //IMPORTANT: decrementing the reference count is what
      // tells the framework the simulation is done with
      // this event.
      m_waitingTask->decrement_ref_count();
      return nullptr;
    }

  };
}

namespace demo {

  class GeantVProducer : public ConfiguredProducer {
  public:
    GeantVProducer(const boost::property_tree::ptree& iConfig);
  private:
    virtual void produce(edm::Event&) override;
    //void runSimulation(tbb::task*);

    /** Functions using new GeantV interface */
    bool RunTransportTask(size_t nevents, const size_t *event_index);
    int BookEvents(int nrequested);

    /** @brief Generate an event set to be processed by a single task.
	Not required as application functionality, the event reading or generation
	can in the external event loop.
    */
    Geant::EventSet* GenerateEventSet(size_t nevents, const size_t *event_index, GeantTaskData *td);

    inline void LockEventGenerator() { while (fEventGeneratorLock.test_and_set(std::memory_order_acquire)) {}; }
    inline void UnlockEventGenerator() { fEventGeneratorLock.clear(std::memory_order_release); }

    int fNevents;
    std::vector<const Getter*> m_getters;
    GeantConfig* fConfig = nullptr;
    GeantRunManager* fRunMgr = nullptr;
    PrimaryGenerator* fPrimaryGenerator = nullptr;
    std::atomic_flag fEventGeneratorLock;              /** Spinlock for event set access locking */
  };

  GeantVProducer::GeantVProducer(const boost::property_tree::ptree& iConfig)
    : ConfiguredProducer(iConfig,kThreadSafeBetweenInstances)
    , fNevents(iConfig.get<int>("Nevents"))
  {
    fEventGeneratorLock.clear();
    registerProduct(demo::DataKey());
    
    for(const boost::property_tree::ptree::value_type& v: iConfig.get_child("toGet")) {
      m_getters.push_back(registerGet(v.second.get<std::string>("label"), 
                                      ""));
    }

    int n_threads = iConfig.get<int>("Nthreads");
    int n_propagators = 1;
    int n_track_max = 500;
    int n_reuse = 100000;
    bool usev3 = true, usenuma = false;
    bool performance = true;

    //e.g. cms2015.root, cms2018.gdml, ExN03.root
    std::string cms_geometry_filename = iConfig.get<std::string>("geometry");

    std::string hepmc_event_filename = iConfig.get<std::string>("hepmc");
    // pp14TeVminbias.root sequence #stable: 608 962 569 499 476 497 429 486 465 619
    // minbias_14TeV.root sequence #stable: 81 84 93 97 87 60 106 91 92 60
    // if empty, use gun generator!

    // instantiate configuration helper
    fConfig = new GeantConfig();

    fConfig->fRunMode = GeantConfig::kExternalLoop;

    fConfig->fGeomFileName = cms_geometry_filename;
    fConfig->fNtotal = fNevents;
    fConfig->fNbuff = n_threads; // iConfig.get<int>("Nbuffered");

    // V3 options
    fConfig->fNstackLanes = 10;
    fConfig->fNmaxBuffSpill = 128;  // New configuration parameter!!!
    fConfig->fUseV3 = usev3;

    fConfig->fUseRungeKutta = false;  // Enable use of RK integration in field for charged particles
    // prop->fEpsilonRK = 0.001;      // Revised / reduced accuracy - vs. 0.0003 default

    fConfig->fUseNuma = usenuma;
    fConfig->fNminThreshold = 5 * n_threads;
    fConfig->fNaverage = 50;

    // Initial vector size, this is no longer an important model parameter,
    // because is gets dynamically modified to accomodate the track flow
    fConfig->fNperBasket = 16; // Initial vector size (tunable)

    // This is now the most important parameter for memory considerations
    fConfig->fMaxPerBasket = n_track_max;  // Maximum vector size (tunable)

    // This is temporarily used as gun energy
    fConfig->fEmax = 10;   // [GeV] used for now to select particle gun energy

    // Activate I/O
    fConfig->fFillTree = false;

    // Set threshold for tracks to be reused in the same volume
    fConfig->fNminReuse = n_reuse;

    // Activate standard scoring   
    fConfig->fUseStdScoring = false;
    if (performance) fConfig->fUseStdScoring = false;

    // Activate vectorized geometry (for now does not work properly with MT)
    fConfig->fUseVectorizedGeom = false;

     // Create run manager
    std::cerr<<"*** GeantRunManager: instantiating with "<< n_propagators <<" propagators and "<< n_threads <<" threads.\n";
    fRunMgr = new GeantRunManager(n_propagators, n_threads, fConfig);

    // Create the tabulated physics process
    std::cerr<<"*** GeantRunManager: setting physics process...\n";
    //fRunMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str()));

    // create the real physics main manager/interface object and set it in the GeantRunManager
    fRunMgr->SetPhysicsInterface(new geantphysics::PhysicsProcessHandler());

    // Create user defined physics list for the full CMS application
    geantphysics::PhysicsListManager::Instance().RegisterPhysicsList(new cmsapp::CMSPhysicsList());

    // Detector construction
    auto detector_construction = new CMSDetectorConstruction(fRunMgr);
    detector_construction->SetGDMLFile(cms_geometry_filename); 
    fRunMgr->SetDetectorConstruction( detector_construction );

    // Setup a primary generator
//    printf("hepmc_event_filename=%s\n", hepmc_event_filename.c_str());
//    if (hepmc_event_filename.empty()) {
      std::cerr<<"*** GeantRunManager: setting up a GunGenerator...\n";
      //double x = rand();
      //double y = rand();
      //double z = rand();
      //double r = sqrt(x*x+y*y+z*z);
      //fPrimaryGenerator = new GunGenerator(fConfig->fNaverage, 11, fConfig->fEmax, 0, 0, 0, x/r, y/r, z/r);
      //fPrimaryGenerator = new GunGenerator(fConfig->fNaverage, 11, fConfig->fEmax, -8, 0, 0, 1, 0, 0);

      CMSParticleGun *cmsgun = new CMSParticleGun();
      cmsgun->SetNumPrimaryPerEvt(fConfig->fNaverage);
      cmsgun->SetPrimaryName("e-");
      cmsgun->SetPrimaryEnergy(fConfig->fEmax);
      double parGunPrimaryDir[3] = {0.1, 0.9, 0.1}; // will get normalized
      cmsgun->SetPrimaryDirection(parGunPrimaryDir);
      fRunMgr->SetPrimaryGenerator(fPrimaryGenerator);
      fPrimaryGenerator = cmsgun;
//    }
    // else {
    //   //.. here for a HepMCGenerator
    //   std::cerr<<"*** GeantRunManager: setting up a HepMCGenerator...\n";
    //   fPrimaryGenerator = new HepMCGenerator(hepmc_event_filename);
    // }

    CMSApplicationTBB *cmsApp = new CMSApplicationTBB(fRunMgr, cmsgun);
    cmsApp->SetPerformanceMode(performance);
    std::cerr<<"*** GeantRunManager: setting up CMSApplicationTBB...\n";
    fRunMgr->SetUserApplication( cmsApp );

    // Start simulation for all propagators
    std::cerr<<"*** GeantRunManager: initializing...\n";
    fRunMgr->Initialize();
    fPrimaryGenerator->InitPrimaryGenerator();

    printf("==========================================================\n");
    printf("= GeantV initialized using maximum %d worker threads ====\n",
	   fRunMgr->GetNthreads());
    printf("==========================================================\n");
  }

  void 
  GeantVProducer::produce(edm::Event& iEvent) {

    // not actually needed because GeantVProducer does not depend on any other cms module
    constexpr size_t kNMax = 100;
    int sum=0;
    for(std::vector<const Getter*>::iterator it = m_getters.begin(), itEnd=m_getters.end();
        it != itEnd;
        ++it) {
      sum += iEvent.get(*it);
    }
    //std::cerr<<"GeantVProducer::produce(): m_getters.size() = "<< m_getters.size() <<" and sum="<< sum <<"\n";

    std::cerr << "GeantVProducer::produce(): *** Run GeantV simulation task ***\n";
    // if first argument is set to >1, then that number of events will be given to GeantV on each EvenSet
    // Provisionally prepare for more events in future
    size_t event_indices[kNMax];
    event_indices[0] = iEvent.index();
    RunTransportTask(1, event_indices);

    std::cerr<<"GeantVProducer <"<< label().c_str() <<"> at "<< this <<": adding to event...\n";
    iEvent.put(this,"",static_cast<int>(iEvent.index()));
    std::cerr<<"GeantVProducer <"<< label().c_str() <<"> at "<< this <<": done!\n";
  }

  /// This is the entry point for the user code to transport as a task a set of events
  bool GeantVProducer::RunTransportTask(size_t nevents, const size_t *event_index)
  {
    // First book a transport task from GeantV run manager
    //int ntotransport = BookEvents(nevents);

    GeantTaskData *td = fRunMgr->BookTransportTask();
    std::cerr<<" RunTransportTask: td= "<< td <<", nevts="<< nevents <<" toy EventID="<<event_index[0]<<"\n";
    if (!td) return false;

    // ... then create the event set
    Geant::EventSet *evset = GenerateEventSet(nevents, event_index, td);

    // ... finally invoke the GeantV transport task
    bool transported = fRunMgr->RunSimulationTask(evset, td);

    // Now we could run some post-transport task
    std::cerr<<" RunTransportTask: task "<< td->fTid <<" nevts="<< nevents
	     <<": transported="<< transported <<"\n";

    return transported;
  }

  int GeantVProducer::BookEvents(int nrequested) {
    static std::atomic_int ntotal(0);
    int nbooked = 0;
    for (int i = 0; i< nrequested; ++i) {
      int ibooked = ntotal.fetch_add(1);
      if (ibooked < fNevents) nbooked++;
    }
    return nbooked;
  }

  Geant::EventSet *GeantVProducer::GenerateEventSet(size_t nevents, const size_t *event_index, Geant::GeantTaskData *td)
  {
    using EventSet = Geant::EventSet;
    using GeantEvent = Geant::GeantEvent;
    using GeantEventInfo = Geant::GeantEventInfo;
    using GeantTrack = Geant::GeantTrack;

    EventSet *evset = new EventSet(nevents);
    LockEventGenerator();
    for (size_t i=0 ; i< nevents; ++i) {
      GeantEvent *event = new GeantEvent();
      GeantEventInfo event_info = fPrimaryGenerator->NextEvent(td);
      while (event_info.ntracks == 0) {
        std::cerr<<"Discarding empty event\n";
        event_info = fPrimaryGenerator->NextEvent(td);
      }
      event->SetEvent(event_index[i]);
      event->SetNprimaries(event_info.ntracks);
      event->SetVertex(event_info.xvert, event_info.yvert, event_info.zvert);
      for (int itr = 0; itr < event_info.ntracks; ++itr) {
        GeantTrack &track = td->GetNewTrack();
        int trackIndex = event->AddPrimary(&track);
        track.SetParticle(trackIndex);
        track.SetPrimaryParticleIndex(itr);
        fPrimaryGenerator->GetTrack(itr, track, td);
      }
      evset->AddEvent(event);
    }
    UnlockEventGenerator();
    return evset;
  }

}
REGISTER_PRODUCER(demo::GeantVProducer);
