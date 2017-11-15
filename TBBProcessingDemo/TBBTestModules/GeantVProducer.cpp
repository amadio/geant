//
//  GeantVProducer.cpp
//  DispatchProcessingDemo
//
//  Created by Chris Jones on 2016/04/25.
//

#include <vector>
#include <memory>
#include "ConfiguredProducer.h"
#include "Event.h"
#include "Waiter.h"

#include "tbb/task.h"

#include "Geant/Config.h"
#include "GeantRunManager.h"
#include "EventSet.h"
#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TTabPhysProcess.h"
#include "CMSApplicationTBB.h"
#include "CMSDetectorConstruction.h"

using namespace Geant;

namespace {

  /* Stand in for work done by GeantV via TBB based tasks*/

  class SimulateTask : public tbb::task {
    tbb::task* m_waitingTask;

  public:
    SimulateTask(tbb::task* iWaiting) :
      m_waitingTask(iWaiting) {}

    tbb::task* execute() override {
      printf("SimTask: iWaiting=%p\n", m_waitingTask);

      //IMPORTANT: decremeting the reference count is what
      // tells the framework the simulation is done with
      // this event.
      m_waitingTask->decrement_ref_count();
      return nullptr;
    }

  };
}

namespace demo {

  static int n_events = 10;
  static int n_propagators = 1;
  static int n_threads = 1;
  static int n_buffered = 5;

  class GeantVProducer : public ConfiguredProducer {
  public:
    GeantVProducer(const boost::property_tree::ptree& iConfig);
  private:
    virtual void produce(edm::Event&) override;
    //void runSimulation(tbb::task*);

    /** Functions using new GeantV interface */
    bool RunTransportTask(size_t nevents);
    int BookEvents(int nrequested);

    /** @brief Generate an event set to be processed by a single task.
	Not required as application functionality, the event reading or generation
	can in the external event loop.
    */
    Geant::EventSet* GenerateEventSet(size_t nevents, GeantTaskData *td);

    std::vector<const Getter*> m_getters;
    GeantConfig* fConfig;
    GeantRunManager* fRunMgr;
    PrimaryGenerator* fPrimaryGenerator;
  };

  GeantVProducer::GeantVProducer(const boost::property_tree::ptree& iConfig)
    : ConfiguredProducer(iConfig,kThreadSafeBetweenInstances)
    , fConfig(0)
    , fRunMgr(0)
    , fPrimaryGenerator(0)
  {
    registerProduct(demo::DataKey());
    
    for(const boost::property_tree::ptree::value_type& v: iConfig.get_child("toGet")) {
      m_getters.push_back(registerGet(v.second.get<std::string>("label"), 
                                      ""));
    }

    static int n_track_max = 500;
    static int n_learn_steps = 0;
    static int n_reuse = 100000;
    static bool monitor = false, score = false, debug = false, coprocessor = false, tbbmode = false, usev3 = true, usenuma = false;
    bool performance = true;

    //std::string cms_geometry_filename("cms2015.root");
    //std::string cms_geometry_filename("cms2018.gdml");
    std::string cms_geometry_filename("ExN03.root");
    std::string xsec_filename("xsec_FTFP_BERT.root");
    std::string fstate_filename("fstate_FTFP_BERT.root");

    //std::string hepmc_event_filename("pp14TeVminbias.root");  // sequence #stable: 608 962 569 499 476 497 429 486 465 619
    //std::string hepmc_event_filename("minbias_14TeV.root"); // sequence #stable: 81 84 93 97 87 60 106 91 92 60
    std::string hepmc_event_filename(""); // use gun generator!

    // instantiate configuration helper
    fConfig = new GeantConfig();

    fConfig->fRunMode = GeantConfig::kExternalLoop;

    fConfig->fGeomFileName = cms_geometry_filename;
    fConfig->fNtotal = n_events;
    fConfig->fNbuff = n_buffered;
    // Default value is 1. (0.1 Tesla)
    fConfig->fBmag = 40.; // 4 Tesla

    // V3 options
    fConfig->fNstackLanes = 10;
    fConfig->fNmaxBuffSpill = 128;  // New configuration parameter!!!
    fConfig->fUseV3 = usev3;

    if (tbbmode) fConfig->fRunMode = GeantConfig::kExternalLoop;
    fConfig->fUseRungeKutta = false;  // Enable use of RK integration in field for charged particles
    // prop->fEpsilonRK = 0.001;      // Revised / reduced accuracy - vs. 0.0003 default

    fConfig->fUseNuma = usenuma;
    fConfig->fNminThreshold = 5 * n_threads;
    fConfig->fNaverage = 5;

    fConfig->fUseMonitoring = monitor;
    fConfig->SetMonitored(GeantConfig::kMonQueue, false);
    fConfig->SetMonitored(GeantConfig::kMonMemory, monitor);
    fConfig->SetMonitored(GeantConfig::kMonBasketsPerVol, false);
    fConfig->SetMonitored(GeantConfig::kMonVectors, false);
    fConfig->SetMonitored(GeantConfig::kMonConcurrency, false);
    fConfig->SetMonitored(GeantConfig::kMonTracksPerEvent, false);

    // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
    // If set to 0 takes the default value of 0.01
    fConfig->fPriorityThr = 0.05;

    // Initial vector size, this is no longer an important model parameter,
    // because is gets dynamically modified to accomodate the track flow
    fConfig->fNperBasket = 16; // Initial vector size (tunable)

    // This is now the most important parameter for memory considerations
    fConfig->fMaxPerBasket = n_track_max;  // Maximum vector size (tunable)

    //if (fConfig) fConfig->fMaxRes = 0;  //??? why the if()!?  and why overwritting previous line?
    fConfig->fEmin = 3.e-6; //  [3 KeV] energy cut
    fConfig->fEmax = 0.3;   // [300MeV] used for now to select particle gun energy

    // Activate debugging using -DBUG_HUNT=ON in your cmake build
    if (debug) {
      fConfig->fUseDebug = true;
      fConfig->fDebugTrk = 1;
      //propagator->fDebugEvt = 0;
      //propagator->fDebugStp = 0;
      //propagator->fDebugRep = 10;
    }

    // Number of steps for learning phase (tunable [0, 1e6])
    // if set to 0 disable learning phase
    fConfig->fLearnSteps = n_learn_steps;
    if (performance) fConfig->fLearnSteps = 0;

    // Activate I/O
    fConfig->fFillTree = false;

    // Set threshold for tracks to be reused in the same volume
    fConfig->fNminReuse = n_reuse;

    // Activate standard scoring   
    fConfig->fUseStdScoring = true;
    if (performance) fConfig->fUseStdScoring = false;

     // Create run manager
    std::cout<<"*** GeantRunManager: instantiating with "<< n_propagators <<" propagators and "<< n_threads <<" threads.\n";
    fRunMgr = new GeantRunManager(n_propagators, n_threads, fConfig);

    // Detector construction
    fRunMgr->SetDetectorConstruction( new CMSDetectorConstruction(cms_geometry_filename.c_str(), fRunMgr) );

    // Create the tabulated physics process
    std::cout<<"*** GeantRunManager: setting physics process...\n";
    fRunMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str()));

#ifdef USE_VECGEOM_NAVIGATOR
#ifdef USE_ROOT
    fRunMgr->LoadVecGeomGeometry();
#else
    fRunMgr->LoadGeometry(cms_geometry_filename.c_str());
#endif
#endif

    // Setup a primary generator
    if (hepmc_event_filename.empty()) {
      std::cout<<"*** GeantRunManager: setting up a GunGenerator...\n";
      fPrimaryGenerator = new GunGenerator(fConfig->fNaverage, 13, fConfig->fEmax, -8, 0, 0, 1, 0, 0);
    } else {
      std::cout<<"*** GeantRunManager: setting up a HepMCGenerator...\n";
      fPrimaryGenerator = new HepMCGenerator(hepmc_event_filename);
    }
    fPrimaryGenerator->InitPrimaryGenerator();

    CMSApplicationTBB *cmsApp = new CMSApplicationTBB(fRunMgr);
    std::cout<<"*** GeantRunManager: setting up CMSApplicationTBB...\n";
    fRunMgr->SetUserApplication( cmsApp );
    if (score) {
      cmsApp->SetScoreType(CMSApplicationTBB::kScore);
    } else {
      cmsApp->SetScoreType(CMSApplicationTBB::kNoScore);
    }

    // Start simulation for all propagators
    std::cout<<"*** GeantRunManager: initializing...\n";
    fRunMgr->Initialize();

    /*
    printf("==========================================================================\n");
    printf("= GeantV run started with %d propagator(s) using %d worker threads each ====\n",
	   fRunMgr->GetNpropagators(), fRunMgr->GetNthreads());

    if (!fConfig->fUsePhysics) printf("  Physics OFF\n");
    else                       printf("  Physics ON\n");

    if (!fConfig->fUseRungeKutta) printf("  Runge-Kutta integration OFF\n");
    else                          printf("  Runge-Kutta integration ON with epsilon= %g\n", fConfig->fEpsilonRK);
    printf("==========================================================================\n");
    */
  }

  void 
  GeantVProducer::produce(edm::Event& iEvent) {

    // just pass pointer to task to be called when this event is simulated
    // send waiting task to CMSApplication for decrementing when event simulation is complete
    //CMSApplicationTBB *cmsApp = static_cast<CMSApplicationTBB*>(fRunMgr->GetUserApplication());
    //cmsApp->SetEventContinuationTask( iEvent.index(), pWaitTask );

    printf("GeantVProducer::produce(): *** Run GeantV simulation task ***\n");
    RunTransportTask(1);

    int sum=0;
    for(std::vector<const Getter*>::iterator it = m_getters.begin(), itEnd=m_getters.end();
        it != itEnd;
        ++it) {
      sum +=iEvent.get(*it);
    }
    printf("GeantVProducer::produce(): m_getters.size() = %lu and sum=%i\n", m_getters.size(), sum);

    // printf("GeantVProducer %s at %p: produce()... runSimulation(%p)\n",label().c_str(), this, pWaitTask);
    // runSimulation(pWaitTask);

    printf("GeantVProducer %s at %p: adding to event...\n",label().c_str(), this);
    iEvent.put(this,"",static_cast<int>(sum));
    printf("GeantVProducer %s at %p: done!\n",label().c_str(), this);
  }

  /* GeantV interface code goes there */
  //void
  //GeantVProducer::runSimulation(tbb::task* iWaitingTask) {
    //auto simTask = new (tbb::task::allocate_root()) SimulateTask{iWaitingTask};
    //tbb::task::spawn(*simTask);
  //}

  /// This is the entry point for the user code to transport as a task a set of events
  bool GeantVProducer::RunTransportTask(size_t nevents)
  {
    // First book a transport task from GeantV run manager
    int ntotransport = 0;
    while ((ntotransport = BookEvents(nevents))) {
      GeantTaskData *td = fRunMgr->BookTransportTask();
      if (!td) return false;

      // ... then create the event set using in this case the user application
      Geant::EventSet *evset = GenerateEventSet(ntotransport, td);

      // ... finally invoke the GeantV transport task
      bool transported = fRunMgr->RunSimulationTask(evset, td);
      // Now we could run some post-transport task
      if (!transported) return false;
    }
    return true;
  }

  int GeantVProducer::BookEvents(int nrequested) {
    static atomic_int ntotal(0);
    int nbooked = 0;
    for (int i = 0; i< nrequested; ++i) {
      int ibooked = ntotal.fetch_add(1);
      if (ibooked < n_events) nbooked++;
    }
    return nbooked;
  }

  Geant::EventSet *GeantVProducer::GenerateEventSet(size_t nevents, Geant::GeantTaskData *td)
  {
    using EventSet = Geant::EventSet;
    using GeantEvent = Geant::GeantEvent;
    using GeantEventInfo = Geant::GeantEventInfo;
    using GeantTrack = Geant::GeantTrack;

    EventSet *evset = new EventSet(nevents);
    for (size_t i=0 ; i< nevents; ++i) {
      GeantEvent *event = new GeantEvent();
      GeantEventInfo event_info = fPrimaryGenerator->NextEvent();
      while (event_info.ntracks == 0) {
	printf("Discarding empty event\n");
	event_info = fPrimaryGenerator->NextEvent();
      }
      event->SetNprimaries(event_info.ntracks);
      event->SetVertex(event_info.xvert, event_info.yvert, event_info.zvert);
      for (int itr = 0; itr < event_info.ntracks; ++itr) {
	GeantTrack &track = td->GetNewTrack();
	track.fParticle = event->AddPrimary(&track);
	track.SetPrimaryParticleIndex(itr);
	fPrimaryGenerator->GetTrack(itr, track);
      }
      evset->AddEvent(event);
    }
    return evset;
  }

}
REGISTER_PRODUCER(demo::GeantVProducer);
