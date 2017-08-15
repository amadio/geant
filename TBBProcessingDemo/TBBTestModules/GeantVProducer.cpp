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
#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TaskMgrTBB.h"
#include "TTabPhysProcess.h"
#include "CMSApplicationTBB.h"

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

  class GeantVProducer : public ConfiguredProducer {
  public:
    GeantVProducer(const boost::property_tree::ptree& iConfig);
  private:
    virtual void produce(edm::Event&) override;

    void runSimulation(tbb::task*);
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

    static int n_events = 10;
    static int n_buffered = 5;
    static int n_threads = 1;
    static int n_track_max = 64;
    static int n_learn_steps = 0;
    static int n_reuse = 4;
    static int n_propagators = 1;
    static int max_memory = 0; /* MB */
    static bool monitor = false, score = false, debug = false, coprocessor = false, tbbmode = true, performance = false;
    std::string cms_geometry_filename("cms2015.root");
    std::string xsec_filename("xsec_FTFP_BERT_G496p02_1mev.root");
    std::string fstate_filename("fstate_FTFP_BERT_G496p02_1mev.root");

    //std::string hepmc_event_filename("pp14TeVminbias.root");  // sequence #stable: 608 962 569 499 476 497 429 486 465 619
    std::string hepmc_event_filename("minbias_14TeV.root"); // sequence #stable: 81 84 93 97 87 60 106 91 92 60

    // instantiate configuration helper
    fConfig = new GeantConfig();

    fConfig->fGeomFileName = cms_geometry_filename;
    fConfig->fNtotal = n_events;
    fConfig->fNbuff = n_buffered;
    // Default value is 1. (0.1 Tesla)
    fConfig->fBmag = 40.; // 4 Tesla

    // Enable use of RK integration in field for charged particles
    fConfig->fUseRungeKutta = false;
    // prop->fEpsilonRK = 0.001;  // Revised / reduced accuracy - vs. 0.0003 default

    fConfig->fNminThreshold=5 * n_threads;
    fConfig->fUseMonitoring = monitor;
    fConfig->fNaverage = 500;

    fConfig->SetMonitored(GeantConfig::kMonQueue, monitor);
    fConfig->SetMonitored(GeantConfig::kMonMemory, monitor);
    fConfig->SetMonitored(GeantConfig::kMonBasketsPerVol, monitor);
    fConfig->SetMonitored(GeantConfig::kMonVectors, monitor);
    fConfig->SetMonitored(GeantConfig::kMonConcurrency, monitor);
    fConfig->SetMonitored(GeantConfig::kMonTracksPerEvent, monitor);
    // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
    // If set to 0 takes the default value of 0.01
    fConfig->fPriorityThr = 0.1;

    // Initial vector size, this is no longer an important model parameter,
    // because is gets dynamically modified to accomodate the track flow
    fConfig->fNperBasket = 16; // Initial vector size

    // This is now the most important parameter for memory considerations
    fConfig->fMaxPerBasket = n_track_max;

    // Maximum user memory limit [MB]
    fConfig->fMaxRes = max_memory;
    //if (fConfig) fConfig->fMaxRes = 0;  //??? why the if()!?  and why overwritting previous line?
    fConfig->fEmin = 0.001; // [1 MeV] energy cut
    fConfig->fEmax = 0.01;  // 10 MeV
    if (debug) {
      fConfig->fUseDebug = true;
      fConfig->fDebugTrk = 1;
      //propagator->fDebugEvt = 0;
      //propagator->fDebugStp = 0;
      //propagator->fDebugRep = 10;
    }
    fConfig->fUseMonitoring = monitor;

    // Set threshold for tracks to be reused in the same volume
    fConfig->fNminReuse = n_reuse;

    // Activate standard scoring   
    fConfig->fUseStdScoring = true;
    if (performance) fConfig->fUseStdScoring = false;
    fConfig->fLearnSteps = n_learn_steps;
    if (performance) fConfig->fLearnSteps = 0;

    // Activate I/O
    fConfig->fFillTree = false;
    fConfig->fTreeSizeWriteThreshold = 100000;
    // Activate old version of single thread serialization/reading
    //   fConfig->fConcurrentWrite = false;

    // Create run manager
    std::cout<<"*** GeantRunManager: instantiating with "<< n_propagators <<" propagators and "<< n_threads <<" threads.\n";
    fRunMgr = new GeantRunManager(n_propagators, n_threads, fConfig);

    // Create the tab. phys process.
    std::cout<<"*** GeantRunManager: setting physics process...\n";
    fRunMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str()));
    std::cout<<"*** GeantRunManager: setting physics process...\n";

// #ifdef USE_VECGEOM_NAVIGATOR
// #ifdef USE_ROOT
//  fRunMgr->LoadVecGeomGeometry();
// #else
//  fRunMgr->LoadGeometry(cms_geometry_filename.c_str());
// #endif
// #endif

    // setup a primary generator
    if (hepmc_event_filename.empty()) {
      std::cout<<"*** GeantRunManager: setting up a GunGenerator...\n";
      //fRunMgr->SetPrimaryGenerator( new GunGenerator(fConfig->fNaverage, 11, fConfig->fEmax, -8, 0, 0, 1, 0, 0) );
      fPrimaryGenerator = new GunGenerator(fConfig->fNaverage, 11, fConfig->fEmax, -8, 0, 0, 1, 0, 0);
    } else {
      std::cout<<"*** GeantRunManager: setting up a HepMCGenerator...\n";
      // fRunMgr->SetPrimaryGenerator( new HepMCGenerator(hepmc_event_filename) );
      fPrimaryGenerator = new HepMCGenerator(hepmc_event_filename);
    }
    fPrimaryGenerator->InitPrimaryGenerator();

    CMSApplicationTBB *CMSApp = new CMSApplicationTBB(fRunMgr);
    std::cout<<"*** GeantRunManager: setting up CMSApplicationTBB...\n";
    fRunMgr->SetUserApplication( CMSApp );
    if (score) {
      CMSApp->SetScoreType(CMSApplicationTBB::kScore);
    } else {
      CMSApp->SetScoreType(CMSApplicationTBB::kNoScore);
    }

//#ifdef GEANT_TBB
    if (tbbmode) {
      std::cout<<"*** GeantRunManager: instantiating TaskMgrTBB...\n";
      fRunMgr->SetTaskMgr( new TaskMgrTBB() );
    }
//#endif

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

    std::shared_ptr<tbb::task> waitTask{new (tbb::task::allocate_root()) tbb::empty_task{},
        [](tbb::task* iTask){tbb::task::destroy(*iTask);} };
    waitTask->set_ref_count(1+1);
    tbb::task* pWaitTask = waitTask.get();

    // just pass pointer to task to be called when this event is simulated
    // send waiting task to CMSApplication for decrementing when event simulation is complete
    CMSApplicationTBB *CMSApp = static_cast<CMSApplicationTBB*>(fRunMgr->GetUserApplication());
    printf("GVProd::produce(): iEvent.index()=%lu and pWaitTask=<%p>\n", iEvent.index(), pWaitTask);
    CMSApp->SetEventContinuationTask( iEvent.index(), pWaitTask );

    // start GeantV simulation here...  needs to be called just once
    static bool geantvStarted = false;
    if(!geantvStarted) {
      geantvStarted = true;
      printf("GeantVProducer::produce(): *** Starting GeantV simulation ***\n");
      fRunMgr->RunSimulation();
    }

    int sum=0;
    for(std::vector<const Getter*>::iterator it = m_getters.begin(), itEnd=m_getters.end();
        it != itEnd;
        ++it) {
      sum +=iEvent.get(*it);
    }
    printf("GeantVProducer::produce(): m_getters.size() = %lu and sum=%i\n", m_getters.size(), sum);

    // printf("GeantVProducer %s at %p: produce()... runSimulation(%p)\n",label().c_str(), this, pWaitTask);
    // runSimulation(pWaitTask);

    printf("GeantVProducer %s at %p: wait_for_all()...\n",label().c_str(), this);
    waitTask->wait_for_all();

    printf("GeantVProducer %s at %p: adding to event...\n",label().c_str(), this);
    iEvent.put(this,"",static_cast<int>(sum));
    printf("GeantVProducer %s at %p: done!\n",label().c_str(), this);
  }

  /* GeantV interface code goes there */
  void
  GeantVProducer::runSimulation(tbb::task* iWaitingTask) {
    //Need to decrement reference count on iWaitingTask when Event has been fully simulated
    // printf("geantv::runSimulation(%p)...\n", iWaitingTask);

    //auto simTask = new (tbb::task::allocate_root()) SimulateTask{iWaitingTask};
    //tbb::task::spawn(*simTask);
  }
  
}
REGISTER_PRODUCER(demo::GeantVProducer);
