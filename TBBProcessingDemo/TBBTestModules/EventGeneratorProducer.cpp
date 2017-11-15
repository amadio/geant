//
//  EventGeneratorProducer.cpp
//  DispatchProcessingDemo
//
//  Created by Chris Jones on 10/3/11.
//  Copyright 2011 FNAL. All rights reserved.
//

#include <vector>
#include "ConfiguredProducer.h"
#include "thread_type_from_config.h"
//#include "busy_wait_scale_factor.h"
//#include "busyWait.h"
#include "Event.h"
//#include "Waiter.h"
#include "GeantEventServer.h"


namespace demo {

  class EventGeneratorProducer : public ConfiguredProducer {
  public:
    EventGeneratorProducer(const boost::property_tree::ptree& iConfig);

  private:
    void doWait(double iSeconds);
    virtual void produce(edm::Event&);

    std::vector<const Getter*> m_getters;
    //std::string hepmc_event_filename;
  };

  EventGeneratorProducer::EventGeneratorProducer(const boost::property_tree::ptree& iConfig):
  ConfiguredProducer(iConfig,thread_type_from_config(iConfig))
  //, Waiter(iConfig)
  {
    registerProduct(demo::DataKey());
    //registerProduct(demo::DataKey("evgenKey"));

    // read primary tracks from file
    std::string filename = iConfig.get<std::string>("fileName");

    // if (m_generator) {
    //   std::cout<<"*** GeantRunManager: setting up a GunGenerator...\n";
    //   fRunMgr->SetPrimaryGenerator( new GunGenerator(fConfig->fNaverage, 11, fConfig->fEmax, -8, 0, 0, 1, 0, 0) );
    // } else {
    //   std::cout<<"*** GeantRunManager: setting up a HepMCGenerator...\n";
    //   // propagator->fPrimaryGenerator->SetEtaRange(-2.,2.);
    //   // propagator->fPrimaryGenerator->SetMomRange(0.,0.5);
    //   // propagator->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");
    //   fRunMgr->SetPrimaryGenerator( new HepMCGenerator(filename) );
    // }

    for(const boost::property_tree::ptree::value_type& v: iConfig.get_child("toGet")) {
      const std::string par{ iConfig.get<std::string>("label") };
      printf("EvtGenProducer: iConfig.get_child(""toGet"") = %s\n", par.c_str()); 
      m_getters.push_back(registerGet(v.second.get<std::string>("label"), ""));
    }
  }
  
  /*
  void
  EventGeneratorProducer::doWait(double iSeconds){
    busyWait(iSeconds*busy_wait_scale_factor);
    printf("EventGeneratorProducer::doWait()... iSeconds=%d\n", iSeconds);
  }
  */

  void
  EventGeneratorProducer::produce(edm::Event& iEvent) {
    printf("Producer %s: produce()...\n",label().c_str());


    int sum=0;
    for(std::vector<const Getter*>::iterator it = m_getters.begin(), itEnd=m_getters.end();
        it != itEnd;
        ++it) {
      sum += iEvent.get(*it);
      printf("%s got %s with value %i\n",label().c_str(), (*it)->label().c_str(), iEvent.get((*it)));
    }
    //wait(iEvent);
    iEvent.put(this,"",static_cast<int>(sum));
  }
}
REGISTER_PRODUCER(demo::EventGeneratorProducer);
